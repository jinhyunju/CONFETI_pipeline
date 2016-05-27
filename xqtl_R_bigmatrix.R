# Master script for comparing eQTL methods
options(warn=-1)
suppressMessages(library(eQTLtools))
suppressMessages(library(icreport))
suppressMessages(require(lrgpr))
suppressMessages(require(rhdf5))
require('parallel')
require('bigmemory')
options(warn=0)

as.character.formula <- function(x){
   Reduce( paste, deparse(x) )  
}

create_idx <- function(n_geno, total_batches){
  div_interval = round(seq(1,n_geno+1, length.out=total_batches+1))
  batch_idx <- cbind(div_interval[-(total_batches+1)], div_interval[-1]-1)
  return(batch_idx)
}

# method is going to be one of the following
# LINEAR, PCA, PEER, PANAMA, PANAMAGS, ICALMM, ICAFIXED, CONFETI
# PANAMAGS is the version of PANAMA with genetic similarity considered

################################################################################
#################         Getting command line inputs          #################
################################################################################

if ( length(commandArgs(TRUE)) < 1 & !exists('input.file')) {
	  stop('Please specify a data file on the command line')
} else if ( !exists('input.file') ) {
	  input.file <- as.character(commandArgs(TRUE)[1])
}


if ( length(commandArgs(TRUE)) < 2 & !exists('MC.CORES')) {
		#  stop('Please specify number of cores on the command line or set the MC.CORES variable')
		MC.CORES <- detectCores()
} else if ( !exists('MC.CORES') ) {
		MC.CORES <- as.integer(commandArgs(TRUE)[2])
}

if ( length(commandArgs(TRUE)) < 3 & !exists('N_BATCH')) {
		# Use default
		N_BATCH <- 20
} else if ( !exists('N_BATCH') ) {
		N_BATCH <- as.integer(commandArgs(TRUE)[3])
}

if ( length(commandArgs(TRUE)) < 4 & !exists('method')) {
 	 stop('Method not defined')
} else {
	  method <- as.character(commandArgs(TRUE)[4])
}


# A meal that costs a lot on a special day
PREFIX <- gsub("[.]h5", "", input.file)
PREFIX <- paste(PREFIX, "_" , method,"r", sep = "")

message(paste("Loading Genotypes and Phenotypes :",input.file))

phenotypes <- h5read(input.file, "/phenotypes/matrix")
colnames(phenotypes) <- h5read(input.file, "/phenotypes/col_info/id")
rownames(phenotypes) <- h5read(input.file, "/phenotypes/row_info/id")

genotypes <- h5read(input.file, "/genotypes/matrix")
colnames(genotypes) <- h5read(input.file, "/genotypes/col_info/id")
rownames(genotypes) <- h5read(input.file, "/genotypes/row_info/id")

# check if covars exists

h5_file <- H5Fopen(input.file)

if(H5Lexists(h5_file, "covars/matrix")){
		covars <- h5read(input.file, "/covars/matrix")
} else {
		covars <- NULL
}

H5close()

# A few checks (not exhaustive!) to see if the data is in the expected format
if ( ! exists('genotypes') ) {stop('Missing matrix: genotypes\n')}
if ( class(genotypes) != 'matrix' ) {stop('Not a matrix: genotypes\n')}
if ( mode(genotypes) != 'numeric' & mode(genotypes) != 'integer' ) {stop('Not numeric: genotypes\n')}

if ( ! exists('covars') ) {stop('Missing matrix: covars\n')}
if ( ! is.null(covars) & mode(covars) != 'numeric' & mode(covars) != 'integer') {stop('Not numeric: covars\n')}

if ( ! exists('phenotypes') ) {stop('Missing matrix: phenotypes\n')}
if ( class(phenotypes) != 'matrix' ) {stop('Not a matrix: phenotypes\n')}
if ( mode(phenotypes) != 'numeric' ) {stop('Not numeric: phenotypes\n')}
if ( nrow(genotypes) != nrow(phenotypes) ) {stop('Number of rows in genotypes & phenotypes does not match\n')}

# get matrix dimensions
n <- nrow(genotypes)
n.genotypes <- ncol(genotypes)
n.phenotypes <- ncol(phenotypes)
pheno.names <- colnames(phenotypes)
geno.names <- colnames(genotypes)
################################################################################
#################         Calculate covariate matrix           #################
################################################################################

message("Correction method = ", method, "\n")
similarity.mx <- h5read(input.file, paste("/K_mx/", method, sep = ""))
K.mx <- svd(similarity.mx)
covars <- NULL

################################################################################
#################             Actual Regression Part           #################
################################################################################

#################               Fixed Corrections             #################


# Set z to the number of columns in the design matrix (intercept + genotype + covars)
if ( is.null(covars) ) {
  # No covariates
  z <- 2
} else if ( is.null(dim(covars)) & length(covars) == n ) {
  # One column/covariate
  z <- 3
} else {
  z <- 2 + ncol(covars)
}

cat('n =', n, '\nk =', z, '(intercept + genotype +', z-2, 'covariate(s))\n')
cat('n.genotypes =', n.genotypes, '\n')
cat('n.phenotypes =', n.phenotypes, '\n')

back.file <- paste(PREFIX,'_pvals.bin', sep = "")
desc.file <- paste(PREFIX,'_pvals.desc', sep = "")

# If we previously started this analysis, the big.matrix files will
# still be there and we will use those (files can be manually removed if required)
if ( file.exists(back.file) ) {
  pvals.bm <- attach.big.matrix(desc.file)
  if ( nrow(pvals.bm) != n.genotypes | ncol(pvals.bm) != n.phenotypes ) {
    stop('Found existing big.matrix files but the dimensions do not match!')
  }
} else {
  pvals.bm <- filebacked.big.matrix(nrow=n.genotypes,
                                    ncol=n.phenotypes,
                                    type='double',
                                    backingfile=back.file,
                                    descriptorfile=desc.file,
                                    dimnames=list(geno.names, pheno.names))
}



# Based on the number of genotypes and batch size we assign blocks of genotypes
# The last job will be smaller (if necessary)
cat('N_BATCH =', N_BATCH, '\n')

geno_batch_idx <- create_idx(n.genotypes, N_BATCH)

message("Entering Linear Mixed Model correction")
message("Method = ", method)

for( b in 1:N_BATCH ){
		start_idx <- geno_batch_idx[b,1]
		end_idx <- geno_batch_idx[b,2]

        message('\nStarting batch ', b, ' of ', N_BATCH, ' : ', start_idx, ' to ', end_idx, '...\n')

		input_geno <- genotypes[,c(start_idx:end_idx)]


		pvals_list <- parallel::mclapply(1:n.phenotypes, function(x) 
										lrgprApply(phenotypes[,x] ~ SNP,
												    features=input_geno,
												    nthreads= 1,
												    decomp = K.mx,
												    progress = FALSE), 
													mc.cores = MC.CORES)

		pvals <- matrix(unlist(pvals_list), ncol = n.phenotypes, byrow = FALSE)

		pvals.sub <- sub.big.matrix(pvals.bm,
		                            firstRow=start_idx,
		                            lastRow=end_idx)
		pvals.sub[, ] <- pvals
		rm(pvals)
		gc()

}



