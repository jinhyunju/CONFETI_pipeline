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
      message("Input File = ", input.file, "\n")
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
PREFIX <- paste(PREFIX, "_" , method, "Rgpca", sep = "")

message(paste("Loading Genotypes and Phenotypes :",input.file))

phenotypes <- h5read(input.file, "/phenotypes/matrix")
colnames(phenotypes) <- h5read(input.file, "/phenotypes/col_info/id")
rownames(phenotypes) <- h5read(input.file, "/phenotypes/row_info/id")

genotypes <- h5read(input.file, "/genotypes/matrix")
colnames(genotypes) <- h5read(input.file, "/genotypes/col_info/id")
rownames(genotypes) <- h5read(input.file, "/genotypes/row_info/id")

geno_pc_matrix <- h5read(input.file, 'GPCA/matrix')


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

pval_h5_file <- paste(PREFIX,'_pvals.h5', sep = "")


h5createFile(pval_h5_file)
h5createDataset(pval_h5_file, "/pvals", c(n.genotypes, n.phenotypes), chunk = c(100,100), level = 0)

# Based on the number of genotypes and batch size we assign blocks of genotypes
# The last job will be smaller (if necessary)
cat('N_BATCH =', N_BATCH, '\n')

geno_batch_idx <- create_idx(n.genotypes, N_BATCH)

message("Entering Linear Mixed Model correction")
message("Method = ", method)
if(method == "LINEAR"){

	for( b in 1:N_BATCH ){
		start_idx <- geno_batch_idx[b,1]
		end_idx <- geno_batch_idx[b,2]

        message('\nStarting batch ', b, ' of ', N_BATCH, ' : ', start_idx, ' to ', end_idx, '...\n')

		input_geno <- genotypes[,c(start_idx:end_idx)]

		pvals_list <- glmApply(phenotypes ~ SNP + geno_pc_matrix[,c(1:2)],
                      features=input_geno,
                      nthreads=MC.CORES)$pValues

		pvals <- matrix(unlist(pvals_list), ncol = n.phenotypes, byrow = FALSE)

        h5write(pvals, pval_h5_file, "/pvals", index = list(c(start_idx:end_idx), NULL))

		rm(pvals)
		gc()

	}


} else {


	for( p in 1:n.phenotypes ){
		cat(".")
		pvals <- vector(length=n.genotypes)
		pvals <- lrgprApply(phenotypes[,p] ~ SNP + geno_pc_matrix[,c(1:2)],
							features=genotypes,
							nthreads= MC.CORES,
							decomp = K.mx,
							progress = FALSE)

		h5write(pvals, pval_h5_file, "/pvals", index = list(NULL,p))

		rm(pvals, pvals_list)
		gc()

	}	
}






