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

if ( length(commandArgs(TRUE)) < 3 & !exists('n_batch')) {
		# Use default
		stop('Batch number not defined')
} else if ( !exists('n_batch') ) {
		n_batch <- as.integer(commandArgs(TRUE)[3])
}

if ( length(commandArgs(TRUE)) < 4 & !exists('method')) {
 	 stop('Method not defined')
} else {
	  method <- as.character(commandArgs(TRUE)[4])
}

if ( length(commandArgs(TRUE)) < 5 & !exists('n_gpcs')) {
 	 stop('Number of Genotype PCs not defined')
} else {
	  n_gpcs <- as.numeric(commandArgs(TRUE)[5])
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


batch_idx <- h5read(input.file, "batch_idx")
current_batch_idx <- batch_idx[,n_batch]

start_idx <- current_batch_idx[1] + 1
end_idx <- current_batch_idx[2]

genotypes <- genotypes[,c(start_idx:end_idx)]

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

cat('n =', n, '\nk =', z, '(intercept + genotype +', n_gpcs, 'genotype PC(s))\n')
cat('n.genotypes =', n.genotypes, '\n')
cat('n.phenotypes =', n.phenotypes, '\n')

pval_h5_file <- paste0(PREFIX,'_pvals_',sprintf("%03d",n_batch),'.h5')


h5createFile(pval_h5_file)
h5createDataset(pval_h5_file, "/pvals", c(n.genotypes, n.phenotypes), chunk = c(100,100), level = 0)
h5write(current_batch_idx, pval_h5_file, "/batch_idx")
# Based on the number of genotypes and batch size we assign blocks of genotypes
# The last job will be smaller (if necessary)
cat('Batch Number =', n_batch, '\n')

message("Entering Linear Mixed Model correction")
message("Method = ", method)
if(method == "LINEAR"){

        message('\nStarting Linear Regression \n')
		pvals_list <- glmApply(phenotypes ~ SNP + geno_pc_matrix[,c(1:n_gpcs)],
                      features=genotypes,
                      nthreads=MC.CORES)$pValues

		pvals <- matrix(unlist(pvals_list), ncol = n.phenotypes, byrow = FALSE)

        h5write(pvals, pval_h5_file, "/pvals")

		rm(pvals)
		gc()


} else {

	pvals_list <- parallel::mclapply(1:n.phenotypes, function(x) 
									lrgprApply(phenotypes[,x] ~ SNP + geno_pc_matrix[,c(1:n_gpcs)],
												features= genotypes,
												nthreads= 1,
												decomp = K.mx,
												progress = FALSE), 
												mc.cores = MC.CORES)

	pvals <- matrix(unlist(pvals_list), ncol = n.phenotypes, byrow = FALSE)

	h5write(pvals, pval_h5_file, "/pvals")

	rm(pvals)
	gc()

	if(FALSE){
		for( p in 1:n.phenotypes ){
			# cat(".")
			pvals <- vector(length=n.genotypes)
			pvals <- lrgprApply(phenotypes[,p] ~ SNP + geno_pc_matrix[,c(1:n_gpcs)],
								features=genotypes,
								nthreads= MC.CORES,
								decomp = K.mx,
								progress = FALSE)

			h5write(pvals, pval_h5_file, "/pvals", index = list(NULL,p))

			rm(pvals, pvals_list)
			gc()

		}	


	}

}






