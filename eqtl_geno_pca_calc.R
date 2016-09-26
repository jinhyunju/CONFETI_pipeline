options(warn=-1)
suppressMessages(library(eQTLtools))
suppressMessages(library(icreport, lib.loc = "/home/jij2009/R_packages/"))
suppressMessages(require(lrgpr))
suppressMessages(require(rhdf5))
options(warn=0)

as.character.formula <- function(x){
   Reduce( paste, deparse(x) )  
}


input_file <- as.character(commandArgs(TRUE)[1])

message("Processing file = ", input_file, "\n")

geno <- h5read(input_file, "genotypes/matrix")
colnames(geno) <- h5read(input_file, "genotypes/col_info/id")
rownames(geno) <- h5read(input_file, "genotypes/row_info/id")


message("Creating Genotype PCA matrix \n")
geno_pca <- prcomp(geno)

var_pct <- (geno_pca$sdev^2) / sum(geno_pca$sdev^2)
n_pc <- which(cumsum(var_pct) > 0.25)[1] 

if(n_pc < 10){
    n_pc <- min(10, dim(geno))
}
geno_x <- geno_pca$x[,1:n_pc]
geno_weight <- diag(geno_pca$sdev[1:n_pc])

geno_K <- geno_x %*% geno_weight %*% geno_weight %*% t(geno_x)

h5createGroup(input_file, "GPCA")
h5write(geno_pca$x[,1:n_pc], input_file, "GPCA/matrix")
h5write(var_pct, input_file, "GPCA/var_percent")
h5write(geno_K, input_file, "GPCA/GPCA_K")
