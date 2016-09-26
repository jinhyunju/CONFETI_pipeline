#library(eQTLtools)
library(icreport)
library(rhdf5)
library(data.table)

file_path <- "/home/jij2009/Genotype_wing_traits_to_Jin_20160526/"
setwd(file_path)

file_name <- "yuxin_eqtl_data.h5"
if(FALSE){
phenotypes <- fread("./phenotype_wing_traits/JBDs_mean_wing_F_2_pass_rm_4Lines.txt", header = TRUE)

genotypes <- fread("./genotype_matrix/dgrp2_maf005_g03_wing_F_132_samples_for_eqtl_20151110.txt", header = TRUE)


phenotypes <- data.frame(phenotypes)

genotypes <- data.frame(genotypes)

rownames(phenotypes) <- phenotypes[,"line.name"]
rownames(genotypes) <- genotypes[,"line_names"]

phenotypes[,"line.names"] <- NULL
genotypes$line_names <- NULL

phenotypes <- as.matrix(phenotypes)
genotypes <- as.matrix(genotypes)
}

	create_h5_file(file_name)

	h5_add_data(file_name = file_name,
	            input_matrix = phenotypes,
	            data_type = "phenotypes")

	h5_add_data(file_name = file_name,
	            input_matrix = genotypes,
	            data_type = "genotypes")

	h5_add_data(file_name = file_name,
	            input_matrix = covars,
	            data_type = "covars")

geno_covar_mx <- cov(t(genotypes))

h5write(geno_covar_mx, file_name, "K_mx/GENOCOV")

	H5close()

