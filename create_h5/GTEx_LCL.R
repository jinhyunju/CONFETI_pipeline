library(icreport)
library(rhdf5)

file_path <- "/zenodotus/dat01/mezeylab_scratch/jij2009/updated_input_files"
output_path <- "/zenodotus/dat01/mezeylab_scratch/jij2009/eqtl_input_h5/"

setwd(file_path)


RData_file <- "GTExLCL.RData"
load(RData_file)
file_name_h5 <- gsub("RData", "h5", RData_file)
file_name <- paste0(output_path, file_name_h5)
colnames(genotypes) <- as.character(snp.info$id)
rownames(phenotypes) <- gsub("[.]", "", rownames(phenotypes))
rownames(genotypes) <- gsub("[.]", "", rownames(genotypes))
# first sort the genes by the chromosome

#	gene.info = gene.info[order(gene.info$chrom.num, gene.info$start, decreasing = FALSE),]

#	phenotypes = phenotypes[,as.character(gene.info$probe)]

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


h5_add_snp_info(file_name = file_name,
				snp_id = snp.info$id,
				snp_chr = as.character(snp.info$chrom),
				snp_pos = snp.info$position)

h5_add_pheno_info(file_name = file_name,
					pheno_id = gene.info$ensembl,
					pheno_chr = as.character(gene.info$chromosome),
					pheno_start = gene.info$start,
					pheno_end = gene.info$end,
					pheno_entrez = gene.info$entrez,
					pheno_symbol = gene.info$symbol)

H5close()
