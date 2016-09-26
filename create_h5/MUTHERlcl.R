library(eQTLtools)
library(icreport)
setwd("/zenodotus/dat01/mezeylab_scratch/jij2009/updated_input_files/")


RData_file <- "MUTHERlcl.RData"
#RData_file <- "GTExAdiposeSubq.RData"
load(RData_file)
file_name <- gsub("RData", "h5", RData_file)

file_name <- paste0("/zenodotus/dat01/mezeylab_scratch/jij2009/eqtl_input_h5/",file_name)

message("Creating h5 File ", file_name)


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
                snp_chr = gsub("chr","",as.character(snp.info$chrom)),
                snp_pos = snp.info$position)

h5_add_pheno_info(file_name = file_name,
                    pheno_id = gene.info$probe,
                    pheno_chr = as.character(gene.info$chromosome),
                    pheno_start = gene.info$start,
                    pheno_end = gene.info$end,
                    pheno_entrez = gene.info$entrez,
                    pheno_symbol = gene.info$symbol)

H5close()