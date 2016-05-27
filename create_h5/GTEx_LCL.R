library(eQTLtools)

RData_file <- "GTExlcl.RData"
#RData_file <- "GTExAdiposeSubq.RData"
load(RData_file)
file_name <- gsub("RData", "h5", RData_file)
colnames(genotypes) <- as.character(snp.info$id)
message("Creating h5 File ", file_name)
create_eqtl_input_h5(file_name)

message("Adding Genotype and phenotype information")
add_geno_pheno_covar_h5(file_name = file_name,
						phenotypes = phenotypes,
						genotypes = genotypes,
						covars = covars)
						
message("Adding SNP information")
add_snp_info_h5(file_name = file_name,
				snp_info = snp.info,
				id_col = "id",
				chr_col = "chrom",
				pos_col = "position")
				
				
message("Adding Gene information")
add_gene_info_h5(file_name = file_name, 
				 gene_info = gene.info,
				 id_col = "ensembl",
				 chr_col = "chromosome",
				 start_col = "start",
				 end_col = "end",
				 entrez_col = "entrez",
				 symbol_col = "symbol")