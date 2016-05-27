library(eQTLtools)

setwd("/zenodotus/dat01/mezeylab_scratch/jij2009/eqtl_input_data/")

RData_files <- grep("MUTHER.+[.]RData", dir(), value = TRUE)

for(r_file in RData_files){
	
	message("Loading file = ", r_file)
	load(r_file)
	file_name <- gsub("RData", "h5", r_file)
	
	message("Creating = ", file_name)
	
	create_eqtl_input_h5(file_name)

	add_geno_pheno_covar_h5(file_name = file_name,
							phenotypes = phenotypes,
							genotypes = genotypes,
							covars = covars)
						
	add_snp_info_h5(file_name = file_name,
					snp_info = snp.info,
					id_col = "id",
					chr_col = "chromosome",
					pos_col = "position")
				
	add_gene_info_h5(file_name = file_name, 
					 gene_info = gene.info,
					 id_col = "probe",
					 chr_col = "chromosome",
					 start_col = "start",
					 end_col = "end",
					 entrez_col = "entrez",
					 symbol_col = "symbol")
	rm(phenotypes, genotypes, covars, snp.info, gene.info )
	
}
