library(eQTLtools)
load("GTExEsophagusMuscle.RData")
file_name <- "GTExEsophagusMuscle.h5"

create_eqtl_input_h5(file_name)

add_geno_pheno_covar_h5(file_name = file_name,
						phenotypes = phenotypes,
						genotypes = genotypes,
						covars = covars)
						
add_snp_info_h5(file_name = file_name,
				snp_info = snp.info,
				id_col = "id",
				chr_col = "chrom",
				pos_col = "position")
				
add_gene_info_h5(file_name = file_name, 
				 gene_info = gene.info,
				 id_col = "ensembl",
				 chr_col = "chromosome",
				 start_col = "start",
				 end_col = "end",
				 entrez_col = "entrez",
				 symbol_col = "symbol")