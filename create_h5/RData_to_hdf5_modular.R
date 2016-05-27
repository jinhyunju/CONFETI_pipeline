library(rhdf5)
library(eQTLtools)

filelist <- grep("[.]RData", dir(), value = TRUE)

transfer_Kmx <- FALSE
create_h5 <- TRUE

for(single.file in filelist){
	
    message("Loading file = ", single.file, "\n")
    load(single.file)

    head(gene.info)
    head(snp.info)
    file.name <- gsub("RData", "h5", single.file)


    if(create_h5){
        message("Creating hdf5 file = ", file.name, "\n")
        create_eqtl_input_h5(file.name)

        if(is.null(colnames(genotypes))){
            colnames(genotypes) <- snp.info$id
        }
		
		add_geno_pheno_covar_h5(file_name = file.name, 
								phenotypes = phenotypes, 
								genotypes = genotypes, 
								covars = covars)


		gene.info <- gene.info[match(colnames(phenotypes), gene.info$probe),]
		snp.info <- snp.info[match(colnames(genotypes), snp.info$id),]
		
		add_snp_info_h5(file_name = file.name, 
						snp_info = snp.info,
						id_col = "id", 
						chr_col = "chrom",
						pos_col = "position")
		
		add_gene_info_h5(file_name = file.name,
		                 gene_info = gene.info,
						 id_col = "ensembl", 
						 chr_col = "chromosome",
						 start_col = "start",
						 end_col = "end", 
						 entrez_col = "entrez",
						 symbol_col = "symbol")
						 


    }

}


if(transfer_Kmx){
    K_mx_input_path <- "/zenodotus/dat01/mezeylab_scratch/jij2009/eqtl_h5_with_Kmx/"
    no_Kmx_files <- grep("[.]h5", dir(eqtl.input.dir), value = TRUE)
    no_Kmx_files <- grep(dataset, no_Kmx_files, value = TRUE)

    for (single_h5 in no_Kmx_files){
        
        Kmx <- h5read(paste(K_mx_input_path, single_h5, sep = ""), "/Kmx/")
        H5close()
        h5createGroup(single_h5, "/K_mx/")
        
        for(j in 1:length(Kmx)){
                list_name <- names(Kmx)[j]
                group_name <- paste("/Kmx/", list_name, sep = "")
                message("Creating group = ", group_name, " in file = ", single_h5, "\n")

                mx_to_save <- Kmx[[list_name]]
                h5createDataset(single_h5, group_name,c(nrow(mx_to_save), ncol(mx_to_save)), 
                                                        chunk = NULL, level = 0)
                h5write(mx_to_save, single_h5, group_name)


        }
    H5close()

    }



}


