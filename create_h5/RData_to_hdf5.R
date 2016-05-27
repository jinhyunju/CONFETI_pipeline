# create file

#h5createFile("filename")


# create groups 

#h5createGroup("filename", "groupname")


#reading hdf5 files in R

#h5read("filename", "group")

# write objects to h5 file

#h5write(object, "filename", "group")

# inspect structure 

#h5ls("filename")


library(rhdf5)

filelist <- grep("[.]RData", dir(), value = TRUE)

transfer_Kmx <- FALSE
create_h5 <- TRUE
level1.groups <- c("phenotypes", "genotypes", "covars")

for(single.file in filelist){
    message("Loading file = ", single.file, "\n")
    load(single.file)

    head(gene.info)
    head(snp.info)
    file.name <- gsub("RData", "h5", single.file)


    if(create_h5){
        message("Creating hdf5 file = ", file.name, "\n")
        h5createFile(file.name)


        for(l1 in 1:length(level1.groups)){
            h5createGroup(file.name, level1.groups[l1])
            h5createGroup(file.name, paste(level1.groups[l1], "col_info", sep = "/"))
            h5createGroup(file.name, paste(level1.groups[l1], "row_info", sep = "/"))
        }

        h5createDataset(file.name, "genotypes/matrix", c(nrow(genotypes), ncol(genotypes)), chunk = NULL, level = 0)
        h5createDataset(file.name, "phenotypes/matrix", c(nrow(phenotypes), ncol(phenotypes)), chunk = NULL, level = 0)


        if(!is.null(covars)){
            h5write(covars, file.name, "covars/matrix") 
        }


        h5write(colnames(phenotypes), file.name, "phenotypes/col_info/id")
        h5write(rownames(phenotypes), file.name, "phenotypes/row_info/id")


        if(is.null(colnames(genotypes))){
            colnames(genotypes) <- snp.info$id
        }
        h5write(colnames(genotypes), file.name, "genotypes/col_info/id")
        h5write(rownames(genotypes), file.name, "genotypes/row_info/id")

        h5write(phenotypes, file.name, "phenotypes/matrix")
        h5write(genotypes, file.name, "genotypes/matrix")
		
        if(!is.null(colnames(covars)) & !is.null(rownames(covars))){
		        h5write(colnames(covars), file.name, "covars/col_info/id")
		        h5write(rownames(covars), file.name, "covars/row_info/id")
        }

        gene.info <- gene.info[match(colnames(phenotypes), gene.info$probe),]
        snp.info <- snp.info[match(colnames(genotypes), snp.info$id),]


        gene.chr <- as.character(gene.info$chromosome)
        gene.chr <- gsub("chr", "", gene.chr)

        h5write(as.character(gene.chr), file.name, "phenotypes/col_info/pheno_chr")
        h5write(as.numeric(gene.info$start), file.name, "phenotypes/col_info/pheno_start")
        h5write(as.numeric(gene.info$end), file.name, "phenotypes/col_info/pheno_end")

        if(!is.null(gene.info$symbol)){
            h5write(as.character(gene.info$symbol), file.name, "phenotypes/col_info/pheno_symbol")
        }
        if(!is.null(gene.info$entrez)){
            h5write(as.character(gene.info$entrez), file.name, "phenotypes/col_info/entrez")
        }

        if(is.null(snp.info$chromosome)){

            snp.chr <- as.character(snp.info[,grep("chr", colnames(snp.info), value = TRUE)])
        } else {

            snp.chr <- as.character(snp.info$chromosome)
        }
        snp.chr <- gsub("chr", "", snp.chr)

        h5write(snp.chr, file.name, "genotypes/col_info/geno_chr")
        h5write(as.numeric(snp.info$position), file.name, "genotypes/col_info/geno_pos")

        H5close()
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

#if(transfer_Kmx){




#}






