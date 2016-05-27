library(eQTLtools) # custom function library
library(icreport)
library(gtools)
library(rhdf5)

# modifying directory list to current directory for qsub script running in $TMPDIR

cis_distance <- 1000000


setwd("/zenodotus/dat01/mezeylab_scratch/jij2009/confeti_paper_results/")
output_dir <- "/zenodotus/dat01/mezeylab_scratch/jij2009/confeti_paper_results//revised_hits/"
# Getting directories which hold eQTL_results

#eqtlmerged.output.path <- paste("/zenodotus/dat01/mezeylab_scratch/jij2009/eqtl_multi_method/", eqtl.merged.directory,"/", sep = "")

dir_list <- grep("pyresult", dir(), value = TRUE)


for(d in dir_list){
    message("Directory = ", d, "\n")
    message("Reading Results from xqtl python output \n")

    file.list <- dir(d)
    cat("Loading Benjamini-Hochberg corrected p-values \n")
    pval.file <- grep("_hits[.]txt", file.list, value = TRUE)

    # change column names
    dataset_prefix <- unlist(strsplit(d,"_"))

    input_file <- grep(dataset_prefix[1], dir("/zenodotus/dat01/mezeylab_scratch/jij2009/eqtl_h5_with_Kmx"), value = TRUE)
    input_file <- paste("/zenodotus/dat01/mezeylab_scratch/jij2009/eqtl_h5_with_Kmx/", input_file, sep = "")

    message("Loading Information from ", input_file)
    gene.info <- as.data.frame(h5read(input_file, "phenotypes/col_info"))
    colnames(gene.info)[grep("id", colnames(gene.info))] <- "pheno_id"

    if(is.null(gene.info$pheno_symbol)){
        # In case the gene symbol is unknown (non-human dataset)
        gene.info$pheno_symbol <- gene.info$pheno_id
    }

    snp.info <- as.data.frame(h5read(input_file, "genotypes/col_info"))
    colnames(snp.info)[grep("id", colnames(snp.info))] <- "geno_id"

    for( i in 1:length(pval.file)){
        single_file <- pval.file[i]
        message("Processing file = ", single_file)
        file_data_info <- unlist(strsplit(single_file,"_"))

        merged_h5_file <- paste(file_data_info[1],"_" ,file_data_info[2],"_eqtldf.h5", sep = "")

        merged_h5_file <- paste(output_dir, merged_h5_file, sep = "")
        #check if merged file already exists
        message("Creating ", merged_h5_file, "\n") 

        eqtl_hits <- read.table(paste("./",d,"/",single_file, sep = ""), header = TRUE, stringsAsFactors = FALSE)

        eqtl_hits <- subset(eqtl_hits, p.bh <= 0.1)
        geno_names <- as.character(h5read(input_file, "genotypes/col_info/id"))
        eqtl_hits$genotype <- geno_names[eqtl_hits$genotype.idx]

        pheno_names <- as.character(h5read(input_file, "phenotypes/col_info/id"))
        eqtl_hits$phenotype <- pheno_names[eqtl_hits$phenotype.idx]


        H5close()

        eqtl_hits <- na.omit(eqtl_hits)

        clean.eqtl <- eQTLtools::cleaning_results(eqtl_hits, snp.info, gene.info, cis_distance)

        clean.eqtl$eqtl_id <- paste(clean.eqtl$pheno_symbol, clean.eqtl$cis_trans, clean.eqtl$geno_chr, sep ="_")

        #clean.eqtl$geno_chr <- factor(clean.eqtl$geno_chr, levels = mixedsort(levels(clean.eqtl$geno_chr)))

        h5createFile(merged_h5_file)
        h5createGroup(merged_h5_file, "eqtl_info")

        nrow_df <- nrow(clean.eqtl)
        chunk_size <- min(nrow_df, 100)
        for(info_col in names(clean.eqtl)){
        #   info_col <- "hit_id"
            # save in chunks to avoid failing to save large datasets
            if(info_col %in% c("geno_pos", "pheno_start", "pheno_end", "pval","p.bh")){
                    h5createDataset(file = merged_h5_file, dataset = paste0("eqtl_info/", info_col),
                            dims = c(length(clean.eqtl[,info_col]),1), chunk = c(chunk_size,1), level = 1,
                            storage.mode = "double")

            } else {
                    h5createDataset(file = merged_h5_file, dataset = paste0("eqtl_info/", info_col),
                            dims = c(length(clean.eqtl[,info_col]),1), chunk = c(chunk_size,1), level = 1,
                            storage.mode = "character", size = 50)

            }

            h5write(clean.eqtl[,info_col], merged_h5_file, paste0("eqtl_info/", info_col))

        }



    }

}


# Read in eQTL results





