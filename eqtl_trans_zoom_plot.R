
# script to prep the data for zoom in plots

create_chr_axis <- function(snp.info){

    chromosomes <- unique(snp.info$geno_chr)
    x.axis <- matrix(0, nrow = 3, ncol = length(chromosomes))
    snp.info$idx <- 1:nrow(snp.info)
    for(k in 1:length(chromosomes)){
        j <- as.character(chromosomes[k])

        if(is.na(j)){
            x.axis[1,k] <- max(x.axis[1,]) + 1
            x.axis[2,k] <- max(snp.info[,"idx"])
        } else {
            x.axis[1,k] <- min(subset(snp.info, geno_chr == j)[,"idx"])
            x.axis[2,k] <- max(subset(snp.info, geno_chr == j)[,"idx"])
        }

    }

    x.axis[3,] <- (x.axis[1,] + x.axis[2,]) / 2

    colnames(x.axis) <- chromosomes

    return(x.axis)
}


library(rhdf5)
library(lrgpr)
# get replicating trans hits info

prefix <- "GTExAdipose"

input_RData <- paste0(prefix, "_trans_inspect.RData")

input_file_path <- "/zenodotus/dat01/mezeylab_scratch/jij2009/confeti_R_GPCA/trans_hit_inspection/"
h5_data_path <- "/zenodotus/dat01/mezeylab_scratch/jij2009/temp_xqtl_split/"

figure_output <- "/zenodotus/dat01/mezeylab_scratch/jij2009/confeti_R_GPCA/figures/manhattan_plots/"


load(paste0(input_file_path, input_RData))

method <- "CONFETI95"

revised_method <- gsub("[0-9]*","",method)

trans_eqtl_ids <- unique(Reduce(c, lapply(rep_trans_data[[revised_method]], function(x) x$eqtl_id)))




dataset_name <- names(rep_trans_data[[revised_method]])

result_list <- list()

for(i in 1:length(dataset_name)){
    single_dataset <- dataset_name[i]
    cat("Processing dataset = ", single_dataset, "\n")
    # get input h5 for corresponding method
    input_h5 <- paste0(h5_data_path,single_dataset,"/",method,"/",single_dataset,".h5")

    phenotypes <- h5read(input_h5, "/phenotypes/matrix")
    colnames(phenotypes) <- h5read(input_h5, "/phenotypes/col_info/id")
    rownames(phenotypes) <- h5read(input_h5, "/phenotypes/row_info/id")


    genotypes <- h5read(input_h5, "/genotypes/matrix")
    colnames(genotypes) <- h5read(input_h5, "/genotypes/col_info/id")
    rownames(genotypes) <- h5read(input_h5, "/genotypes/row_info/id")

    snp.info <- as.data.frame(h5read(input_h5, "genotypes/col_info"))
    geno_pc_matrix <- h5read(input_h5, 'GPCA/matrix')
    message("Correction method = ", method, "\n")
    similarity.mx <- h5read(input_h5, paste("/K_mx/", method, sep = ""))
    K.mx <- svd(similarity.mx)

    single_dataset_df <- rep_trans_data[[revised_method]][[single_dataset]]

    x.axis <- create_chr_axis(snp.info)

    unique_chr <- unique(snp.info$geno_chr)
    chr_colors <- rep(c("darkblue","darkcyan"),length.out = length(unique_chr))
    names(chr_colors) <- unique_chr
    snp.info$color <- chr_colors[snp.info$geno_chr] 
    for(j in 1:length(trans_eqtl_ids)){
        single_trans_id <- trans_eqtl_ids[j]
        cat("Processing hit = ", single_trans_id, "\n")
        # for each trans hit
        single_dataset_single_trans <- subset(single_dataset_df, 
                                            single_dataset_df$eqtl_id == single_trans_id)


        
        pheno_chr <- single_dataset_single_trans[1,"pheno_chr"]
        geno_chr <- single_dataset_single_trans[1,"geno_chr"]
        single_genotype <- single_dataset_single_trans[1,"genotype"]
        single_phenotype <- single_dataset_single_trans[1,"phenotype"]
        single_pheno_id <- single_dataset_single_trans[1,"pheno_symbol"]

        #single_cytoband <- single_dataset_single_trans[1,"cytoband"]
        #single_cytoband_info <- chrom.bands[match(single_cytoband, chrom.bands$cyto_lookup),]
        # save pvals

  #      chrom_idx <- snp.info$geno_chr == single_cytoband_info$chrom.num
  #      within_range <- snp.info$geno_pos >= single_cytoband_info$chromStart & snp.info$geno_pos <= single_cytoband_info$chromEnd

        phenotype_subset <- phenotypes[,single_phenotype]

        n_gpcs <- 2

        pvals <- lrgprApply(phenotype_subset ~ SNP + geno_pc_matrix[,c(1:n_gpcs)],
                            features= genotypes,
                            nthreads= 1,
                            decomp = K.mx,
                            progress = FALSE)
        
        trans_geno_idx <- which(names(pvals) == single_genotype)

        png(file = paste0(figure_output, single_dataset,"_", single_trans_id, "_trans_manhattan.png"), width = 1200, height = 400)

            plot(-log10(pvals), col = as.character(snp.info$color),
                main = paste(single_dataset, single_pheno_id, pheno_chr,
                              "eQTL_on",geno_chr,sep = "_"), 
                 xlab = "chromosomes", ylab = "-log10(p-value)",xaxt='n')
            text(trans_geno_idx, -log10(pvals)[trans_geno_idx], single_genotype, cex=1, pos=4, col="red")
            axis(1, at=x.axis[3,],labels=colnames(x.axis),
                 col.axis="black", las=1)
        dev.off()
        

        result_list[[single_trans_id]][[single_dataset]] <- list("sig_pheno" =single_pheno_id,
                                                                 "sig_snp" = single_genotype,
                                                                 "pvals" = pvals)

    }


}                                                            

prefix <- gsub("[.]Rdata", "", input_RData)
save(result_list, file = paste0(input_file_path, prefix, "_results.RData"))