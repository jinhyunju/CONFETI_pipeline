# Individual dataset plotting 

# lamda plot generation script
library(rhdf5)
library(plyr)
library(ggplot2)
library(reshape2)
library(icreport)

read_in_h5_eqtl <- function(input_h5){

    h5_file <- H5Fopen(input_h5)
    updated_file <- H5Lexists(h5_file, "/eqtl_info")
    if(updated_file){
        temp.df <- as.data.frame(h5read(input_h5, "/eqtl_info"), stringsAsFactors = FALSE)
    } else {
        temp.df <- h5read(input_h5, "/eqtl_df")
    }
    
    # in case you are reading in bonferroni results
    if(is.null(temp.df$p.bh)){ 
            temp.df$p.bh <- temp.df$pval 
    }
    #temp.df <- na.omit(temp.df)
    # fix for hit ids that have a "_" in genotype or phenotype ids
    fixed_hit <- regmatches(temp.df$hit_id, regexpr("_", temp.df$hit_id), invert = TRUE)
    temp.df$fixed_hit_id <- sapply(fixed_hit, function(x) paste(x[1], x[2], sep = "::"))
    temp.df$info_id <- with(temp.df, paste(fixed_hit_id,geno_chr, geno_pos, pheno_chr, pheno_start, pheno_symbol,cis_trans, sep = "::"))
    #temp.df <- temp.df[!duplicated(temp.df$eqtl_id),]
    #temp.df <- dplyr::select(temp.df, eqtl_id, pval, p.bh)
    H5close()
    return(temp.df)    
}

cis_trans_counter <- function(list.entry){
  if(is.null(list.entry)){
    cis_trans_result <- c(0,0)

  } else {
    cis_trans <- sapply(strsplit(list.entry, "_"), function(x) x[2])
    cis_trans_result <- c(sum(cis_trans == "cis"), sum(cis_trans == "trans"))
  }
   names(cis_trans_result) <- c("cis","trans")
  return(cis_trans_result)
}

eqtl_hits_fdr_plots <- function(eqtl.merged.list, 
                                   max.fdr = 0.1, 
                                   figure.path, 
                                   data_prefix, 
                                   factor_sort = NULL){

    today_date <- gsub("-", "_", Sys.Date())
    plot_prefix <- paste(data_prefix, today_date, sep = "_")

    max.pval.data <- sapply(eqtl.merged.list, function(a) max(range(a$p.bh), na.rm = TRUE))
    max.pval.data <- max(unlist(max.pval.data), na.rm = TRUE)

    if(max.pval.data < max.fdr){

        max.fdr <- max.pval.data

    }

    FDR_range <- seq(0, max.fdr, length.out = 100)

    cis_trans_list <- list()


    for( c in 1:length(FDR_range)) {

            significant_unique <- lapply(eqtl.merged.list, function(a) unique(subset(a, a$p.bh < FDR_range[c])$eqtl_id)) # don't forget to change it to eqtl_id
            ct_count_df <- data.frame(melt(sapply(significant_unique, cis_trans_counter)))
            ct_count_df$FDR <- FDR_range[c]
            colnames(ct_count_df) <- c("cis_trans", "method", "count", "FDR")
            cis_trans_list[[c]] <- ct_count_df
    }

    cis_trans_df <- Reduce(rbind, lapply(cis_trans_list, function(x) x))

    if(!is.null(factor_sort)){
        cis_trans_df$method <- factor(cis_trans_df$method, levels = factor_sort)
    }

    cis_trans_rep_plot <- ggplot(cis_trans_df, aes(x = FDR, y = count, col = method, group = method)) + facet_grid(cis_trans~.,scales = "free_y") +
                  geom_line(size = 1)  + scale_color_manual(values = custom_pal)
                  
    cis_trans_rep_plot <- ggplot_add_theme(cis_trans_rep_plot)

    pdf(file = paste(figure.path, plot_prefix, "_eqtl_hits_per_FDR_plot.pdf", sep = ""), height = 8, width = 5)
        cis_trans_rep_plot <- cis_trans_rep_plot+ theme(legend.text = element_text(size = 10), 
                                                        legend.key.size = unit(1, "cm"),
                                                        strip.background = element_rect(colour = "white", 
                                                                                        fill = "white",
                                                                                        size = 0.8))
        print(cis_trans_rep_plot)
    dev.off()


}
sig_hits_map <- function(hits.info, prefix = "test", output_path = "./", organism = "human",
                          chr_info_file = '/zenodotus/dat01/mezeylab_scratch/jij2009/General_information/chromInfo.txt.gz'){

    # input file was downloaded from the following url
    # http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/chromInfo.txt.gz
    chrom.lengths <- read.table(chr_info_file, col.names=c('chromosome', 'size', 'file'))
    use.chromosomes <-  paste0('chr', c(1:22, 'X', 'Y', 'M'))
    chrom.lengths <- chrom.lengths[match(use.chromosomes, chrom.lengths$chromosome), "size"]

    if(organism == "yeast"){
        chr_info_file = '/zenodotus/dat01/mezeylab_scratch/jij2009/General_information/yeast_chr_length.txt'
        chrom_df <- read.table(chr_info_file, header = FALSE, 
                               col.names = c('chromosome', 'id', 'size'), stringsAsFactors = FALSE)
        chrom.lengths <- as.numeric(chrom_df$size)
        snp.chroms <- chrom_df$chromosome
        gene.chroms <- snp.chroms
    }
    snp.chroms <- gsub("chr", "", use.chromosomes)
    gene.chroms <- snp.chroms

    snp.lengths <- as.numeric(chrom.lengths)
    gene.lengths <- as.numeric(chrom.lengths)
    snp.mins <- gene.mins <- 0


    snp.offsets <- (cumsum(c(0, snp.lengths))[-length(snp.chroms)]-snp.mins)[1:length(snp.chroms)]
    gene.offsets <- (cumsum(c(0, gene.lengths))[-length(gene.chroms)]-gene.mins)[1:length(gene.chroms)]

    names(snp.offsets) <- snp.chroms
    names(gene.offsets) <- gene.chroms

    max.snp <- sum(snp.lengths)
    max.gene <- sum(gene.lengths)

    snp.ticks <- c(cumsum(c(0, snp.lengths)))
    gene.ticks <- c(cumsum(c(0, gene.lengths)))

    snp.labs <- (snp.ticks[-1] + snp.ticks[-length(snp.ticks)])/2
    gene.labs <- (gene.ticks[-1] + gene.ticks[-length(gene.ticks)])/2

    rbPal <- colorRampPalette(c('red','black'))
    # add minimum representative value to pvalues to avoid 0 -> INF conversion
    hits.info$col <- ifelse(hits.info$cis_trans == "cis", "darkgreen", "blue")
    #rbPal(2)[as.numeric(cut(c(1:nrow(hits.info$cis_trans)),breaks = 10))]
    hits.info$shape <- 20

#    message("Creating Significant Hits Map \n")
#    message("File prefix = ", prefix, "\n")
#    message("Output directory = ", output_path, "\n")
    # Set up plot
    pdf(paste(output_path,prefix, "_replicating_hits_map.pdf", sep = ""))
        old.par <- par(mar=c(3, 3, 2, 0), mgp=c(1.75, 0.5, 0))

        plot(0, type='n', xlab='SNP chromosome', ylab='Gene chromosome', xaxt='n', yaxt='n',
           xlim=c(0, max.snp), ylim=c(0, max.gene),
           bty='n', main = prefix)
        # Plot axes
        axis(1, snp.ticks, labels=FALSE)
        axis(2, gene.ticks, labels=FALSE)
        axis(1, snp.labs, snp.chroms[1:length(snp.labs)], lwd=0)
        axis(2, gene.labs, gene.chroms[1:length(gene.labs)], lwd=0)
        # Grid lines
        abline(h=gene.ticks, col='#dddddd')
        abline(v=snp.ticks, col='#dddddd')

        # Plot hits
        with(hits.info, 
            points(snp.offsets[geno_chr]+geno_pos, 
                gene.offsets[pheno_chr]+pheno_start, col =  col,pch = shape, cex = 0.6))
        par(old.par)
    dev.off()


}

get_lambda <- function(pval_median){
    lambda <- qchisq(1 - pval_median,1) / qchisq(0.5, 1) 
    return(lambda)

}
hexconvert <- function(colorname){
    hex_code <- do.call(rgb, as.list(col2rgb(colorname)/255))
    return(hex_code)
}

# Function to create custom palette for plotting
confeti_palette <- function(){
    plot_pal <- RColorBrewer::brewer.pal(n = 8, name = "Set1")[1:5]

    plot_pal <- c(plot_pal, hexconvert("gold"),hexconvert("skyblue"), hexconvert("deeppink3"))
    return(plot_pal)
}

custom_pal <- confeti_palette()
#input_path <- "/zenodotus/dat01/mezeylab_scratch/jij2009/confeti_paper_results/POP_corrected/pval_medians/"

dataset_name <- as.character(commandArgs(TRUE))[1]

master_path <- "/zenodotus/dat01/mezeylab_scratch/jij2009/confeti_R_GPCA/"

if(length(commandArgs(TRUE)) < 1){

  dataset_name <- "GTExSkinLeg"

  figure_output <- paste0(master_path, "figures/Indv_datasets/",dataset_name,"/")
  dir.create(figure_output, showWarnings = FALSE)
  pval_median_file_path <- paste0(master_path,dataset_name, "/pval_medians")

  eqtl_info_input <- paste0(master_path, dataset_name, "/eqtl_hits")

}

message("Generating genomic inflation factor plots\n")
message("Figures generated in = ", figure_output, "\n")

method_factors = c("CONFETI", "CONPANA", "PANAMA", "PCALMM", "ICE", "LINEAR")
#method_factors <- c("PARTHYBRIDpy", "PARTICApy", "PANAMApy", "ICEpy", "PEER", "LINEAR")
pval_median_files <- list.files(pval_median_file_path, full.names= TRUE, pattern = "median.h5")
file_name_only <- list.files(pval_median_file_path, full.names= FALSE, pattern = "median.h5")

pval_median_list <- list()
unique.methods <- unique(sapply(strsplit(file_name_only,"_"), function(x) x[2]))

for(i in 1:length(unique.methods)){
    revised_method <- gsub("Rgpca","",gsub("py", "", gsub("[0-9]","", unique.methods[i])))
    revised_method <- gsub("gpca", "", revised_method)
    if(revised_method %in% method_factors){
        message("Reading results for method = ", unique.methods[i])

        method.files <- grep(paste(unique.methods[i],"_",sep=""), pval_median_files, value = TRUE)

        pval_median_list[[revised_method]] <- h5read(method.files, "median_pvals")

    }
}


lambda_list <- lapply(pval_median_list, get_lambda)
lambda_methods <- as.data.frame(lambda_list)
lambda_df <- melt(lambda_methods)

colnames(lambda_df) <- c("method", "lambda")

if(!is.null(method_factors)){
    lambda_df$method <- factor(lambda_df$method, levels = method_factors )  
}

lambda_df$lambda_diff <- with(lambda_df, as.numeric(lambda) - 1)

pdf(paste0(figure_output,dataset_name, "_lambda_boxplot.pdf"), width = 6, height = 4)

    p <- ggplot(lambda_df, aes(x = method, y = lambda_diff, group = method, fill = method)) + geom_boxplot(outlier.shape = 16)  + theme_bw()
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none', plot.title = element_text(face = "bold")) + ylim(c(-1,1))+ scale_fill_manual(values = custom_pal)
    p <- p + geom_hline(aes(yintercept=0)) + ggtitle(paste0(dataset_name, "_Genomic_Inflation_Factor"))
    #p <- p + scale_fill_manual(values = custom_pal)
    print(p)
dev.off()


#### Generating hits / FDR plot

eqtl_info_files <- list.files(eqtl_info_input, pattern = ".h5", full.names = TRUE)

eqtl_hit_list <- list()

for(i in 1:length(unique.methods)){
    revised_method <- gsub("Rgpca","",gsub("py", "", gsub("[0-9]","", unique.methods[i])))
    revised_method <- gsub("gpca", "", revised_method)
    if(revised_method %in% method_factors){
        message("Reading results for method = ", unique.methods[i])

        single_info_file <- grep(unique.methods[i], eqtl_info_files, value = TRUE)

        eqtl_hit_list[[revised_method]] <- read_in_h5_eqtl(single_info_file)
        H5close()
    }
}

eqtl_hits_fdr_plots(eqtl_hit_list, 
                    max.fdr = 0.1, 
                    figure.path = figure_output, 
                    data_prefix = dataset_name, 
                    factor_sort = NULL)


for( m in names(eqtl_hit_list)){

    sig_hits_map(eqtl_hit_list[[m]], 
              prefix = paste0(dataset_name,"_",m), 
              organism="human", output_path = figure_output)

}
