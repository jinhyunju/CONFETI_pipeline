suppressPackageStartupMessages(library(ggplot2))

suppressPackageStartupMessages(library(icreport))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(rhdf5))

# Function to convert colornames to hexcode 
hexconvert <- function(colorname){
    hex_code <- do.call(rgb, as.list(col2rgb(colorname)/255))
    return(hex_code)
}

# Function to create custom palette for plotting
confeti_palette <- function(){
    plot_pal <- brewer.pal(n = 8, name = "Set1")[1:5]

    plot_pal <- c(plot_pal, hexconvert("gold"),hexconvert("skyblue"), hexconvert("deeppink3"))
    return(plot_pal)
}

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
eqtl_fdr_cis_trans_plots <- function(eqtl.merged.list, 
	                                   max.fdr = 0.1, 
	                                   figure.path, 
	                                   factor_sort = NULL){


    max.pval.data <- sapply(eqtl.merged.list, function(a) sapply(a, function(b) max(range(b$p.bh), na.rm = TRUE)))
    max.pval.data <- max(unlist(max.pval.data), na.rm = TRUE)

    if(max.pval.data < max.fdr){

        max.fdr <- max.pval.data

    }

    FDR_range <- seq(0, max.fdr, length.out = 100)

    cis_trans_list <- list()

	for( l in 1:length(eqtl.merged.list)){

    	for( c in 1:length(FDR_range)) {

            	significant_unique <- lapply(eqtl.merged.list[[l]], function(a) unique(subset(a, a$p.bh < FDR_range[c])$eqtl_id) ) # don't forget to change it to eqtl_id
            	ct_count_df <- data.frame(melt(sapply(significant_unique, cis_trans_counter)))
 			
	 			
	            ct_count_df$FDR <- FDR_range[c]
	            colnames(ct_count_df) <- c("cis_trans", "method", "count", "FDR")
	            cis_trans_list[[c]] <- ct_count_df
 		}

 		today_date <- gsub("-", "_", Sys.Date())
 		data_prefix <- names(eqtl.merged.list)[l]
    	plot_prefix <- paste(data_prefix, today_date, sep = "_")

 		cis_trans_df <- Reduce(rbind, lapply(cis_trans_list, function(x) x))

	    plot_df <- dcast(cis_trans_df, method + FDR ~ cis_trans, value.var = "count")

	    plot_df$total_rep <- with(plot_df, cis + trans)

	
	#    plot_df$method <- gsub("py", "", as.character(plot_df$method))
	#    plot_df$method <- gsub("PARTICA", "CONFETI", as.character(plot_df$method))


	    #plot_df <- subset(plot_df, method %in% methods.to.plot)
	    plot_df[is.na(plot_df)] <- 0
	    if(!is.null(factor_sort)){
	        plot_df$method <- factor(plot_df$method, levels = factor_sort)
	    }

	    # Creating individual cis and trans replication plots
	    total_hit_plot <- ggplot(plot_df, aes(x = FDR, y = total_rep, col = method, group = method)) + 
	                  geom_line(size = 1) + labs(title = paste(data_prefix, "All")) + scale_color_manual(values = custom_pal)
	                  
	    total_hit_plot <- ggplot_add_theme(total_hit_plot)


	    cis_hit_plot <- ggplot(plot_df, aes(x = FDR, y = cis, col = method, group = method)) + 
	                  geom_line(size = 1) + labs(title = paste(data_prefix, "cis")) + scale_color_manual(values = custom_pal)
	                  

	    cis_hit_plot <- ggplot_add_theme(cis_hit_plot)

	    trans_hit_plot <-  ggplot(plot_df, aes(x = FDR, y = trans, col = method, group = method)) + 
	                  geom_line(size = 1) + labs(title = paste(data_prefix, "trans")) + scale_color_manual(values = custom_pal)

	    trans_hit_plot <- ggplot_add_theme(trans_hit_plot)

	    pdf(file = paste(figure.path, plot_prefix, "_combined_FDR_replication_plot.pdf", sep = ""), height = 5, width = 12)
	        total_hit_plot <- total_hit_plot + theme(legend.position = "none")
	        cis_hit_plot <- cis_hit_plot + theme(legend.position = "none")
	        trans_hit_plot <- trans_hit_plot + theme(legend.text = element_text(size = 10), legend.key.size = unit(1, "cm"))
	        multiplot(plotlist = list(total_hit_plot, cis_hit_plot,trans_hit_plot), layout = matrix(c(1,1,2,2,3,3,3), nrow = 1, ncol = 7))
	    dev.off()
    }





}
custom_pal <- confeti_palette()
#+ facet_grid(group~.)
# 
#+ 
#dataset.name <- "GTExEso"
custom_pal <- confeti_palette()
#+ facet_grid(group~.)
# 
#+ 
#dataset.name <- "GTExEso"
dataset.name <- as.character(commandArgs(TRUE)[1])

if(is.na(dataset.name)){

    dataset.name <- "GTEx"

}
replication_plots <- FALSE
pseudo_trans_plot = FALSE
pca_geno_correlations <- FALSE

merged.file.path <- "/zenodotus/dat01/mezeylab_scratch/jij2009/confeti_paper_results/eqtl_hits/GTExSkin_hits/"
eqtl_data_path <- "/zenodotus/dat01/mezeylab_scratch/jij2009/eqtl_input_h5/"
#merged.file.path <- "/zenodotus/dat01/mezeylab_scratch/jij2009/eqtl_multi_method/eqtl_merged_results/GTExSkin_eqtl_hits/"
#eqtl_data_file <- "/zenodotus/dat01/mezeylab_scratch/jij2009/eqtl_h5_with_Kmx/GTExEsophagusJunction.h5"
#eqtl_data_file <- "/zenodotus/dat01/mezeylab_scratch/jij2009/eqtl_h5_with_Kmx/GTExAdiposeSubq.h5"
eqtl_data <- grep(dataset.name, dir(eqtl_data_path), value = TRUE)
eqtl_data_file <- paste0(eqtl_data_path,eqtl_data[1])
#figure_path <- "/zenodotus/dat01/mezeylab_scratch/jij2009/confeti_paper_results/Figures/"
figure_path <- "/zenodotus/dat01/mezeylab_scratch/jij2009/confeti_paper_results/test_plots/indv_dataset/"
method_factors <- c("CONFETI", "CONPANA", "PANAMA", "ICE", "LINEAR")
#method_factors <- c("PARTICApy","CONPANApy", "PANAMApy","CONFETIpy", "ICEpy", "PCALMM2")
#method_factors <- c("CONPANApy","CONFETIpy","PARTICApy", "PANAMApy", "PCALMM2py","ICEpy", "PEER", "LINEAR")
#method_factors <- c("CONPANArawmx","PARTICArawmx", "PANAMArawmx", "ICErawmx","PCALMM2rawmx")
#method_factors <- c("CONPANAcoeff","PARTICAcoeff", "PANAMAcoeff", "ICEcoeff","PCALMM2coeff")
#method_factors <- c("CONFETIlrgk","CONPANAlrgk","PARTICAlrgk", "PANAMAlrgk", "ICElrgk","PCALMM2lrgk")
#method_factors <- c("CONFETIgpca","CONPANAgpca","PARTICAgpca", "PANAMAgpca", "ICEgpca","PCALMM2gpca")
bh_cutoff <- 0.025

#methods.to.plot <- c("FULLICAr", "ICALMMr", "PANAMAr", "ICEr","PCALMMr","LINEARr")

today_date <- gsub("-", "_", Sys.Date())
prefix <- paste("_", today_date, sep = "")

setwd(merged.file.path)

dataset.files <- grep(dataset.name, dir(), value = TRUE)
dataset.files <- grep("[.]h5", dataset.files, value = TRUE)

eqtl.merged.list <- list()

unique.datasets <- unique(sapply(strsplit(dataset.files,"_"), function(x) x[1]))
unique.methods <- unique(sapply(strsplit(dataset.files,"_"), function(x) x[2]))

if(is.null(method_factors)){

    method_factors <- unique.methods

}
##############################################################################
########################### Import eQTL Data  ################################
##############################################################################

for(i in 1:length(unique.methods)){
    revised_method <- gsub("py", "", gsub("[0-9]","",unique.methods[i]))
    if(revised_method %in% method_factors){
        message("Processing for method = ", unique.methods[i])
  
        method.files <- grep(paste(unique.methods[i],"_",sep=""), dataset.files, value = TRUE)

        for(f in 1:length(method.files)){
            message("Loading file = ", method.files[f])
            subdataset <- unlist(strsplit(method.files[f], "_"))[1]
            eqtl.merged.list[[subdataset]][[revised_method]] <- read_in_h5_eqtl(method.files[f])
            eqtl.merged.list[[subdataset]][[revised_method]]$method <- unique.methods[i]
            H5close()
        }  

    }

}
cat("Read in complete\n")

eqtl_fdr_cis_trans_plots(eqtl.merged.list, max.fdr =bh_cutoff, figure.path = figure_path, factor_sort = method_factors)