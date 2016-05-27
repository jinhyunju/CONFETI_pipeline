# lamda plot generation script
library(rhdf5)
library(plyr)
library(ggplot2)

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

dataset_name <- as.character(commandArgs(TRUE)[1])
input_path <- as.character(commandArgs(TRUE))[2]
figure_path <- as.character(commandArgs(TRUE))[3]

message("Generating genomic inflation factor plots\n")

message("Figures generated in = ", figure_path, "\n")

if(length(commandArgs(TRUE)) < 1){

#  dataset_name <- "GTExAdiposeSubq"
  input_path <- "/zenodotus/dat01/mezeylab_scratch/jij2009/confeti_paper_results/pval_medians/"
  figure_path <- "/zenodotus/dat01/mezeylab_scratch/jij2009/confeti_paper_results/test_plots/lambda_plots/"

}
setwd(input_path)
method_factors = c("CONFETI", "CONPANA", "PANAMA", "ICE", "LINEAR")
#method_factors <- c("PARTHYBRIDpy", "PARTICApy", "PANAMApy", "ICEpy", "PEER", "LINEAR")
all_pval_median_files <- grep("pvals_median.h5", dir(), value = TRUE)


unique.datasets <- unique(sapply(strsplit(all_pval_median_files,"_"), function(x) x[1]))

for(dataset_name in unique.datasets){
    pval_median_files <- grep(dataset_name, all_pval_median_files, value = TRUE)

    pval_median_list <- list()

    unique.methods <- unique(sapply(strsplit(pval_median_files,"_"), function(x) x[2]))


    for(i in 1:length(unique.methods)){
      revised_method <- gsub("py", "", gsub("[0-9]","", unique.methods[i]))

      if(revised_method %in% method_factors){
          message("Processing for method = ", unique.methods[i])
      
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

    pdf(paste0(figure_path, "/",dataset_name, "_lambda_boxplot.pdf"), width = 6, height = 4)

      p <- ggplot(lambda_df, aes(x = method, y = lambda_diff, group = method, fill = method)) + geom_boxplot(outlier.shape = NA)  + theme_bw()
      p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none', plot.title = element_text(face = "bold")) + ylim(c(-1,1))+ scale_fill_manual(values = custom_pal)
      p <- p + geom_hline(aes(yintercept=0)) + ggtitle(paste0(dataset_name, "_Genomic_Inflation_Factor"))
      #p <- p + scale_fill_manual(values = custom_pal)
     print(p)
    dev.off()

}



#lambda_medians <- lapply(lambda_list, function(x) lapply(x, function(a) median(a, na.rm = TRUE)))
#lambda_median_plot <- reshape2::melt(ldply(lambda_medians, data.frame))
#colnames(lambda_median_plot) <- c("method", "dataset", "median") 

#lambda_median_plot$method <- factor(lambda_median_plot$method, levels = method_factors )

#lambda_median_plot$diff <- lambda_median_plot$median - 1

#pdf(paste0(figure_path, dataset_name, "_lambda_median.pdf"), width = 6, height = 3)

#  p <- ggplot(lambda_median_plot, aes(x = method, y = diff, group = method, fill = method)) + geom_bar(stat = 'identity') + facet_grid(.~dataset) + scale_fill_manual(values = custom_pal) + theme_bw()
#  p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none') + ylim(c(-0.25,0.25))
# print(p)
#dev.off()
