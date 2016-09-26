suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(icreport))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(rhdf5))


source("~/eqtl_multi_method/eqtl_hits_plotting_functions.R")
get.cband <- function(chrom, position) {
  if ( is.na(position) ) {
    return ( chrom )
  } else {
   return ( paste(chrom, chrom.bands$name[which.max(chrom.bands$chrom.num==chrom & position < chrom.bands$chromEnd)], sep='') )
  }
}
##############################################################################
############################ Start of Script #################################
##############################################################################

custom_pal <- confeti_palette()
#+ facet_grid(group~.)
# 
#+ 
#dataset.name <- "GTExEso"

dataset.name <- as.character(commandArgs(TRUE)[1])

dataset.name <- "GTEx"

pval_cutoff <- 0.01
#sig_hits_criteria <-  "bon"
sig_hits_criteria <-  "BH"

replication_plots <- TRUE
pseudo_trans_plot = TRUE
pca_geno_correlations <- FALSE
master_path <- "/zenodotus/dat01/mezeylab_scratch/jij2009/confeti_R_GPCA/"
eqtl_data_path <- "/zenodotus/dat01/mezeylab_scratch/jij2009/eqtl_input_h5/"
master_figure_path <- "/zenodotus/dat01/mezeylab_scratch/jij2009/confeti_R_GPCA/figures/"

figure_path <- paste0(master_figure_path, dataset.name, "_", sig_hits_criteria,"/")
eqtl_data_file <- list.files(eqtl_data_path, pattern = dataset.name, full.names = TRUE)
dir.create(figure_path, showWarnings = FALSE)

#figure_path <- paste0(master_path,"plots_all_methods_160912/")
method_factors <- c("CONFETI", "CONPANA", "PANAMA", "PCAKMX", "ICE", "LINEAR")


today_date <- gsub("-", "_", Sys.Date())

plot_prefix <- paste0(dataset.name, "_", today_date)


dataset.files <- list.files(master_path, recursive = TRUE, 
                            pattern = paste0(dataset.name,".*eqtldf[.]h5"), 
                            full.names = TRUE)

file_names_only <- sapply(strsplit(dataset.files, "/"), function(x) x[length(x)])


eqtl_hit_list <- list()

unique.datasets <- unique(sapply(strsplit(file_names_only,"_"), function(x) x[1]))
unique.methods <- unique(sapply(strsplit(file_names_only,"_"), function(x) x[2]))

if(is.null(method_factors)){

    method_factors <- unique.methods

}

total_tests <- list()

for(single_dataset in unique.datasets){
    data_file <- grep(single_dataset, eqtl_data_file, value = TRUE)
    n_geno <- length(h5read(data_file, "genotypes/col_info/id"))
    n_pheno <- length(h5read(data_file, "phenotypes/col_info/id"))
    total_tests[[single_dataset]] <- as.numeric(n_geno) * as.numeric(n_pheno)

}

##############################################################################
########################### Import eQTL Data  ################################
##############################################################################

for(i in 1:length(unique.methods)){
    revised_method <- gsub("Rgpca", "", unique.methods[i])
    revised_method <- gsub("py", "", gsub("[0-9]", "", revised_method))
    revised_method <- gsub("gpca", "", revised_method)

    if(revised_method %in% method_factors){

        message("Processing for method = ", unique.methods[i])
  
        method.files <- grep(paste(unique.methods[i],"_",sep=""), dataset.files, value = TRUE)

        for(f in 1:length(method.files)){
            
            file_prefix <- sapply(strsplit(method.files[f], "/eqtl_hits/"), function(x) x[2])
            message("Loading file = ", file_prefix)
            subdataset <- unlist(strsplit(file_prefix, "_"))[1]
            eqtl_hit_list[[revised_method]][[subdataset]] <- read_in_h5_eqtl(method.files[f])
        }


    }
      
}
cat("Read in complete\n")
# get unique genotypes
cat("Matching cytobands with genotypes \n")
load('~/eqtl_multi_method/cytoband/chrom_bands.RData')

geno_only <- plyr::ldply(lapply(eqtl_hit_list, function(x) 
                          plyr::ldply(lapply(x, function(a) 
                          subset(a, select = c(genotype, geno_chr, geno_pos))))))

geno_only <- geno_only[!duplicated(geno_only$genotype),]

geno_only$chrom_band <- apply(geno_only, 1, function(x) get.cband(x["geno_chr"], x["geno_pos"]))

for(m in names(eqtl_hit_list)){
    for (d in names(eqtl_hit_list[[m]])){

      eqtl_hit_list[[m]][[d]]$cytoband <- geno_only[match(eqtl_hit_list[[m]][[d]]$genotype, geno_only$genotype), "chrom_band"]
      eqtl_hit_list[[m]][[d]]$eqtl_id <- with(eqtl_hit_list[[m]][[d]], 
                                                 paste(pheno_symbol, cis_trans, cytoband, sep ="_"))
    }


}




##############################################################################
####################### Create Replication Plots #############################
##############################################################################
message("Creating replication plots by FDR, bonferroni and rank")

if(replication_plots){
    eqtl_rep_fraction_plots(eqtl_hit_list, max.fdr = 0.02, 
                   figure.path = figure_path, 
                   plot_prefix = plot_prefix, factor_sort = method_factors)


    eqtl_rep_bonferroni_plots(eqtl_hit_list, 
                          total_tests = total_tests,
                          max.bon = 0.2,
                          figure.path = figure_path,  
                          plot_prefix = plot_prefix, 
                          factor_sort = method_factors)

    eqtl_rep_rank_plots(eqtl_hit_list, max.rank = 10000, 
                    figure.path = figure_path, 
                    plot_prefix = plot_prefix,
                    factor_sort = method_factors)

}

message("Replication plots complete\n")


##############################################################################
####################### Choosing Significant Hits ############################
##############################################################################

message("Subsetting hits by ", sig_hits_criteria, " with a cutoff of ",pval_cutoff)

if(sig_hits_criteria == "BH"){

    sig_hit_list <- lapply(eqtl_hit_list, 
                       function(x) lapply(x, function(a) subset(a, a$p.bh < pval_cutoff)))


} else if (sig_hits_criteria == "bon"){

    sig_hit_list <- list()
    for(m in names(eqtl_hit_list)){
        for (d in names(eqtl_hit_list[[m]])){
            bon_threshold <- as.numeric(0.05 / total_tests[[d]])
            #cat("Total tests = ", total_tests[[d]], "\n")
            #cat("bon_threshold = ", -log10(bon_threshold),"\n")
            sig_hit_list[[m]][[d]] <- subset(eqtl_hit_list[[m]][[d]], eqtl_hit_list[[m]][[d]]$pval < bon_threshold)
        }

    }

}


# part where we extract the original hit information based on replicating hits
unique_rep_hits <- get_intersecting_hits(sig_hit_list)

extract_columns <- c("genotype", "phenotype", "geno_chr", 
                     "geno_pos", "pheno_chr", "pheno_start", 
                     "pheno_symbol", "cis_trans", "eqtl_id")
rep_info_df <- list()
rep_hit_list <- list()
for(single_method in names(unique_rep_hits)){
    cat("Processing method =", single_method, "\n")
#    rep_info_df[[single_method]] <- lapply(sig_hit_list[[single_method]], 
#                                    function(x) subset(x, x$eqtl_id %in% unique_rep_hits[[single_method]]))

#    rep_info_ids[[single_method]] <- lapply(rep_info_df[[single_method]], function(x) x[,"info_id"])
    rep_info_df[[single_method]] <- lapply(sig_hit_list[[single_method]], 
                                           function(x) subset(x, x$eqtl_id %in% unique_rep_hits[[single_method]])[, extract_columns])
    rep_hit_list[[single_method]] <- Reduce(rbind, rep_info_df[[single_method]])
}


if(pseudo_trans_plot){
    source("~/eqtl_multi_method/identify_pseudo_trans_qtl.r")
    #source("~/eqtl_multi_method/identify_pseudo_trans_qtl_blat.r")
    rep_hits_pseudo_labeled <- lapply(rep_hit_list, function(x) detect_pseudo_trans(x, entrez2genome))

    for(i in names(rep_hits_pseudo_labeled)){
        rep_hits_pseudo_labeled[[i]]$eqtl_pseudo_id <- with(rep_hits_pseudo_labeled[[i]], paste(eqtl_id, pseudo.cis.final, sep = "_")) 

    }

    rep_hit_reduced <- lapply(rep_hits_pseudo_labeled, function(x) x[!duplicated(x$eqtl_pseudo_id),])

    ##############################################################################
    #################### Generating Significant hit map ##########################
    ##############################################################################

    for( i in names(rep_hit_reduced)){
        sig_hits_map(rep_hit_reduced[[i]], 
                      prefix = paste0(plot_prefix,"_",i),
                      output_path = figure_path)
    }

    ##############################################################################
    #################### Generating Bar Plots for Hits  ##########################
    ##############################################################################

    trans_hits_only <- lapply(rep_hit_reduced, function(x) subset(x, x$cis_trans == "trans"))

    #trans_plot_data <- rep_hits_plot_data(input_list = trans_hits_only, plot_factors = method_factors)
    trans_plot_data <- rep_hits_plot_data_trans(input_list = trans_hits_only, plot_factors = method_factors)

    rep_hits_barplot_trans(trans_plot_data, 
                           figure_path =figure_path, 
                           prefix = plot_prefix, 
                           color_set = custom_pal, no_axis_ticks = FALSE)

    trans_plot_dummy <- subset(trans_plot_data, variable == "hits")
    trans_plot_dummy$hits <- 1

    gene_sort <- trans_plot_dummy$hits
    names(gene_sort) <- trans_plot_dummy$gene
    unique_genes <- unique(names(gene_sort))

    # merging counts for unique genes
    count_vec <- c()
    for( u in unique_genes){
        total_count <- sum(gene_sort[which(names(gene_sort) == u)])
        names(total_count) <- u
        count_vec <- c(count_vec, total_count)
    }

    trans_plot_dummy$gene <- factor(trans_plot_dummy$gene, 
                                    levels = unique(names(sort(count_vec, decreasing = TRUE))))

    rep_hits_matrix_plot(trans_plot_dummy, figure_path = figure_path, 
                     dataset_name =  plot_prefix,  cis_trans = "trans")

    rep_hits_barplot(trans_plot_dummy, figure_path = figure_path, 
                     prefix = plot_prefix, 
                     color_set = custom_pal)


    cis_hits_only <- lapply(rep_hit_reduced, function(x) subset(x, x$cis_trans == "cis"))

    cis_plot_data <- rep_hits_plot_data(input_list = cis_hits_only, plot_factors = method_factors)

    cis_plot_data$hits <- 1
    rep_hits_matrix_plot(cis_plot_data, figure_path = figure_path, 
                     dataset_name =  plot_prefix, cis_trans = "cis")

}


rep_trans_hit_id <- lapply(trans_hits_only, function(x) subset(x, x$pseudo.cis.final == FALSE)[,"eqtl_id"])

# extract trans hits that are replicating from original dataset
rep_trans_data <- list()
for(i in names(eqtl_hit_list)){

    for(d in names(eqtl_hit_list[[i]])){

        temp_df <- subset(eqtl_hit_list[[i]][[d]], eqtl_hit_list[[i]][[d]][,"eqtl_id"] %in% rep_trans_hit_id[[i]])
        rep_trans_data[[i]][[d]] <- temp_df[!duplicated(temp_df$eqtl_id),]
    }

}


trans_info_save_path <- "/zenodotus/dat01/mezeylab_scratch/jij2009/confeti_R_GPCA/trans_hit_inspection/"

save(rep_trans_data, file = paste0(trans_info_save_path, dataset.name, "_trans_inspect.RData"))
