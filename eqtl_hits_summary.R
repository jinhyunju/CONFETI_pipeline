suppressPackageStartupMessages(library(ggplot2))

suppressPackageStartupMessages(library(icreport))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(rhdf5))
##############################################################################
########################## Defining Functions ################################
##############################################################################
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

find_intersecting_hits <- function(input_list, n.hits){
  #sig.df <- input.df[c(1:n.hits),]
  method.data.split <- lapply(input_list, function(a) lapply(a, function(x) x[1:n.hits,]))
  unique.hit.list <- lapply(method.data.split, function(x) lapply(x, function(a) unique(a$eqtl_id)))
  intersecting_hits <- lapply(unique.hit.list, function(x) Reduce(intersect, x))
  return(intersecting_hits)
}

eqtl_rep_fraction_plots <- function(eqtl.merged.list, 
                                   max.fdr = 0.1, 
                                   figure.path, 
                                   data_prefix, 
                                   factor_sort = NULL){

    today_date <- gsub("-", "_", Sys.Date())
    plot_prefix <- paste(data_prefix, today_date, sep = "_")

    max.pval.data <- sapply(eqtl.merged.list, function(a) sapply(a, function(b) max(range(b$p.bh), na.rm = TRUE)))
    max.pval.data <- max(unlist(max.pval.data), na.rm = TRUE)

    if(max.pval.data < max.fdr){

        max.fdr <- max.pval.data

    }

    FDR_range <- seq(0, max.fdr, length.out = 100)

    cis_trans_list <- list()
    union_count_list <- list()

    for( c in 1:length(FDR_range)) {

            significant_unique <- lapply(eqtl.merged.list, function(a) lapply(a, function(b) unique(subset(b, b$p.bh < FDR_range[c])$eqtl_id) ) ) # don't forget to change it to eqtl_id
            intersecting_hits <- lapply(significant_unique, function(x) Reduce(intersect, x))
            union_hits <- lapply(significant_unique, function(x) unique(Reduce(c, x)) )
            ct_count_df <- data.frame(melt(sapply(intersecting_hits, cis_trans_counter)))
            ct_count_df$FDR <- FDR_range[c]
            union_count_df <- data.frame(melt(sapply(union_hits, cis_trans_counter)))
            union_count_df$FDR <- FDR_range[c]
            colnames(ct_count_df) <- c("cis_trans", "method", "count", "FDR")
            colnames(union_count_df) <- c("cis_trans", "method", "union", "FDR")
            union_count_list[[c]] <- union_count_df
            cis_trans_list[[c]] <- ct_count_df
    }

    cis_trans_df <- Reduce(rbind, lapply(cis_trans_list, function(x) x))
    union_count_df <- Reduce(rbind, lapply(union_count_list, function(x) x))

    cis_trans_cast <- dcast(cis_trans_df, method + FDR ~ cis_trans, value.var = "count")
    union_cast <- dcast(union_count_df, method + FDR ~ cis_trans, value.var = "union")
    colnames(union_cast) <- c("method", "FDR", "cis_union", "trans_union")
    plot_df <- merge(cis_trans_cast, union_cast)
    #plot_df <- subset(plot_df, FDR > 0)


    plot_df$total_rep <- with(plot_df, cis + trans)
    plot_df$total_union <- with(plot_df, cis_union + trans_union)
    plot_df$total_rep_ratio <- with(plot_df, total_rep / total_union)

    plot_df$cis_rep_ratio <- with(plot_df, cis / cis_union )
    plot_df$trans_rep_ratio <- with(plot_df, trans / trans_union)

#    plot_df$method <- gsub("py", "", as.character(plot_df$method))
#    plot_df$method <- gsub("PARTICA", "CONFETI", as.character(plot_df$method))


    #plot_df <- subset(plot_df, method %in% methods.to.plot)
    plot_df[is.na(plot_df)] <- 0
    if(!is.null(factor_sort)){
        plot_df$method <- factor(plot_df$method, levels = factor_sort)
    }

    # Creating individual cis and trans replication plots
    total_rep_plot <- ggplot(plot_df, aes(x = FDR, y = total_rep, col = method, group = method)) + 
                  geom_line(size = 1) + labs(title = paste(data_prefix, "All")) + scale_color_manual(values = custom_pal)
                  
    total_rep_plot <- ggplot_add_theme(total_rep_plot)


    cis_rep_plot <- ggplot(plot_df, aes(x = FDR, y = cis, col = method, group = method)) + 
                  geom_line(size = 1) + labs(title = paste(data_prefix, "cis")) + scale_color_manual(values = custom_pal)
                  

    cis_rep_plot <- ggplot_add_theme(cis_rep_plot)

    trans_rep_plot <-  ggplot(plot_df, aes(x = FDR, y = trans, col = method, group = method)) + 
                  geom_line(size = 1) + labs(title = paste(data_prefix, "trans")) + scale_color_manual(values = custom_pal)

    trans_rep_plot <- ggplot_add_theme(trans_rep_plot)

    pdf(file = paste(figure.path, plot_prefix, "_combined_FDR_replication_plot.pdf", sep = ""), height = 5, width = 12)
        total_rep_plot <- total_rep_plot + theme(legend.position = "none")
        cis_rep_plot <- cis_rep_plot + theme(legend.position = "none")
        trans_rep_plot <- trans_rep_plot + theme(legend.text = element_text(size = 10), legend.key.size = unit(1, "cm"))
        multiplot(plotlist = list(total_rep_plot, cis_rep_plot,trans_rep_plot), layout = matrix(c(1,1,2,2,3,3,3), nrow = 1, ncol = 7))
    dev.off()


    total_ratio_plot <- ggplot(plot_df, aes(x = FDR, y = total_rep_ratio, col = method, group = method)) + 
                  geom_line(size = 1) + labs(title = paste(data_prefix, "All")) + scale_color_manual(values = custom_pal)
                  
    total_ratio_plot <- ggplot_add_theme(total_ratio_plot)


    cis_ratio_plot <- ggplot(plot_df, aes(x = FDR, y = cis_rep_ratio, col = method, group = method)) + 
                  geom_line(size = 1) + labs(title = paste(data_prefix, "cis")) + scale_color_manual(values = custom_pal)
                  

    cis_ratio_plot <- ggplot_add_theme(cis_ratio_plot)

    trans_ratio_plot <-  ggplot(plot_df, aes(x = FDR, y = trans_rep_ratio, col = method, group = method)) + 
                  geom_line(size = 1) + labs(title = paste(data_prefix, "trans")) + scale_color_manual(values = custom_pal)

    trans_ratio_plot <- ggplot_add_theme(trans_ratio_plot)

    pdf(file = paste(figure.path, plot_prefix, "_combined_FDR_rep_ratio_plot.pdf", sep = ""), height = 5, width = 12)
        total_ratio_plot <- total_ratio_plot + theme(legend.position = "none")
        cis_ratio_plot <- cis_ratio_plot + theme(legend.position = "none")
        trans_ratio_plot <- trans_ratio_plot + theme(legend.text = element_text(size = 10), legend.key.size = unit(1, "cm"))
        multiplot(plotlist = list(total_ratio_plot, cis_ratio_plot,trans_ratio_plot), layout = matrix(c(1,1,2,2,3,3,3), nrow = 1, ncol = 7))
    dev.off()
    
    total_union_plot <- ggplot(plot_df, aes(x = FDR, y = total_union, col = method, group = method)) + 
                  geom_line(size = 1) + labs(title = paste(data_prefix, "All")) + scale_color_manual(values = custom_pal)
                  
    total_union_plot <- ggplot_add_theme(total_union_plot)


    cis_union_plot <- ggplot(plot_df, aes(x = FDR, y = cis_union, col = method, group = method)) + 
                  geom_line(size = 1) + labs(title = paste(data_prefix, "cis")) + scale_color_manual(values = custom_pal)
                  

    cis_union_plot <- ggplot_add_theme(cis_union_plot)

    trans_union_plot <-  ggplot(plot_df, aes(x = FDR, y = trans_union, col = method, group = method)) + 
                  geom_line(size = 1) + labs(title = paste(data_prefix, "trans")) + scale_color_manual(values = custom_pal)

    trans_union_plot <- ggplot_add_theme(trans_union_plot)

    pdf(file = paste(figure.path, plot_prefix, "_combined_FDR_union_plot.pdf", sep = ""), height = 5, width = 12)
        total_union_plot <- total_union_plot + theme(legend.position = "none")
        cis_union_plot <- cis_union_plot + theme(legend.position = "none")
        trans_union_plot <- trans_union_plot + theme(legend.text = element_text(size = 10), legend.key.size = unit(1, "cm"))
        multiplot(plotlist = list(total_union_plot, cis_union_plot,trans_union_plot), layout = matrix(c(1,1,2,2,3,3,3), nrow = 1, ncol = 7))
    dev.off()

}

eqtl_rep_bonferroni_plots <- function(eqtl.merged.list, 
                                      total_tests = NULL,
                                      max.bon = 0.2,
                                      figure.path, 
                                      data_prefix, 
                                      factor_sort = NULL){
    cat("Generate bonferroni rep plots\n")
    today_date <- gsub("-", "_", Sys.Date())
    plot_prefix <- paste(data_prefix, today_date, sep = "_")

    max.pval.data <- sapply(eqtl.merged.list, function(a) sapply(a, function(b) max(range(b$pval), na.rm = TRUE)))
    max.pval.data <- max(unlist(max.pval.data), na.rm = TRUE)

    if(max.pval.data < max.bon){

        max.bon <- max.pval.data

    }

    bon_range <- seq(0, max.bon, length.out = 100)
    bon_range <- bon_range[-1]

    cis_trans_list <- list()
    union_count_list <- list()

    for( c in 1:length(bon_range)) {
            significant_unique <- list()
            for(m in names(eqtl.merged.list)){
                for (d in names(eqtl.merged.list[[m]])){
                    bon_threshold <- as.numeric(bon_range[c] / total_tests[[d]])
                    #cat("Total tests = ", total_tests[[d]], "\n")
                    #cat("bon_threshold = ", -log10(bon_threshold),"\n")
                    significant_unique[[m]][[d]] <- unique(subset(eqtl.merged.list[[m]][[d]], eqtl.merged.list[[m]][[d]]$pval < bon_threshold)$eqtl_id)
                }

            }
            #significant_unique <- lapply(eqtl.merged.list, function(a) lapply(a, function(b) unique(subset(b, b$p.bh < FDR_range[c])$eqtl_id) ) ) # don't forget to change it to eqtl_id
            intersecting_hits <- lapply(significant_unique, function(x) Reduce(intersect, x))
            union_hits <- lapply(significant_unique, function(x) unique(Reduce(c, x)) )
            ct_count_df <- data.frame(melt(sapply(intersecting_hits, cis_trans_counter)))
            ct_count_df$bon <- bon_range[c]
            union_count_df <- data.frame(melt(sapply(union_hits, cis_trans_counter)))
            union_count_df$bon <- bon_range[c]
            colnames(ct_count_df) <- c("cis_trans", "method", "count", "bon")
            colnames(union_count_df) <- c("cis_trans", "method", "union", "bon")
            union_count_list[[c]] <- union_count_df
            cis_trans_list[[c]] <- ct_count_df
    }

    cis_trans_df <- Reduce(rbind, lapply(cis_trans_list, function(x) x))
    union_count_df <- Reduce(rbind, lapply(union_count_list, function(x) x))

    cis_trans_cast <- dcast(cis_trans_df, method + bon ~ cis_trans, value.var = "count")
    union_cast <- dcast(union_count_df, method + bon ~ cis_trans, value.var = "union")
    colnames(union_cast) <- c("method", "bon", "cis_union", "trans_union")
    plot_df <- merge(cis_trans_cast, union_cast)
    #plot_df <- subset(plot_df, FDR > 0)


    plot_df$total_rep <- with(plot_df, cis + trans)
    plot_df$total_union <- with(plot_df, cis_union + trans_union)
    plot_df$total_rep_ratio <- with(plot_df, total_rep / total_union)

    plot_df$cis_rep_ratio <- with(plot_df, cis / cis_union )
    plot_df$trans_rep_ratio <- with(plot_df, trans / trans_union)

#    plot_df$method <- gsub("py", "", as.character(plot_df$method))
#    plot_df$method <- gsub("PARTICA", "CONFETI", as.character(plot_df$method))


    #plot_df <- subset(plot_df, method %in% methods.to.plot)
    plot_df[is.na(plot_df)] <- 0
    if(!is.null(factor_sort)){
        plot_df$method <- factor(plot_df$method, levels = factor_sort)
    }

    # Creating individual cis and trans replication plots
    total_rep_plot <- ggplot(plot_df, aes(x = bon, y = total_rep, col = method, group = method)) + 
                  geom_line(size = 1) + labs(title = paste(data_prefix, "All hits")) + scale_color_manual(values = custom_pal)
                  
    total_rep_plot <- ggplot_add_theme(total_rep_plot)


    cis_rep_plot <- ggplot(plot_df, aes(x = bon, y = cis, col = method, group = method)) + 
                  geom_line(size = 1) + labs(title = paste(data_prefix, "cis hits")) + scale_color_manual(values = custom_pal)
                  

    cis_rep_plot <- ggplot_add_theme(cis_rep_plot)

    trans_rep_plot <-  ggplot(plot_df, aes(x = bon, y = trans, col = method, group = method)) + 
                  geom_line(size = 1) + labs(title = paste(data_prefix, "trans hits")) + scale_color_manual(values = custom_pal)

    trans_rep_plot <- ggplot_add_theme(trans_rep_plot)

    pdf(file = paste(figure.path, plot_prefix, "_combined_bon_replication_plot.pdf", sep = ""), height = 5, width = 10)
        total_rep_plot <- total_rep_plot + theme(legend.position = "none")
        cis_rep_plot <- cis_rep_plot + theme(legend.position = "none")
        trans_rep_plot <- trans_rep_plot + theme(legend.text = element_text(size = 12), legend.key.size = unit(1, "cm"))
        multiplot(plotlist = list(total_rep_plot, cis_rep_plot,trans_rep_plot), layout = matrix(c(1,1,2,2,3,3,3), nrow = 1, ncol = 7))
    dev.off()


    total_ratio_plot <- ggplot(plot_df, aes(x = bon, y = total_rep_ratio, col = method, group = method)) + 
                  geom_line(size = 1) + labs(title = paste(data_prefix, "All hits")) + scale_color_manual(values = custom_pal)
                  
    total_ratio_plot <- ggplot_add_theme(total_ratio_plot)


    cis_ratio_plot <- ggplot(plot_df, aes(x = bon, y = cis_rep_ratio, col = method, group = method)) + 
                  geom_line(size = 1) + labs(title = paste(data_prefix, "cis hits")) + scale_color_manual(values = custom_pal)
                  

    cis_ratio_plot <- ggplot_add_theme(cis_ratio_plot)

    trans_ratio_plot <-  ggplot(plot_df, aes(x = bon, y = trans_rep_ratio, col = method, group = method)) + 
                  geom_line(size = 1) + labs(title = paste(data_prefix, "trans hits")) + scale_color_manual(values = custom_pal)

    trans_ratio_plot <- ggplot_add_theme(trans_ratio_plot)

    pdf(file = paste(figure.path, plot_prefix, "_combined_bon_rep_ratio_plot.pdf", sep = ""), height = 5, width = 10)
        total_ratio_plot <- total_ratio_plot + theme(legend.position = "none")
        cis_ratio_plot <- cis_ratio_plot + theme(legend.position = "none")
        trans_ratio_plot <- trans_ratio_plot + theme(legend.text = element_text(size = 12), legend.key.size = unit(1, "cm"))
        multiplot(plotlist = list(total_ratio_plot, cis_ratio_plot,trans_ratio_plot), layout = matrix(c(1,1,2,2,3,3,3), nrow = 1, ncol = 7))
    dev.off()
    
    total_union_plot <- ggplot(plot_df, aes(x = bon, y = total_union, col = method, group = method)) + 
                  geom_line(size = 1) + labs(title = paste(data_prefix, "All hits")) + scale_color_manual(values = custom_pal)
                  
    total_union_plot <- ggplot_add_theme(total_union_plot)


    cis_union_plot <- ggplot(plot_df, aes(x = bon, y = cis_union, col = method, group = method)) + 
                  geom_line(size = 1) + labs(title = paste(data_prefix, "cis hits")) + scale_color_manual(values = custom_pal)
                  

    cis_union_plot <- ggplot_add_theme(cis_union_plot)

    trans_union_plot <-  ggplot(plot_df, aes(x = bon, y = trans_union, col = method, group = method)) + 
                  geom_line(size = 1) + labs(title = paste(data_prefix, "trans hits")) + scale_color_manual(values = custom_pal)

    trans_union_plot <- ggplot_add_theme(trans_union_plot)

    pdf(file = paste(figure.path, plot_prefix, "_combined_bon_union_plot.pdf", sep = ""), height = 5, width = 10)
        total_union_plot <- total_union_plot + theme(legend.position = "none")
        cis_union_plot <- cis_union_plot + theme(legend.position = "none")
        trans_union_plot <- trans_union_plot + theme(legend.text = element_text(size = 12), legend.key.size = unit(1, "cm"))
        multiplot(plotlist = list(total_union_plot, cis_union_plot,trans_union_plot), layout = matrix(c(1,1,2,2,3,3,3), nrow = 1, ncol = 7))
    dev.off()

}

eqtl_rep_rank_plots <- function(eqtl.merged.list, 
                               max.rank = 10000, 
                               figure.path, 
                               data_prefix, 
                               factor_sort = NULL){
    today_date <- gsub("-", "_", Sys.Date())
    plot_prefix <- paste(data_prefix, today_date, sep = "_")
    
     # getting max rank of data
    data.max.rank <- max(unlist(lapply(eqtl.merged.list, function(x) lapply(x, nrow))))

    if(is.null(max.rank)){
        message("Using maximum number of hits as number to plot \n")
        max.rank <- data.max.rank
    } else {

        if( data.max.rank < max.rank){

            message("Specified max.rank is too large, using data.max.rank instead \n")
            max.rank <- data.max.rank
        }


    }
    #hit.range <- round(seq(100, min.number.hits, length.out = 100))
    hit.range <- round(seq(10, max.rank, length.out = 100))

    message("Creating Rank replication plot \n")
    intersect.short.list <- lapply(hit.range, function(x) find_intersecting_hits(eqtl.merged.list, x))

    short_ct_list <- NULL

    for(i in 1:length(intersect.short.list)){
      cutoff <- hit.range[i]
      temp <- NULL
      for(j in 1:length(intersect.short.list[[i]])){
        method <- names(intersect.short.list[[i]])[j]
        rep_count <- length(intersect.short.list[[i]][[j]])
        cis_trans <- sapply(strsplit(unique(intersect.short.list[[i]][[j]]), "_"), function(x) x[2])
        cis_count <- length(which(cis_trans == "cis"))
        trans_count <- length(which(cis_trans == "trans"))
        temp <- rbind(temp,c(cutoff, method, rep_count, cis_count, trans_count))
        
      }
      short_ct_list <- rbind(short_ct_list, temp)
      
    }

    message("Passed checkpoint 3 \n")
    colnames(short_ct_list) <- c("rank", "method", "total_rep", "cis_rep", "trans_rap")

    rank_plot <- data.frame("rank" = as.numeric(short_ct_list[,1]), 
                           "method" = as.character(short_ct_list[,2]), 
                           "total_rep" = as.numeric(short_ct_list[,3]),
                           "cis_rep" = as.numeric(short_ct_list[,4]), 
                           "trans_rep" = as.numeric(short_ct_list[,5]))

    if(!is.null(factor_sort)){
        rank_plot$method <- factor(rank_plot$method, levels = factor_sort)
    }
    method.number.of.hits = lapply(eqtl.merged.list, function(x) max(sapply(x, function(a) nrow(a) )))

    # to stop the rank plot if rank exceedes maximum number of hits
    for(method.names in names(method.number.of.hits)){
            rank_plot[(rank_plot$method == method.names & rank_plot$rank > method.number.of.hits[method.names]),"rank"] <- method.number.of.hits[method.names] 

    }
    rank_plot$rep_ratio <- (rank_plot$total_rep / rank_plot$rank)
    #rank_plot$method <- factor(rank_plot$method, levels = methods.to.plot)


    rank_total_replication <- ggplot(rank_plot, aes(x = rank, y = total_rep, col = method, group = method)) + 
                geom_line(size = 1) + labs(title = paste(data_prefix, "all rank rep"))+ scale_color_manual(values = custom_pal)

    rank_total_replication <- ggplot_add_theme(rank_total_replication)



    rank_cis_replication <- ggplot(rank_plot, aes(x = rank, y = cis_rep, col = method, group = method)) + 
                geom_line(size = 1) + labs(title = paste(data_prefix, "cis rank rep"))+ scale_color_manual(values = custom_pal)

    rank_cis_replication <- ggplot_add_theme(rank_cis_replication)



    rank_trans_replication <- ggplot(rank_plot, aes(x = rank, y = trans_rep, col = method, group = method)) + 
                geom_line(size = 1) + labs(title = paste(data_prefix, "trans rank rep"))+ scale_color_manual(values = custom_pal)

    rank_trans_replication <- ggplot_add_theme(rank_trans_replication)


    rank_rep_ratio <- ggplot(rank_plot, aes(x = rank, y = rep_ratio, col = method, group = method)) + 
                geom_line(size = 1) + labs(title = paste(data_prefix, "rank rep ratio"))+ scale_color_manual(values = custom_pal)

    rank_rep_ratio <- ggplot_add_theme(rank_rep_ratio)

    pdf(file = paste(figure.path, plot_prefix, "_combined_rank_replication_plot.pdf", sep = ""), height = 5, width = 10)
    rank_total_replication <- rank_total_replication + theme(legend.position = "none")
    rank_cis_replication <- rank_cis_replication + theme(legend.position = "none")
    rank_trans_replication <- rank_trans_replication + theme(legend.text = element_text(size = 12), legend.key.size = unit(1, "cm"))

    multiplot(plotlist = list(rank_total_replication, rank_cis_replication,rank_trans_replication), layout = matrix(c(1,1,2,2,3,3,3), nrow = 1, ncol = 7))
    dev.off()
 
}
# Function for creating significant hits plot
sig_hits_plot <- function(hits.info, prefix = "test", output_path = "./",
                          chr_info_file = '/zenodotus/dat01/mezeylab_scratch/jij2009/General_information/chromInfo.txt.gz'){

    # input file was downloaded from the following url
    # http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/chromInfo.txt.gz
    chrom.lengths <- read.table(chr_info_file, col.names=c('chromosome', 'size', 'file'))
    use.chromosomes <-  paste0('chr', c(1:22, 'X', 'Y', 'M'))
    chrom.lengths <- chrom.lengths[match(use.chromosomes, chrom.lengths$chromosome), "size"]

    snp.chroms <- gsub("chr", "", use.chromosomes)
    gene.chroms <- snp.chroms

    snp.lengths <- as.numeric(chrom.lengths)
    gene.lengths <- as.numeric(chrom.lengths)
    snp.mins <- gene.mins <- 0


    snp.offsets <- (cumsum(c(0, snp.lengths))[-length(snp.chroms)]-snp.mins)[1:length(snp.chroms)]
    gene.offsets <- (cumsum(c(0, gene.lengths))[-length(gene.chroms)]-gene.mins)[1:length(gene.chroms)]

    max.snp <- sum(snp.lengths)
    max.gene <- sum(gene.lengths)

    snp.ticks <- c(cumsum(c(0, snp.lengths)))
    gene.ticks <- c(cumsum(c(0, gene.lengths)))

    snp.labs <- (snp.ticks[-1] + snp.ticks[-length(snp.ticks)])/2
    gene.labs <- (gene.ticks[-1] + gene.ticks[-length(gene.ticks)])/2

    rbPal <- colorRampPalette(c('red','black'))
    # add minimum representative value to pvalues to avoid 0 -> INF conversion
    hits.info$col <- ifelse(hits.info$cis_trans == "cis", "black", "blue")
    #rbPal(2)[as.numeric(cut(c(1:nrow(hits.info$cis_trans)),breaks = 10))]
    hits.info$shape <- ifelse(hits.info$pseudo.cis.final, 8, 20)
    hits.info$col[hits.info$pseudo.cis.final] <- "red"

#    message("Creating Significant Hits Map \n")
#    message("File prefix = ", prefix, "\n")
#    message("Output directory = ", output_path, "\n")
    # Set up plot
    pdf(paste(output_path,prefix, "_rep_hits_map.pdf", sep = ""))
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

# Load data saved in hdf5 format 
load_data <- function(input_file = NULL, data_type = c("phenotypes","genotypes", "covars")){
    data_matrix <- h5read(input_file, paste0(data_type,"/matrix"))
    rownames(data_matrix) <- h5read(input_file, paste0(data_type,"/row_info/id"))
    colnames(data_matrix) <- h5read(input_file, paste0(data_type,"/col_info/id"))
    return(data_matrix)
}

# parse input ids and create a dataframe 
create_rep_hit_df <- function(input_ids){
    split_id_list <- strsplit(input_ids, "::")
    genotype <- sapply(split_id_list, function(x) x[1])
    phenotype <- sapply(split_id_list, function(x) x[2])
    geno_chr <- sapply(split_id_list, function(x) x[3])
    geno_pos <- sapply(split_id_list, function(x) x[4])
    pheno_chr <- sapply(split_id_list, function(x) x[5])
    pheno_start <- sapply(split_id_list, function(x) x[6])
    pheno_symbol <- sapply(split_id_list, function(x) x[7])
    cis_trans <- sapply(split_id_list, function(x) x[8])
    eqtl_id <- paste(pheno_symbol, cis_trans, pheno_chr, geno_chr, sep = "_")
    output_rep_df <- data.frame("genotype" = genotype,
                                 "phenotype" = phenotype, 
                                 "geno_chr" = as.numeric(geno_chr), 
                                 "geno_pos" = as.numeric(geno_pos), 
                                 "pheno_chr" = as.numeric(pheno_chr), 
                                 "pheno_start" = as.numeric(pheno_start), 
                                 "pheno_symbol" = pheno_symbol, 
                                 "cis_trans" = cis_trans, 
                                 "eqtl_id" = eqtl_id, stringsAsFactors = FALSE)
    return(output_rep_df)
    
}

# read in eqtl results saved in hdf5 format
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

get_intersecting_hits <- function(input_list){
    
    unique_id_list <- lapply(input_list, function(a) lapply(a, function(b) as.character(unique(b$eqtl_id)) ) )

    intersecting_unique <- lapply(unique_id_list, function(x) Reduce(intersect, x))
    
    return(intersecting_unique)
    
}

rep_hits_plot_data <- function(input_list, plot_factors){
    gene_trans_count <- lapply(input_list, function(x) sort(table(x$pheno_symbol), decreasing = TRUE))

    hit_count_list <- lapply(gene_trans_count , function(x) cbind(hits = x, gene = names(x)))
    
    list_names <- names(hit_count_list)
    lns <- sapply(hit_count_list, nrow)
    plot_data <- as.data.frame(do.call("rbind", hit_count_list))
    plot_data$group <- rep(list_names, lns)
    plot_data$hits <- as.numeric(plot_data$hits)
    plot_data$group <- factor(plot_data$group, 
                              levels = plot_factors)
    gene_sort <- plot_data$hits
    names(gene_sort) <- rownames(plot_data)
    unique_genes <- unique(names(gene_sort))

    # merging counts for unique genes
    count_vec <- c()
    for( u in unique_genes){
        total_count <- sum(gene_sort[which(names(gene_sort) == u)])
        names(total_count) <- u
        count_vec <- c(count_vec, total_count)
    }

    plot_data$gene <- factor(rownames(plot_data), levels = unique(names(sort(count_vec, decreasing = TRUE))))

    return(plot_data)
    
}

rep_hits_barplot <- function(plot_data, figure_path, prefix, color_set = custom_pal, no_axis_ticks = TRUE, no_gene_names = FALSE){
    pdf(paste(figure_path, prefix,"_rep_hits.pdf", sep = ""), height = 10, width = 30)
    #plot(1:100, 1:100)
    p <- ggplot(plot_data, aes(x = gene, y = hits, fill = group)) + geom_bar(stat = 'identity', position = 'dodge') + theme_bw()
    p <- p + facet_grid(group~.) + scale_fill_manual(values = color_set) 
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size=8), axis.title.y=element_blank()) 
    if(no_axis_ticks){
         p <- p + theme(axis.ticks.y=element_blank() ,axis.text.y = element_blank())
    }
    if(no_gene_names){

        p <- p + theme(axis.ticks.x=element_blank() ,axis.text.x = element_blank())
    }
    print(p)
    dev.off()
}

count_pseudo_hits <- function(input_list){
    count_table_list <- lapply(input_list, function(x) table(x$pheno_symbol, x$pseudo.cis.final))
    temp_list <- list()
    for( i in names(count_table_list)){
#        colnames(gene_trans_count[[i]]) <- c("hits", "pseudo_hits")

        if(dim(count_table_list[[i]])[1] == 0  ){
            temp_df <- data.frame("hits" = 0, "pseudo_hits" = 0, "gene" = NA)
            temp_list[[i]] <- temp_df

        } else if (dim(count_table_list[[i]])[2] == 1 ){
            
            if(colnames(count_table_list[[i]]) == "FALSE"){
                temp_df <- data.frame("hits" = 1, "pseudo_hits" = 1, "gene" = rownames(count_table_list[[i]]))
            } else {
                temp_df <- data.frame("hits" = 1, "pseudo_hits" = 0, "gene" = rownames(count_table_list[[i]]))
            }
            temp_list[[i]] <- temp_df

        } else {
            gene_names <- rownames(count_table_list[[i]])
            hit_count <- count_table_list[[i]][,1]
            pseudo_count <- count_table_list[[i]][,2]
            total_count <- hit_count + pseudo_count
            temp_df <- data.frame("hits" = hit_count, "pseudo_hits" = pseudo_count, "total_hits"=total_count, "gene" = gene_names)
            temp_df <- temp_df[order(temp_df$total_hits, decreasing = TRUE),]
            temp_df$total_hits <- NULL
            temp_list[[i]] <- temp_df

        }

    }

    return(temp_list)
}

rep_hits_plot_data_trans <- function(input_list, plot_factors){

    hit_count_list <- count_pseudo_hits(input_list)

    list_names <- names(hit_count_list)
    lns <- sapply(hit_count_list, nrow)
    plot_data <- as.data.frame(do.call("rbind", hit_count_list))
    plot_data$group <- rep(list_names, lns)
    plot_data$hits <- as.numeric(plot_data$hits)
    plot_data$group <- factor(plot_data$group, 
                              levels = plot_factors)
    #plot_data <- na.omit(plot_data)
    gene_sort <- plot_data$hits
    names(gene_sort) <- as.character(plot_data$gene)
    unique_genes <- unique(na.omit(plot_data$gene))

    # merging counts for unique genes
    count_vec <- c()
    for( u in unique_genes){
        
        total_count <- sum(gene_sort[which(names(gene_sort) == u)])
        names(total_count) <- u
        count_vec <- c(count_vec, total_count)
    }

    plot_data$gene <- factor(plot_data$gene, levels = unique(names(sort(count_vec, decreasing = TRUE))))

    melt_plot_data <- melt(plot_data, id.var = c("gene", "group"))
    melt_plot_data <- na.omit(melt_plot_data)
    return(melt_plot_data)
    
}



rep_hits_barplot_trans <- function(plot_data, figure_path, prefix, color_set = custom_pal, no_axis_ticks = TRUE){
    pdf(paste(figure_path, prefix,"_rep_hits.pdf", sep = ""), height = 10, width = 20)
    #plot(1:100, 1:100)
    p <- ggplot(plot_data, aes(x = gene, y = value, fill = variable)) + geom_bar( stat = 'identity') + theme_bw()
    p <- p + facet_grid(group~.) + scale_fill_manual(values = color_set) 
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size=8), axis.title.y=element_blank()) 
    if(no_axis_ticks){
         p <- p + theme(axis.ticks.y=element_blank() ,axis.text.y = element_blank())
    }
    print(p)
    dev.off()
}

convert_data_to_matrix <- function(input_df){

    plot_matrix <- reshape2::acast(input_df, group ~ gene, value.var = "hits")
    plot_matrix[is.na(plot_matrix)] <- 0
    plot_matrix <- plot_matrix * row(plot_matrix)

    return(plot_matrix)
}


rep_hits_matrix_plot <- function(plot_data, figure_path,
                                 color_set = custom_pal, dataset_name, bh_cutoff, cis_trans){
    file_prefix = paste(dataset_name,"fdr",gsub("[.]","",bh_cutoff),cis_trans,sep = "_")

    plot_matrix <- convert_data_to_matrix(plot_data)
    png(paste(figure_path, file_prefix,"_rep_hits.png", sep = ""), width = 1000, height = 250)
    par(mar = c(2,10,4,6))
    image(1:ncol(plot_matrix), 1:nrow(plot_matrix), t(plot_matrix[nrow(plot_matrix):1, ]), 
          axes = FALSE, col = c("#FFFFFF", color_set)[1:(nrow(plot_matrix)+1)], ann = FALSE) 
        # ann suppresses axis labels 
        # las flips labels
    axis(2, nrow(plot_matrix):1, rownames(plot_matrix), las = 1, lwd = 0, lwd.ticks = 1)
    axis(4, nrow(plot_matrix):1, apply(plot_matrix>0,1,sum), las = 1, lwd = 0)
    axis(4, nrow(plot_matrix) + 1, "Total", las = 1, lwd = 0, xpd = NA) #xpd clipping plot device
    title(main = paste(dataset_name, "Replicating",cis_trans,"hits"))
    dev.off()
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

if(is.na(dataset.name)){

    dataset.name <- "GTExHeart"

}
replication_plots <- FALSE
pseudo_trans_plot = TRUE
pca_geno_correlations <- FALSE

merged.file.path <- "/zenodotus/dat01/mezeylab_scratch/jij2009/confeti_paper_results/eqtl_hits/GTExHeart_hits/"
eqtl_data_path <- "/zenodotus/dat01/mezeylab_scratch/jij2009/eqtl_input_h5/"
#merged.file.path <- "/zenodotus/dat01/mezeylab_scratch/jij2009/eqtl_multi_method/eqtl_merged_results/GTExSkin_eqtl_hits/"
#eqtl_data_file <- "/zenodotus/dat01/mezeylab_scratch/jij2009/eqtl_h5_with_Kmx/GTExEsophagusJunction.h5"
#eqtl_data_file <- "/zenodotus/dat01/mezeylab_scratch/jij2009/eqtl_h5_with_Kmx/GTExAdiposeSubq.h5"
eqtl_data <- grep(dataset.name, dir(eqtl_data_path), value = TRUE)
eqtl_data_file <- paste0(eqtl_data_path,eqtl_data[1])
#figure_path <- "/zenodotus/dat01/mezeylab_scratch/jij2009/confeti_paper_results/Figures/"
figure_path <- "/zenodotus/dat01/mezeylab_scratch/jij2009/confeti_paper_results/test_plots/"
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

total_tests <- list()

for(single_dataset in unique.datasets){
    data_file <- grep(single_dataset, eqtl_data, value = TRUE)
    n_geno <- length(h5read(paste0(eqtl_data_path,data_file), "genotypes/col_info/id"))
    n_pheno <- length(h5read(paste0(eqtl_data_path,data_file), "phenotypes/col_info/id"))
    total_tests[[single_dataset]] <- as.numeric(n_geno) * as.numeric(n_pheno)

}

##############################################################################
########################### Import eQTL Data  ################################
##############################################################################

for(i in 1:length(unique.methods)){
    revised_method <- gsub("py", "", gsub("[0-9]", "", unique.methods[i]))
    if(revised_method %in% method_factors){

        message("Processing for method = ", unique.methods[i])
  
        method.files <- grep(paste(unique.methods[i],"_",sep=""), dataset.files, value = TRUE)

        for(f in 1:length(method.files)){
        message("Loading file = ", method.files[f])
            subdataset <- unlist(strsplit(method.files[f], "_"))[1]
            eqtl.merged.list[[revised_method]][[subdataset]] <- read_in_h5_eqtl(method.files[f])
        }


    }
      
}
cat("Read in complete\n")
##############################################################################
####################### Create Replication Plots #############################
##############################################################################

if(replication_plots){
    eqtl_rep_fraction_plots(eqtl.merged.list, max.fdr = bh_cutoff, 
                   figure.path = paste0(figure_path,"replication_fdr_plots/"), 
                   data_prefix = dataset.name, factor_sort = method_factors)

    eqtl_rep_rank_plots(eqtl.merged.list, max.rank = 10000, 
                        figure.path = paste0(figure_path,"replication_fdr_plots/"), 
                        data_prefix = dataset.name, factor_sort = method_factors)

    eqtl_rep_bonferroni_plots(eqtl.merged.list, 
                          total_tests = total_tests,
                          max.bon = 0.2,
                          figure.path = paste0(figure_path,"replication_fdr_plots/"), 
                          data_prefix = dataset.name, 
                          factor_sort = method_factors)

}

cat("Replication plots complete\n")
bh_sig_list <- lapply(eqtl.merged.list, 
                       function(x) lapply(x, function(a) subset(a, a$p.bh < bh_cutoff)))

#hit_id_list <- lapply(bh_sig_list, function(a) lapply(a, function(b) as.character(unique(b$info_id)) ) )

#intersecting_hit_ids <- lapply(hit_id_list, function(x) Reduce(intersect, x))

unique_rep_hits <- get_intersecting_hits(bh_sig_list)

rep_info_ids <- list()
intersecting_info_ids <- list()
for(single_method in names(unique_rep_hits)){
    cat("Processing method =", single_method, "\n")
    rep_info_ids[[single_method]] <- lapply(bh_sig_list[[single_method]], 
                                            function(x) subset(x, x$eqtl_id %in% unique_rep_hits[[single_method]])[,"info_id"])
    intersecting_info_ids[[single_method]] <- unique(unlist(rep_info_ids[[single_method]]))
    
}


rep_hit_list <- list()

for( single_method in names(intersecting_info_ids)){

    cat("Processing method = ", single_method, "\n")
    rep_hit_list[[single_method]] <- create_rep_hit_df(intersecting_info_ids[[single_method]])
}


##############################################################################
#################### Identifying pseudo trans hits ###########################
##############################################################################

if(pseudo_trans_plot){
    source("~/eqtl_multi_method/identify_pseudo_trans_qtl.r")
    #source("~/eqtl_multi_method/identify_pseudo_trans_qtl_blat.r")
    pseudo_cis <- lapply(rep_hit_list, function(x) detect_pseudo_trans(x, entrez2genome))

    for(i in names(pseudo_cis)){
        pseudo_cis[[i]]$eqtl_pseudo_id <- with(pseudo_cis[[i]], paste(eqtl_id, pseudo.cis.final, sep = "_")) 

    }


    rep_hits_pseudo_labeled <- pseudo_cis

    rep_hit_reduced <- lapply(rep_hits_pseudo_labeled, function(x) x[!duplicated(x$eqtl_pseudo_id),])

    ##############################################################################
    #################### Generating Significant hit map ##########################
    ##############################################################################

    for( i in names(rep_hit_reduced)){
        sig_hits_plot(rep_hit_reduced[[i]], 
                      prefix = paste(dataset.name, i, sep = "_"),
                      output_path = paste0(figure_path,"rep_sig_map/"))
    }

    ##############################################################################
    #################### Generating Bar Plots for Hits  ##########################
    ##############################################################################

    trans_hits_only <- lapply(rep_hit_reduced, function(x) subset(x, x$cis_trans == "trans"))

    #trans_plot_data <- rep_hits_plot_data(input_list = trans_hits_only, plot_factors = method_factors)
    trans_plot_data <- rep_hits_plot_data_trans(input_list = trans_hits_only, plot_factors = method_factors)

    rep_hits_barplot_trans(trans_plot_data, figure_path = paste0(figure_path,"cis_trans_rep_hits/"), 
                     prefix = paste(dataset.name,"fdr",gsub("[.]","",bh_cutoff),"trans",sep = "_"), 
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

    rep_hits_matrix_plot(trans_plot_dummy, figure_path = paste0(figure_path,"cis_trans_rep_hits/"), 
                     dataset_name = dataset.name, bh_cutoff = bh_cutoff, cis_trans = "trans")

    rep_hits_barplot(trans_plot_dummy, figure_path = paste0(figure_path,"cis_trans_rep_hits/"), 
                     prefix = paste(dataset.name,"fdr",gsub("[.]","",bh_cutoff),"trans_dummy",sep = "_"), 
                     color_set = custom_pal)


    cis_hits_only <- lapply(rep_hit_reduced, function(x) subset(x, x$cis_trans == "cis"))

    cis_plot_data <- rep_hits_plot_data(input_list = cis_hits_only, plot_factors = method_factors)

    cis_plot_data$hits <- 1
    rep_hits_matrix_plot(cis_plot_data, figure_path = paste0(figure_path,"cis_trans_rep_hits/"), 
                     dataset_name = dataset.name, bh_cutoff = bh_cutoff, cis_trans = "cis")
    #rep_hits_barplot(cis_plot_data, figure_path = paste0(figure_path,"cis_trans_rep_hits/"), 
    #                 prefix = paste(dataset.name,"fdr",gsub("[.]","",bh_cutoff),"cis",sep = "_"), 
    #                 color_set = custom_pal, no_gene_names = TRUE)


}

##############################################################################
################# #Include p-value lambda plots here  ########################
##############################################################################


##############################################################################
############ Inspect Genotype Correlation among trans hits  ##################
##############################################################################

plot_cor_list <- list()
reduced_rep_list <- lapply(rep_hit_list, function(x) x[!duplicated(x$eqtl_id),])
# check first for only trans
if(pca_geno_correlations){

    geno_data <- load_h5_data(eqtl_data_file, "genotypes")
    geno_pca <- prcomp(geno_data)
    
    for(i in names(reduced_rep_list)){
        # get genotypes for each method
        trans_geno <- subset(reduced_rep_list[[i]], cis_trans == "trans")[,"genotype"]
        # get genotype values from h5 file
        geno_subset <- geno_data[,trans_geno]
        # test correlation
        geno_trans_pc1 <- as.numeric(apply(geno_subset, 2, function(x) cor(x, geno_pca$x[,1])))
        geno_trans_pc2 <- as.numeric(apply(geno_subset, 2, function(x) cor(x, geno_pca$x[,2])))

        # get genotypes for each method
        cis_geno <- subset(reduced_rep_list[[i]], cis_trans == "cis")[,"genotype"]
        # get genotype values from h5 file
        geno_cis <- geno_data[,cis_geno]
        # test correlation
        geno_cis_pc1 <- apply(geno_cis, 2, function(x) cor(x, geno_pca$x[,1]))
        geno_cis_pc2 <- apply(geno_cis, 2, function(x) cor(x, geno_pca$x[,2]))

        plot_trans <- data.frame("cis_trans" = "trans", 
                                 "geno_pc1" = as.numeric(geno_trans_pc1), 
                                 "geno_pc2" = as.numeric(geno_trans_pc2), "method" = i)
        plot_cis <- data.frame("cis_trans" = "cis", 
                               "geno_pc1" = as.numeric(geno_cis_pc1), 
                               "geno_pc2" = as.numeric(geno_cis_pc2), "method" = i)

        plot_cor_list[[i]] <- rbind(plot_trans, plot_cis)
    }
    plot_cor <- Reduce(rbind, plot_cor_list)

    plot_title <- paste(dataset.name, "Genotype PC correlation", sep = " ")

    plot_cor$method <- factor(plot_cor$method, levels = method_factors)

    pdf(file= paste0(figure_path,"pca_correlation/",dataset.name,"_PC1_cor_combined.pdf"), width = 10, height = 4 )
        p <- ggplot(plot_cor, aes(x = geno_pc1, fill = cis_trans)) + geom_histogram(aes(y=..count..)) + xlab("Correlation with PC1") + facet_grid(.~method)
        p <- ggplot_add_theme(p) 
        p <- p + ggtitle(plot_title) + theme(strip.text = element_text( face = "bold"),
                                             strip.background = element_rect(colour = "grey", fill = "white",size = 1.5))
        print(p)
    dev.off()

    pdf(file= paste0(figure_path,"pca_correlation/",dataset.name,"_PC1_abscor_combined.pdf"), width = 10, height = 4 )
        p <- ggplot(plot_cor, aes(x = abs(geno_pc1), fill = cis_trans)) + geom_histogram(aes(y=..count..)) + xlab("Correlation with PC1") + facet_grid(.~method)
        p <- ggplot_add_theme(p) 
        p <- p + ggtitle(plot_title) + theme(strip.text = element_text( face = "bold"),
                                             strip.background = element_rect(colour = "grey", fill = "white",size = 1.5))
        print(p)
    dev.off()
}


##############################################################################
############ Inspect Genotype Correlation among trans hits  ##################
##############################################################################

if(FALSE){
    top_trans_hits <- lapply(gene_trans_count, function(x) names(x)[1])

    trans_gene_snp_list <- list()

    for( single_method in names(top_trans_hits)){
        cat("Processing method = ", single_method)
        trans_gene_snp_list[[single_method]] <- subset(trans_hits_only[[single_method]], 
            trans_hits_only[[single_method]]$pheno_symbol == top_trans_hits[[single_method]])
        
    }

    phenotypes <- load_data(eqtl_data_file, data_type = "phenotypes")
    genotypes <- load_data(eqtl_data_file, data_type = "genotypes")

    geno_cor <- function(input_genotypes, geno_mx){
        geno_subset <- geno_mx[,input_genotypes]
        geno_cor_mx <- cor(geno_subset)
        return(geno_cor_mx)
        
    }

    geno_cor_list <- lapply(trans_gene_snp_list, function(x) geno_cor(x$genotype, genotypes))
    geno_cor_vec <- lapply(geno_cor_list, function(x) abs(as.vector(x[upper.tri(x)])))


    test <- lapply(geno_cor_vec, function(x) cbind(x = seq_along(x), y = x))

    list.names <- names(test)
    lns <- sapply(test, nrow)
    dat <- as.data.frame(do.call("rbind", test))
    dat$group <- rep(list.names, lns)


    png(paste0("/zenodotus/dat01/mezeylab_scratch/jij2009/confeti_paper_results/Figures/geno_cor_plots/","test.png"))
    ggplot(dat, aes(x = group, y = y, colour = group)) +
        theme_bw() +
        geom_boxplot()
    dev.off()
        
    boxplot_cor_mx <- function(input_mx, prefix){
        png(file = paste0(prefix,"_cor_matrix.png"))
        upper_mx <- input_mx[upper.tri(input_mx)]
        boxplot(as.vector(abs(upper_mx)))
        dev.off()
        
    }
    # set plotting directory to organize plots
    for(i in names(geno_cor_list)){
        boxplot_cor_mx(geno_cor_list[[i]],
            paste0("/zenodotus/dat01/mezeylab_scratch/jij2009/confeti_paper_results/Figures/geno_cor_plots/",i))
        
    }



    # Retrieving unique hits (removing LD genotype hits)

    for(single_method in names(rep_hit_list)){
        rep_hit_list[[single_method]]$eqtl_id <- with(rep_hit_list[[single_method]],
                                                      paste(pheno_symbol, cis_trans,
                                                            pheno_chr, geno_chr, sep = "_"))
        
    }

    rep_unique <- lapply(rep_hit_list, function(x) x[!duplicated(x$eqtl_id),])

    trans_hits_unique <- lapply(rep_unique, function(x) subset(x, x$cis_trans == "trans"))

    unique_trans_count <- lapply(trans_hits_unique, function(x) sort(table(x$pheno_symbol), decreasing = TRUE))


    for( i in names(rep_unique)){
        sig_hits_plot(rep_unique[[i]], 
                      prefix = paste(dataset.name, i, "unique",sep = "_"),
                      output_path = figure_path)
    }

}

