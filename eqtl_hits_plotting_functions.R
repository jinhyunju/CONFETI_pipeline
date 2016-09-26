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
#    fixed_hit <- regmatches(temp.df$hit_id, regexpr("_", temp.df$hit_id), invert = TRUE)
#    temp.df$fixed_hit_id <- sapply(fixed_hit, function(x) paste(x[1], x[2], sep = "::"))
#    temp.df$info_id <- with(temp.df, paste(fixed_hit_id,geno_chr, geno_pos, pheno_chr, pheno_start, pheno_symbol,cis_trans, sep = "::"))

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
    pdf(paste(figure_path, prefix,"_trans_hit_by_method_label.pdf", sep = ""), height = 10, width = 30)
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
                temp_df <- data.frame("hits" = 1, "pseudo_hits" = 0, "gene" = rownames(count_table_list[[i]]))
            } else {
                temp_df <- data.frame("hits" = 0, "pseudo_hits" = 1, "gene" = rownames(count_table_list[[i]]))
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
    pdf(paste(figure_path, prefix,"_trans_hit_by_method_count.pdf", sep = ""), height = 10, width = 20)
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
                                 color_set = custom_pal, dataset_name, cis_trans){

    plot_matrix <- convert_data_to_matrix(plot_data)
    png(paste0(figure_path, dataset_name,"_replicating_",cis_trans,"_hits_overlap.png"), width = 1000, height = 250)
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

eqtl_rep_fraction_plots <- function(eqtl.merged.list, 
                                   max.fdr = 0.1, 
                                   figure.path, 
                                   plot_prefix, 
                                   factor_sort = NULL){

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
            ct_count_df <- data.frame(melt(sapply(intersecting_hits, cis_trans_counter)))
            ct_count_df$FDR <- FDR_range[c]
            colnames(ct_count_df) <- c("cis_trans", "method", "count", "FDR")
            cis_trans_list[[c]] <- ct_count_df
    }

    cis_trans_df <- Reduce(rbind, lapply(cis_trans_list, function(x) x))
  
    cis_trans_cast <- dcast(cis_trans_df, method + FDR ~ cis_trans, value.var = "count")
    if(!is.null(factor_sort)){
        cis_trans_df$method <- factor(cis_trans_df$method, levels = factor_sort)
    }

    cis_trans_rep_plot <- ggplot(cis_trans_df, aes(x = FDR, y = count, col = method, group = method)) + facet_grid(cis_trans~.,scales = "free_y") +
                  geom_line(size = 1)  + scale_color_manual(values = custom_pal)
                  
    cis_trans_rep_plot <- ggplot_add_theme(cis_trans_rep_plot)

    pdf(file = paste(figure.path, plot_prefix, "_replication_FDR.pdf", sep = ""), height = 8, width = 5)
        cis_trans_rep_plot <- cis_trans_rep_plot+ theme(legend.text = element_text(size = 10), 
                                                        legend.key.size = unit(1, "cm"),
                                                        strip.background = element_rect(colour = "white", 
                                                                                        fill = "white",
                                                                                        size = 0.8))
        print(cis_trans_rep_plot)
    dev.off()


}



eqtl_rep_bonferroni_plots <- function(eqtl.merged.list, 
                                      total_tests = NULL,
                                      max.bon = 0.2,
                                      figure.path, 
                                      plot_prefix, 
                                      factor_sort = NULL){
    cat("Generate bonferroni rep plots\n")

    max.pval.data <- sapply(eqtl.merged.list, function(a) sapply(a, function(b) max(range(b$pval), na.rm = TRUE)))
    max.pval.data <- max(unlist(max.pval.data), na.rm = TRUE)

    if(max.pval.data < max.bon){

        max.bon <- max.pval.data

    }

    bon_range <- seq(0, max.bon, length.out = 100)
#    bon_range <- bon_range[-1]

    cis_trans_list <- list()

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
            ct_count_df <- data.frame(melt(sapply(intersecting_hits, cis_trans_counter)))
            ct_count_df$bon <- bon_range[c]
            colnames(ct_count_df) <- c("cis_trans", "method", "count", "bon")
            cis_trans_list[[c]] <- ct_count_df
    }

    cis_trans_df <- Reduce(rbind, lapply(cis_trans_list, function(x) x))
 

    cis_trans_rep_plot <- ggplot(cis_trans_df, aes(x = bon, y = count, col = method, group = method)) + facet_grid(cis_trans~.,scales = "free_y") +
                  geom_line(size = 1)  + scale_color_manual(values = custom_pal)
                  
    cis_trans_rep_plot <- ggplot_add_theme(cis_trans_rep_plot)

    pdf(file = paste(figure.path, plot_prefix, "_replication_bonferroni.pdf", sep = ""), height = 8, width = 5)
        cis_trans_rep_plot <- cis_trans_rep_plot+ theme(legend.text = element_text(size = 10), 
                                                        legend.key.size = unit(1, "cm"),
                                                        strip.background = element_rect(colour = "white", 
                                                                                        fill = "white",
                                                                                        size = 0.8))
        print(cis_trans_rep_plot)
    dev.off()



}


eqtl_rep_rank_plots <- function(eqtl.merged.list, 
                               max.rank = 10000, 
                               figure.path, 
                               plot_prefix, 
                               factor_sort = NULL){

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
    hit.range <- round(seq(1, max.rank, length.out = 100))

    message("Creating Rank replication plot \n")
    intersecting_hits <- lapply(hit.range, function(x) find_intersecting_hits(eqtl.merged.list, x))

    cis_trans_list <- list()

    for(i in 1:length(hit.range)){
        ct_count_df <- data.frame(melt(sapply(intersecting_hits[[i]], cis_trans_counter)))
        ct_count_df$rank <- hit.range[i]
        colnames(ct_count_df) <- c("cis_trans", "method", "count", "rank")
        cis_trans_list[[i]] <- ct_count_df
            
    }

    cis_trans_df <- Reduce(rbind, lapply(cis_trans_list, function(x) x))
 

    cis_trans_rep_plot <- ggplot(cis_trans_df, aes(x = rank, y = count, col = method, group = method)) + 
                  facet_grid(cis_trans~.,scales = "free_y") +
                  geom_line(size = 1)  + scale_color_manual(values = custom_pal)
                  
    cis_trans_rep_plot <- ggplot_add_theme(cis_trans_rep_plot)

    pdf(file = paste(figure.path, plot_prefix, "_replication_rank.pdf", sep = ""), height = 8, width = 5)
        cis_trans_rep_plot <- cis_trans_rep_plot+ theme(legend.text = element_text(size = 10), 
                                                        legend.key.size = unit(1, "cm"),
                                                        strip.background = element_rect(colour = "white", 
                                                                                        fill = "white",
                                                                                        size = 0.8))
        print(cis_trans_rep_plot)
    dev.off()
}



# Function for creating significant hits plot
sig_hits_map <- function(hits.info, prefix = "test", output_path = "./",
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
    hits.info$shape <- ifelse(hits.info$cis_trans == "trans" & hits.info$pseudo.cis.final, 4, 1)
    hits.info$col[hits.info$cis_trans == "trans" & hits.info$pseudo.cis.final] <- "red"

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