library(rhdf5)
library(ggplot2)
result_directory <- "/zenodotus/dat01/mezeylab_scratch/jij2009/eqtl_multi_method/BH_sig_hits_by_dataset/"
h5_directory <- "/zenodotus/dat01/mezeylab_scratch/jij2009/eqtl_h5_with_Kmx/"
directory_list <- dir(result_directory)

python_list <- grep("pyresult", directory_list, value = TRUE)

R_list <- grep("Rresult", directory_list, value = TRUE)

file_info_list <- c()

for(single_dir in python_list){
    message("Processing = ", single_dir)
    current_path <- paste(result_directory, single_dir, "/", sep = "")
    file_list <- grep("[.]txt", dir(current_path), value = TRUE)
    bh_file_list <- grep("_bh_", file_list, value = TRUE)
    for(single_file in bh_file_list){
        
        file_in_process <- paste0(current_path, single_file)
        
        single_file_info <- file.info(file_in_process)
        size_formated <- utils:::format.object_size(single_file_info["size"], "MB")
        message("File = ", single_file, " / size = ", size_formated)
        
        parsed_file <- unlist(strsplit(single_file, "_"))
        
        data_set <- parsed_file[1]
        method <- parsed_file[2]
        pval_correction <- parsed_file[3]
        cutoff <- parsed_file[4]
        
        # Get phenotype / genotype information from h5 file
        
        h5_file <- paste0(h5_directory,data_set, ".h5")
        h5_info <- h5ls(h5_file)
        n_pheno <- as.numeric(h5_info[(h5_info$group == "/phenotypes/col_info" & 
        h5_info$name == "id"), "dim"])
        n_geno <- as.numeric(h5_info[(h5_info$group == "/genotypes/col_info" & 
        h5_info$name == "id"), "dim"])
        n_sample <- as.numeric(h5_info[(h5_info$group == "/phenotypes/row_info" & 
        h5_info$name == "id"), "dim"])
        bon_pval <- 0.05 / (n_geno * n_pheno)
        
        pipe_file <- pipe(paste0("cut -f3 -d' ' ", file_in_process), 'r') 
        raw_pvals <- read.table(pipe_file, header = TRUE)
        bon_sig_hits <- sum(raw_pvals < bon_pval)
        close(pipe_file)
        #bon_pval = 0
        #bon_sig_hits = 0
        temp <- data.frame("dataset" = data_set, 
                   "method" = method, 
                   "pval_method" = pval_correction, 
                   "bh_cutoff" = cutoff, 
                   "bon_cutoff" = as.numeric(bon_pval), 
                   "bon_sig_hits" = as.numeric(bon_sig_hits),
                   "size" = size_formated, 
                   "size_raw" = as.numeric(single_file_info["size"]),
                   "filename" = single_file, row.names = NULL)
        file_info_list <- rbind(file_info_list, temp)
        rm(raw_pvals)
    }    
}

#python_info_df = data.frame(Reduce(rbind, lapply(file_info_list, function(x) as.character(x))), row.names = NULL)

#colnames(python_info_df) <- c("dataset" , "method" , "pval_method" , "bh_cutoff", 
#                              "bon_cutoff", "size","size_raw","filename")

python_info_df <- file_info_list
write.csv(python_info_df, file = "~/eqtl_multi_method/Python_analysis_table.csv", row.names = FALSE)

comment_out = TRUE
data_groups <- c("Geuvadis", "GTEx", "HapMap", "NHSAE", "Smith2008")
if(comment_out){

pdf(file = paste0("Test","_test_dataset_summary.pdf"), width = 12, height = 5)
for(single_dataset in data_groups){
    cat("Processing datset = ", single_dataset, "\n")
    subset_info <- python_info_df[grep(single_dataset, python_info_df$dataset),]
    subset_info$size_raw <- as.numeric(as.character(subset_info$size_raw))
    subset_info$size_mb <- subset_info$size_raw / 10^6
    subset_info$log_size <- log10(subset_info$size_raw)
    unique_datsets <- length(unique(as.character(subset_info$dataset)))
    print(ggplot(subset_info, aes( x = method, y = size_mb, fill = method))+ geom_bar(stat = 'identity') + 
    facet_grid(.~dataset) + scale_fill_brewer(palette="Set1"))
    
}
dev.off()

pdf(file = paste0("Test","_bon_dataset_summary.pdf"), width = 12, height = 5)
for(single_dataset in data_groups){
    cat("Processing datset = ", single_dataset, "\n")
    subset_info <- python_info_df[grep(single_dataset, python_info_df$dataset),]
#    subset_info$size_raw <- as.numeric(as.character(subset_info$size_raw))
#    subset_info$size_mb <- subset_info$size_raw / 10^6
#    subset_info$log_size <- log10(subset_info$size_raw)
#    unique_datsets <- length(unique(as.character(subset_info$dataset)))
    print(ggplot(subset_info, aes( x = method, y = bon_sig_hits, fill = method))+ geom_bar(stat = 'identity') + 
    facet_grid(.~dataset) + scale_fill_brewer(palette="Set1"))
    
}
dev.off()

}
# Dataset - Method - BH cutoff - file size - file name 

#single_file_info <- file.info(input_file)

#size_formated <- utils:::format.object_size(single_file_info["size"], "auto")
