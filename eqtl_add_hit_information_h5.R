library(eQTLtools) # custom function library
library(icreport)
library(gtools)
library(rhdf5)

# modifying directory list to current directory for qsub script running in $TMPDIR

cis_distance <- as.numeric(commandArgs(TRUE)[1])

input_file <- as.character(commandArgs(TRUE)[2])

output_dir <- as.character(commandArgs(TRUE)[3])


# Getting directories which hold eQTL_results

#eqtlmerged.output.path <- paste("/zenodotus/dat01/mezeylab_scratch/jij2009/eqtl_multi_method/", eqtl.merged.directory,"/", sep = "")


message("Reading Results from xqtl python output \n")

file.list <- dir("./")

# Read in eQTL results

cat("Loading Benjamini-Hochberg corrected p-values \n")
pval.file <- grep("_hits[.]txt", file.list, value = TRUE)

# change column names

gene.info <- as.data.frame(h5read(input_file, "phenotypes/col_info"))
colnames(gene.info)[grep("id", colnames(gene.info))] <- "pheno_id"
if(is.null(gene.info$pheno_symbol)){
		# In case the gene symbol is unknown (non-human dataset)
		gene.info$pheno_symbol <- gene.info$pheno_id
}



snp.info <- as.data.frame(h5read(input_file, "genotypes/col_info"))
colnames(snp.info)[grep("id", colnames(snp.info))] <- "geno_id"


file_data_info <- unlist(strsplit(pval.file,"_"))

merged_h5_file <- paste(file_data_info[1],"_" ,file_data_info[2],"_eqtldf.h5", sep = "")

merged_h5_file <- paste(output_dir, merged_h5_file, sep = "")
#check if merged file already exists

eqtl_hits <- read.table(pval.file, header = TRUE, stringsAsFactors = FALSE)


geno_names <- as.character(h5read(input_file, "genotypes/col_info/id"))
eqtl_hits$genotype <- geno_names[eqtl_hits$genotype.idx]

pheno_names <- as.character(h5read(input_file, "phenotypes/col_info/id"))
eqtl_hits$phenotype <- pheno_names[eqtl_hits$phenotype.idx]


H5close()

eqtl_hits <- na.omit(eqtl_hits)

clean.eqtl <- eQTLtools::cleaning_results(eqtl_hits, snp.info, gene.info, cis_distance)
#clean.eqtl$geno_chr <- as.character(clean.eqtl$geno_chr)
clean.eqtl$eqtl_id <- paste(clean.eqtl$pheno_symbol, clean.eqtl$cis_trans, as.character(clean.eqtl$geno_chr), sep ="_")

#clean.eqtl$geno_chr <- factor(clean.eqtl$geno_chr, levels = mixedsort(levels(clean.eqtl$geno_chr)))

message("Creating ", merged_h5_file, "\n") 

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

#h5createFile(merged_h5_file)

#h5write(clean.eqtl, merged_h5_file, "/eqtl_df")

#with(clean.eqtl, cbind(pheno_symbol, cis_trans, geno_chr, as.character(geno_chr), eqtl_id))



