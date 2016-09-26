# Splitting Geuvadis populations including YRI

# Load expression data including all 462 samples

load("~/Datasets/Geuvadis/Geuvadis_462_ICA_input.RData")
# loads info.mx, phenotype.mx

# Load gene information and snp information 

load("/zenodotus/dat01/mezeylab_scratch/jij2009/updated_input_files/GeuvEU.RData")
# loads genotypes, phenotypes, snp.info, gene.info

yri_geno <- read.table("~/Datasets/Geuvadis/YRI.traw", header = TRUE)

yri_geno$CHR <- NULL
yri_geno$X.C.M <- NULL
yri_geno$POS <- NULL
yri_geno$COUNTED <- NULL
yri_geno$ALT <- NULL

rownames(yri_geno) <- as.character(yri_geno$SNP)

yri_geno$SNP <- NULL

colnames(yri_geno) <- sapply(strsplit(colnames(yri_geno),"_"),function(x) x[2])

# Load genotype matrix 