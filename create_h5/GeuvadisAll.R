#library(eQTLtools)

# get all genotypes for geuvadis 

# files are located in "/home/sas2030/jin/"

library(icreport)
library(rhdf5)
library(data.table)
library(bigmemory)
# get genotypes from European population
load("/zenodotus/dat01/mezeylab_scratch/jij2009/updated_input_files/GeuvEU.RData")


if(FALSE){

# of course all the genotypes were not included...

# after asking sushila AGAIN, I got the file "new_geno_info.RData"
# load("new_geno_info.RData")
# this loads geno.info
# of course this was the SAME list... so no difference was made
# now she tells me that "converted_snp_ids.RData" should have the information
load("/home/sas2030/jin/converted_snp_ids.RData")
# this loads new.names and old.names (index matched rsIDs and snpIDs respectively)

geuvadis_rsID <- as.character(snp.info$id)
geuvadis_snpID <- old.names[match(geuvadis_rsID, new.names)]


all_geno_snpID = fread("/home/sas2030/jin/genotypes.cols",
                             sep="auto", sep2="auto", header=FALSE, na.strings="NA",
                             stringsAsFactors=FALSE, verbose=FALSE)


idx_in_geno_mx <- match(geuvadis_snpID, all_geno_snpID[,get("V1")])

geno_mx <- attach.big.matrix("/home/sas2030/jin/genotypes.desc")

# create a copy of the genotype big matrix with the selected genotypes
geno_subset <- deepcopy(geno_mx, cols = idx_in_geno_mx)
geno_subset_mx <- as.matrix(geno_subset)
colnames(geno_subset_mx) <- geuvadis_rsID

geno_no_na <- t(na.omit(t(geno_subset_mx)))
}
# there seems to be a lot of NA values in the dataset. need to figure out how to deal with those. 

# for now I'll just omit any entries with NA which brings the genotype count down to 831016 from 845099




# Load expression data + covariate matrix


Geuvadis.660.expr.dt = fread("~/Datasets/Geuvadis/GD660.GeneQuantRPKM.txt",
                             sep="auto", sep2="auto", nrows=-1L, header="auto", na.strings="NA",
                             stringsAsFactors=FALSE, verbose=FALSE)
setkey(Geuvadis.660.expr.dt, Chr)
Geuvadis.660.expr.mx = as.matrix(Geuvadis.660.expr.dt[,!c("TargetID","Gene_Symbol","Chr","Coord"), with=FALSE])
rownames(Geuvadis.660.expr.mx) = Geuvadis.660.expr.dt[,Gene_Symbol]



short.QC.462 = read.table("~/Datasets/Geuvadis/E-GEUV-1.sdrf.txt", 
						  sep = "\t", header = TRUE, stringsAsFactors = FALSE)
filenames.462 = as.character(unique(short.QC.462$Assay.Name))
# Extended QC file
Geuvadis.667.info.dt = fread("~/Datasets/Geuvadis/GD667.QCstats.masterfile.txt",
                             sep="\t", nrows=-1L,na.strings="NA",header = TRUE,
                             stringsAsFactors=FALSE)

setkey(Geuvadis.667.info.dt, Sample_ID)

Geuvadis.462.info.dt = Geuvadis.667.info.dt[Geuvadis.667.info.dt$FILE %in% filenames.462,
                                            c("FILE","Sample_ID","Population","SeqLabNumber","RNAConcentration_ng.ul","LibraryPrepPlate","LibraryConcentrationMethod",
                                              "LibraryConcentration_ng.ul","INSERT_SIZE_MODE","Mapped","Multiple_mapping","Total_read","Unique_mapping"),with = FALSE]


filenames.462 <- as.character(filenames.462)

# Extract samples used in analysis by file name
Geuvadis.462.RPKM.mx = Geuvadis.660.expr.mx[,filenames.462]
# Set column names to id
colnames(Geuvadis.462.RPKM.mx) <- substr(colnames(Geuvadis.462.RPKM.mx),1,7)


setkey(Geuvadis.462.info.dt,Sample_ID)
colnames(Geuvadis.462.info.dt)

# Gender Processing with 1000G ped file
genome1000.ped <- read.table("~/Datasets/Geuvadis/1000G_sampleinfo.ped", 
						     sep = "\t", header = TRUE, stringsAsFactors = F)
sex.info <- genome1000.ped[,c(2,5)]
rownames(sex.info) <- sex.info$Individual.ID
sex.info$Individual.ID <- NULL


# 1 is MALE 2 is FEMALE

Geuvadis.462.info.dt[,Sex := sex.info[Geuvadis.462.info.dt$Sample_ID,]]
Geuvadis.462.info.dt$Sex <- sub(1,"MALE",Geuvadis.462.info.dt$Sex)
Geuvadis.462.info.dt$Sex <- sub(2,"FEMALE",Geuvadis.462.info.dt$Sex)

phenotype.mx <- as.matrix(Geuvadis.462.RPKM.mx)

info.df <- as.data.frame(Geuvadis.462.info.dt)
rownames(info.df) <- Geuvadis.462.info.dt[,get("Sample_ID")]
# save a version of geuvadis data with all samples in one place. 


if(FALSE){

all_pop_file <- "/zenodotus/dat01/mezeylab_scratch/jij2009/eqtl_input_h5/GeuvadisALL.h5"
all_pop_pheno <- t(Geuvadis.462.RPKM.mx[gene.info$probe,])
# all_pop_geno <- geno_no_na[Geuvadis.462.info.dt[,get("Sample_ID")],]

all_snp_info <- snp.info[snp.info$id %in% colnames(geno_no_na),]
# Also create datasets by spliting populations
gene.info = gene.info[order(gene.info$chrom.num, gene.info$start, decreasing = FALSE),]

create_h5_file(all_pop_file)

h5_add_data(file_name = all_pop_file,
			input_matrix = all_pop_pheno,
			data_type = "phenotypes")

h5_add_data(file_name = all_pop_file,
			input_matrix = all_pop_geno,
			data_type = "genotypes")

h5_add_data(file_name = all_pop_file,
			input_matrix = as.matrix(info.df),
			data_type = "covars")

h5_add_snp_info(file_name = all_pop_file,
				snp_id = all_snp_info$id,
				snp_chr = as.character(all_snp_info$chromosome),
				snp_pos = all_snp_info$position)

h5_add_pheno_info(file_name = all_pop_file,
					pheno_id = gene.info$probe,
					pheno_chr = as.character(gene.info$chromosome),
					pheno_start = gene.info$start,
					pheno_end = gene.info$end,
					pheno_entrez = gene.info$entrez,
					pheno_symbol = gene.info$symbol)

H5close()

unique_pop <- unique(info.df$Population)
for(single_pop in unique_pop){

	#RData_file <- "GTExAdiposeSubq.RData"

	file_name <- paste0("/zenodotus/dat01/mezeylab_scratch/jij2009/eqtl_input_h5/Geuvadis",single_pop,".h5")

	sample_id_pop <- info.df[info.df$Population == single_pop, "Sample_ID"]
	# first sort the genes by the chromosome
	phenotypes = all_pop_pheno[sample_id_pop,]
	genotypes = all_pop_geno[sample_id_pop,]
	covars = as.matrix(info.df[info.df$Population == single_pop,])
	cat("Number of samples in ", single_pop, "=", length(sample_id_pop),"\n")
	cat("phenotypes dimension = ", dim(phenotypes)[1], "x", dim(phenotypes)[2],"\n")
	cat("genotypes dimension = ", dim(genotypes)[1], "x", dim(genotypes)[2],"\n")
	create_h5_file(file_name)

	h5_add_data(file_name = file_name,
	            input_matrix = phenotypes,
	            data_type = "phenotypes")

	h5_add_data(file_name = file_name,
	            input_matrix = genotypes,
	            data_type = "genotypes")

	h5_add_data(file_name = file_name,
	            input_matrix = covars,
	            data_type = "covars")


	h5_add_snp_info(file_name = file_name,
	                snp_id = all_snp_info$id,
	                snp_chr = as.character(all_snp_info$chromosome),
	                snp_pos = all_snp_info$position)

	h5_add_pheno_info(file_name = file_name,
	                  pheno_id = gene.info$probe,
	                  pheno_chr = as.character(gene.info$chromosome),
	                  pheno_start = gene.info$start,
	                  pheno_end = gene.info$end,
	                  pheno_entrez = gene.info$entrez,
	                  pheno_symbol = gene.info$symbol)

	H5close()
	rm(genotypes, phenotypes, covars)

}
}