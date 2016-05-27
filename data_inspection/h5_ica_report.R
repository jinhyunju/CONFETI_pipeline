
library(lrgpr)

input_path <- "/zenodotus/dat01/mezeylab_scratch/jij2009/eqtl_h5_with_Kmx/"

file_list <- grep("Geuvadis", dir(input_path), value = TRUE)



for(i in file_list){
	single_file <- i
	prefix <- gsub("[.]h5", "",single_file)
	h5file <- paste0(input_path, single_file)

	ica_result <- h5_ica(h5_file = h5file,var.cutoff = 95)
	#ica_result <- ica_covar_check(ica_list = ica_result,h5_file = h5file)
	ica_result <- ica_genotype_association(ica_list = ica_result, h5_file = h5file)
	ica_result <- get_gene_info(ica_result, h5file)

	ica_report(ica_result, h5_file = h5file, prefix = prefix, output.path = "/zenodotus/dat01/mezeylab_scratch/jij2009/eqtl_reports/")



}
 