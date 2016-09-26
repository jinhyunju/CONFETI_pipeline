options(warn=-1)
suppressMessages(require(picaplot))

suppressMessages(require(rhdf5))
options(warn=0)


file_name <- as.character(commandArgs(TRUE)[1])
prefix <- gsub("[.]h5", "", file_name)
output_dir <- "/zenodotus/dat01/mezeylab_scratch/jij2009/eqtl_genotype_pca/"
input_path <- "/zenodotus/dat01/mezeylab_scratch/jij2009/eqtl_input_h5/"


dir.create(file.path(output_dir, prefix), showWarnings = FALSE)
figure_output <- paste0(output_dir, prefix, "/")


input_file <- paste0(input_path, file_name)

message("Processing file = ", input_file, "\n")

geno <- h5read(input_file, "genotypes/matrix")
colnames(geno) <- h5read(input_file, "genotypes/col_info/id")
rownames(geno) <- h5read(input_file, "genotypes/row_info/id")


message("Creating Genotype PCA matrix \n")
geno_pca <- run_pca(t(geno))


png(file = paste0(figure_output, prefix, "_gpca_summary.png"), 
    width = 2400, height = 800, res = 150)
    summary_plot = plot_summary(geno_pca, plot_here = TRUE)
dev.off()


for (i in 1:10){

    png(file = paste0(figure_output, prefix, "_gpca_",i,".png"), 
        width = 900, height = 400)

    p = plot_component(geno_pca, comp_idx = i)

    dev.off()


}
