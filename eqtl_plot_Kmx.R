library(rhdf5)
library(NMF)
library(ggplot2)

input_h5 <- as.character(commandArgs(TRUE)[1])

message("Processing file = ", input_h5, "\n")

# get files
prefix <- gsub("[.]h5", "", input_h5)
K_mx_list <- h5read(input_h5, "/K_mx")


for (i in 1:length(K_mx_list)){
	
	method <- names(K_mx_list)[i]
	message("Creating file for method = ", method, "\n")
	single_Kmx <- K_mx_list[[method]]
	png(file = paste(prefix, "_" ,method, "_Kmx.png", sep = ""), width = 900, height = 800)
	NMF::aheatmap(single_Kmx, Rowv = NA, Colv = NA, main = paste(prefix,"_", method, "_Kmx", sep = ""))
	dev.off()
	
}

