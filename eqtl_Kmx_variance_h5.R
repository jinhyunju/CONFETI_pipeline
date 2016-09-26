options(warn=-1)
suppressMessages(library(eQTLtools))
suppressMessages(library(icreport, lib.loc = "/home/jij2009/R_packages/"))
suppressMessages(require(lrgpr))
suppressMessages(require(rhdf5))
options(warn=0)

as.character.formula <- function(x){
   Reduce( paste, deparse(x) )  
}

input_file <- as.character(commandArgs(TRUE)[1])
method <- as.character(commandArgs(TRUE)[2])
var_data <- as.numeric(commandArgs(TRUE)[3])


message("Processing file = ", input_file, "\n")

pheno <- load_h5_data(input_file, "phenotypes")
geno <- load_h5_data(input_file, "genotypes")

if(method %in% c("CONFETI", "CONPANA")){
	method <- "CONFETI"
    message("Creating ", method ," similarity matrix \n")
	message("Correction method = ",method,"\n")
	message("Running ICA\n")
	group_to_save <- paste0(method, var_data)
	h5createGroup(input_file, group_to_save)
    pca.pheno <- prcomp(pheno)
    percent <- (cumsum(pca.pheno$sdev^2) /sum(pca.pheno$sdev^2)) * 100
    k.est <- which(percent > var_data)[1]
	
	ICA.obj <- h5_ica(h5_file = input_file, k.est = k.est)

	ICA.obj <- ica_genotype_association(ica_list = ICA.obj, h5_file = input_file)
	n.covars <- length(ICA.obj$non_genetic)
	
	message(n.covars, " out of ", k.est, " ICs used as confounding factors \n")

	if(n.covars == 0){
		stop("No ICs estimated as hidden factors")
	}

	non_genetic_factors <- match(ICA.obj$non_genetic, colnames(ICA.obj$S))
	recon.mx <- ICA.obj$S[,non_genetic_factors] %*% ICA.obj$A[non_genetic_factors,]
	similarity.mx <- cov(recon.mx)
	
	
    h5write(similarity.mx, input_file, paste0("K_mx/", group_to_save))
	h5write(ICA.obj$S, input_file, paste0(group_to_save,"/S"))
	h5write(ICA.obj$A, input_file, paste0(group_to_save,"/A"))
	h5write(ICA.obj$non_genetic, input_file, paste0(group_to_save,"/hf"))
	h5write(recon.mx, input_file, paste0(group_to_save,"/pheno_mx"))
	rm(ICA.obj, similarity.mx, recon.mx, n.covars)
	gc()
	
		
}

if(method == "ICAAMX"){

    message("Creating ICAAMX similarity matrix \n")
	message("Correction method = ICAAMX \n")
	message("Running ICA\n")

	group_to_save <- paste0(method, var_data)
	h5createGroup(input_file, group_to_save)
    pca.pheno <- prcomp(pheno)
    percent <- (cumsum(pca.pheno$sdev^2) /sum(pca.pheno$sdev^2)) * 100
    k.est <- which(percent > var_data)[1]
	
	ICA.obj <- h5_ica(h5_file = input_file, k.est = k.est)

	ICA.obj <- ica_genotype_association(ica_list = ICA.obj, h5_file = input_file)
	n.covars <- length(ICA.obj$non_genetic)
	
	message(n.covars, " out of ", k.est, " ICs used as confounding factors \n")

	if(n.covars == 0){
		stop("No ICs estimated as hidden factors")
	}

	non_genetic_factors <- match(ICA.obj$non_genetic, colnames(ICA.obj$S))

	weighted_comps <- t(ICA.obj$A[non_genetic_factors,]) %*% diag(ICA.obj$percent_var[non_genetic_factors])
	similarity.mx <- cov(t(weighted_comps))
	
    h5write(similarity.mx, input_file, paste0("K_mx/", group_to_save))
	h5write(ICA.obj$S, input_file, paste0(group_to_save,"/S"))
	h5write(ICA.obj$A, input_file, paste0(group_to_save,"/A"))
	h5write(ICA.obj$non_genetic, input_file, paste0(group_to_save,"/hf"))
	
	rm(ICA.obj, similarity.mx,  n.covars)
	gc()
	
		
}


if(method == "PCALMM"){

    message("Correction method = PCALMM\n")
    message("Running PCA\n")
	group_to_save <- paste0(method, var_data)

    pca.pheno <- gene_expr_pca(phenotype.mx = t(pheno))
    k.est <- which(cumsum(pca.pheno$var.percent) > var_data)[1]
    message(k.est," PCs explained more than ", var_data," percent of the variance \n")
   
    pca.pheno$hf <- pca_genotype_test(pca.pheno, geno, n.cores = 1)$hf

    hf_idx <- match(pca.pheno$hf, colnames(pca.pheno$x))
    hf_idx <- hf_idx[hf_idx <= k.est]

    n.covars <- length(hf_idx)

    message(n.covars, " out of ",k.est," PCs used as confounding factors \n")

    weighted.pc.proj <- pca.pheno$x[,hf_idx] %*% diag(pca.pheno$var.percent[hf_idx])

    Kmx_PCA <- cov(t(weighted.pc.proj))
#	    Kmx_PCA <- Kmx_PCA / sum(diag(Kmx_PCA))
    h5write(Kmx_PCA, input_file, paste0("K_mx/",group_to_save))
}

if(method == "PCAKMX"){

    message("Correction method = PCALMM\n")
    message("Running PCA\n")
	group_to_save <- paste0(method, var_data)

    pca.pheno <- gene_expr_pca(phenotype.mx = t(pheno))
   
    pca.pheno$hf <- pca_genotype_test(pca.pheno, geno, n.cores = 1)$hf

    hf_idx <- match(pca.pheno$hf, colnames(pca.pheno$x))

    n.covars <- length(hf_idx)

    message(n.covars, " out of ",length(pca.pheno$var.percent)," PCs used as confounding factors \n")

    weighted.pc.proj <- pca.pheno$x[,hf_idx] %*% diag(pca.pheno$var.percent[hf_idx])

    Kmx_PCA <- cov(t(weighted.pc.proj))
	Kmx_PCA <- Kmx_PCA / sum(diag(Kmx_PCA))
    h5write(Kmx_PCA, input_file, paste0("K_mx/",group_to_save))
}


if(method == "ICE"){

    message("Creating ICE similarity matrix \n")
	message("Generating sample similarity matrix based on expression values\n")

	message("Normalizing expression values")

	phenotypes.norm <- apply(pheno, 2, function(x) (x - mean(x)) / sd(x) )
	expr.cov <- cov(t(phenotypes.norm))

    h5write(expr.cov, input_file, "K_mx/ICE")
}

if(method == "LINEAR"){

    message("Creating LINEAR similarity matrix \n")
	message("Generating Identity matrix\n")

	linear_Kmx <- diag(rep(1, nrow(pheno)))

    h5write(linear_Kmx, input_file, "K_mx/LINEAR")
}

H5close()
