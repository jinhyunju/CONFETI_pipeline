import sys, getopt    # to get command line arguments
import limix.deprecated.modules.panama as PANAMA # get panama from limix
import limix.deprecated.modules.qtl as qtl
import numpy as np    # for creating arrays and saving csv files
import re             # regular expression for string modification
import os             # retrieving current path and modifying directories 
import scipy as sp
import h5py
import gc
from sklearn import preprocessing

# function for getting command line arguments 
def get_options(argv):
    inputfile = ''
    try:
        opts, args = getopt.getopt(argv, "i:m:v:",["input=","method=","var="])
    except getopt.GetoptError:
        print 'Incorrect input for options \n'
        print '-i <input data> -m <method>'
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-i", "--input"):
            inputfile = arg
        elif opt in ("-m", "--method"):
            method = arg
        elif opt in ("-v", "--var"):
            var_data = arg

    print "Input file =", inputfile
    print "Method = ", method
    print "Variance = ", var_data
    return inputfile, method, var_data



file_name, method,var_data  = get_options(sys.argv[1:])

# Input file location

# open hdf5 file for processing 
h5_file = h5py.File(file_name,'r+')

##############################################################################
################## read in genotype and phenotype data #######################
##############################################################################
# read in genotype and phenotype data
pheno_read = np.array(h5_file['phenotypes/matrix'])
geno_read = np.array(h5_file['genotypes/matrix'])

print "Centering Phenotypes \n"
pheno_mean = np.mean(pheno_read, 1)
pheno_centered = pheno_read - pheno_mean[:, np.newaxis]

pheno = np.transpose(pheno_centered)

pheno_raw = np.transpose(pheno_read)
geno = np.transpose(geno_read)
genetic_similarity = np.cov(geno)

N_pheno = pheno.shape[1]
N_geno = geno.shape[1]
N_sample = geno.shape[0]


if method == "PANAMA":
    print "Running SVD to get number of factors\n"
    U, s, V = np.linalg.svd(pheno_centered, full_matrices= False)
    var = s ** 2
    percent_var = var / sum(var)
    cumsum_var = np.cumsum(percent_var)
    # find number of components needed to explain more than 90% of variance
    k_est = [ n for n,i in enumerate(cumsum_var) if i> (float(var_data) / float(100)) ][0] + 1
    print "%d factors used to train PANAMA model (explaining more than %f percent of variance)\n" % (k_est,float(var_data))

    print "No pre calculated PANAMA matrix detected, running PANAMA\n"
    panama_model = PANAMA.PANAMA(Y = pheno, Kpop = genetic_similarity)
    panama_model.train(k_est)
    Kpanama = panama_model.get_Kpanama()
    Ktot = panama_model.get_Ktot()
    vcomps = panama_model.get_varianceComps()

    group_to_save = method + str(var_data)
    print "Saving PANAMA matrix to input hdf5 file \n"
    Kmx_group = str("/K_mx/" + group_to_save)

    save_PANAMA_mx = h5_file.create_dataset(Kmx_group, (N_sample, N_sample), dtype = 'float64')
    save_PANAMA_mx[:,:] = Kpanama

    PANAMA_info = h5_file.create_group(group_to_save)
    h5_file.create_dataset(str(group_to_save + "/vcomp"), data = np.array(vcomps.values()))
    h5_file.create_dataset(str(group_to_save + "/names"), data = np.array(vcomps.keys()))
else :
    print "> Skipping PANAMA calculation, matrix already exists \n"


if method == "CONPANA":
    print "No pre calculated CONPANA matrix detected, running CONPANA\n"
    
    print "Importing CONFETI recon matrix for PANAMA calculation\n"  
    group_to_save = method + str(var_data)
    Kmx_group = str("/K_mx/" + group_to_save)
    group_to_load = "CONFETI" + str(var_data) + "/pheno_mx"
    ph_recon_mx = h5_file[group_to_load]
    CONFETI_recon_mx = ph_recon_mx[:]
    print "Running SVD to get number of factors\n"
    U, s, V = np.linalg.svd(CONFETI_recon_mx, full_matrices= False)
    var = s ** 2
    percent_var = var / sum(var)
    cumsum_var = np.cumsum(percent_var)
    # find number of components needed to explain more than 90% of variance
    k_est = [ n for n,i in enumerate(cumsum_var) if i> (float(var_data) / float(100)) ][0] + 1
    print "%d factors used to train PANAMA model\n" % k_est
       
    part_hybrid_model = PANAMA.PANAMA(Y = CONFETI_recon_mx, Kpop = genetic_similarity)
    part_hybrid_model.train(k_est)
    Kpanama_ph = part_hybrid_model.get_Kpanama()
    Ktot_ph = part_hybrid_model.get_Ktot()
    vcomps_ph = part_hybrid_model.get_varianceComps()
    
    print "Saving PANAMA matrix to input hdf5 file \n"
    save_PANAMA_ph = h5_file.create_dataset(Kmx_group, (N_sample, N_sample), dtype = 'float64')
    save_PANAMA_ph[:,:] = Kpanama_ph
    
    PANAMA_info = h5_file.create_group(group_to_save)
    h5_file.create_dataset(str(group_to_save + "/vcomp"), data = np.array(vcomps_ph.values()))
    h5_file.create_dataset(str(group_to_save + "/names"), data = np.array(vcomps_ph.keys()))
else :
    print "> Skipping CONPANA matrix calculation, matrix already exists \n"

h5py.File.close(h5_file)
