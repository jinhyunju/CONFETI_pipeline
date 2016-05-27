#!/usr/bin/python

import sys, getopt    # to get command line arguments
import limix.deprecated.modules.panama as PANAMA # get panama from limix
import limix.deprecated.modules.qtl as qtl
import numpy as np    # for creating arrays and saving csv files
import re             # regular expression for string modification
import os             # retrieving current path and modifying directories 
import scipy as sp
import h5py
import gc
from multiprocessing import Pool, Process, Value, Array
from sklearn import preprocessing

# function for getting command line arguments 
def get_options(argv):
    inputfile = ''
    try:
        opts, args = getopt.getopt(argv, "i:m:c:b:",["input=","--method", "--cores", "--batches"])
    except getopt.GetoptError:
        print 'Incorrect input for options \n'
        print '-i <input data> -m <Method> -p <n.cores> -b <batch_numbers>'
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-i", "--input"):
            inputfile = arg
        elif opt in ("-m", "--method"):
            method = arg
        elif opt in ("-c", "--cores"):
            cores = int(arg)
        elif opt in ("-b", "--batches"):
            N_batch = int(arg)

    print "Input file =", inputfile
    print "Method for eQTL = ", method
    print "Number of cores = ", cores
    print "Number of Batches = ", N_batch
    return inputfile, method, cores, N_batch


# function for running the linear mixed model with multiprocessing
def worker(i):
    global pheno, input_geno, K_mx
    #print "Processing phenotype # %d" % (i +1)
    single_pheno = pheno[:,i]
    lmm_result = qtl.test_lmm(snps = input_geno, pheno = single_pheno, K = K_mx, covs = None, test = 'lrt', verbose = False)
    single_pvals = np.array(lmm_result.getPv())
    return single_pvals

def create_idx(n_genotype, total_batches):
    div_interval = range(0,n_genotype, n_genotype / total_batches)
    batch_idx = []
    for i in range(0, total_batches):
        if (i == (total_batches - 1)):
            batch_idx.append([div_interval[i], n_genotype])
        else :
            batch_idx.append([div_interval[i], div_interval[i+1]])    
    return batch_idx


if __name__ == '__main__':
    
    file_name, method, cores, N_batch = get_options(sys.argv[1:])
    # Input file location
    
    prefix = re.sub('\.h5', '',file_name)

    if "LINEAR" in method:
        method = "LINEAR"
    elif "ICE" in method:
        method = "ICE"
    
    # The script will be run in $TMPDIR so file should be created in current directory
    pval_path = './'
    pval_file = prefix + "_" + method + "py_pvals.h5"

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
    geno = np.transpose(geno_read)


    N_pheno = pheno.shape[1]
    N_geno = geno.shape[1]
    N_sample = geno.shape[0]
    ##############################################################################
    ####################### calculate K matrix based on method ###################
    ##############################################################################
    print "Confounding Factor Method = %s \n" % method
    print "Loading %s similarity matrix \n" % method
    
    Kmx_location = "/K_mx/" + method
    K_mx = np.array(h5_file[Kmx_location]).astype('float64')

    h5py.File.close(h5_file)

    Kmx_file = prefix + "_" + method + "_Kmx.h5"
    Kmxh5 = h5py.File(Kmx_file, 'w')
    Kmx_save = Kmxh5.create_dataset("Kmx", (N_sample,N_sample), dtype = 'float32', chunks = True)
    Kmx_save[:,:] = K_mx
    h5py.File.close(Kmxh5)
    
    ##############################################################################
    ######################### Generate Batch Indexes #############################
    ##############################################################################
    parallel_h5 = h5py.File(pval_file, 'w')

    pvals = parallel_h5.create_dataset("pvals", (N_pheno, N_geno), dtype = 'float32', chunks = True)

    batch_idx = create_idx(N_geno, N_batch)

    for b in range(0,N_batch):
        print "Processing Batch Number = %d \n" % (b+1)
        subset_idx = batch_idx[b]
        input_geno = geno[:, range(subset_idx[0], subset_idx[1])]
        N_subset_geno = input_geno.shape[1]  
        ##############################################################################
        ############################### eQTL Part ####################################
        ##############################################################################
        p = Pool(cores)
        parallel_pvals = p.map(worker, range(0,N_pheno))
    
        p.close()
        p.join() 
        
        print "Saving P-values for batch %d \n" % (b+1)
        parallel_pvals = np.array(parallel_pvals).reshape(N_pheno,N_subset_geno)
        pvals[:,range(subset_idx[0], subset_idx[1])] = parallel_pvals		
        parallel_pvals = None
        gc.collect()

    h5py.File.close(parallel_h5)
    
    



