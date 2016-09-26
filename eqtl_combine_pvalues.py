import h5py
import numpy as np
import sys, getopt
import os
import re

def get_options(argv):
    inputfile = ''
    try:
        opts, args = getopt.getopt(argv, "i:p:o:",["input=", "pvalpath=", "outputpath="])
    except getopt.GetoptError:
        print 'Incorrect input for options \n'
        print '-i <input data> -p <pvalue_path> -o <output_path>'
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-i", "--input"):
            inputfile = str(arg)
        elif opt in ("-p", "--pvalpath"):
            pvalue_path = str(arg)
        elif opt in ("-o", "--outputpath"):
            output_path = str(arg)

    print "Input file =", inputfile
    print "Path to p-value files = ", pvalue_path
    print "Output path to save combined p-value file = ", output_path
    return inputfile, pvalue_path, output_path

if __name__ == '__main__':

    file_name,  pvalue_path, output_path = get_options(sys.argv[1:])
    # find all files in p-value directory

    file_list = os.listdir(pvalue_path)
    file_list.sort()

    # get original dimensions 
    h5_file = h5py.File(file_name, 'a')

    geno_mx = h5_file['genotypes/matrix']
    pheno_mx = h5_file['phenotypes/matrix']

    n_geno = geno_mx.shape[0]
    n_pheno = pheno_mx.shape[0]
    h5py.File.close(h5_file)

    # create master p_value matrix
    pval_prefix = file_list[0].split("_pvals_", 1)[0]

    pval_all_file_name = output_path + pval_prefix + "_pvals.h5"

    pval_combine = h5py.File(pval_all_file_name, 'a')
    all_pvals = pval_combine.create_dataset("pvals",(n_pheno,n_geno),dtype= 'float64', chunks = True)

    #saved_geno = 0
    # combine all files 
    for single_file in file_list:
        single_file_path = pvalue_path + single_file
        pval_subset_file = h5py.File(single_file_path, "r")
        pval_subset_mx = np.array(pval_subset_file["pvals"])
        subset_batch_idx = np.array(pval_subset_file["batch_idx"])
        #geno_count = pval_subset_mx.shape[1]
        #all_pvals[:,range(saved_geno, saved_geno + geno_count)] = pval_subset_mx
        print "Processing file =", single_file
        print "P-value indexes = [%d, %d]" %(subset_batch_idx[0], subset_batch_idx[1])
        all_pvals[:,range(int(subset_batch_idx[0]), int(subset_batch_idx[1]))] = pval_subset_mx
        #saved_geno = saved_geno + geno_count
        h5py.File.close(pval_subset_file)

    # Close script
    h5py.File.close(pval_combine)