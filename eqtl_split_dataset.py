import h5py
import numpy as np
import sys, getopt
import os
import re

def get_options(argv):
    inputfile = ''
    try:
        opts, args = getopt.getopt(argv, "i:b:",["input=", "batches="])
    except getopt.GetoptError:
        print 'Incorrect input for options \n'
        print '-i <input data> -b <batch_numbers>'
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-i", "--input"):
            inputfile = arg
        elif opt in ("-b", "--batches"):
            N_batch = int(arg)

    print "Input file =", inputfile
    print "Number of Batches = ", N_batch
    return inputfile,  N_batch


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

    file_name,  N_batch = get_options(sys.argv[1:])

    prefix = re.sub('\.h5', '' ,file_name)

    h5_file = h5py.File(file_name, 'a')

    geno_mx = np.array(h5_file['genotypes/matrix'])

    n_sample = geno_mx.shape[1]
    n_geno = geno_mx.shape[0]

    batch_idx = np.array(create_idx(n_geno, N_batch))

    save_batch_idx = h5_file.create_dataset("batch_idx",(batch_idx.shape[0],batch_idx.shape[1]),data = batch_idx ,dtype= 'float32', chunks = True)

    h5py.File.close(h5_file)