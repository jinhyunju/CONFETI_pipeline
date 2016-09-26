# python script for saving values below a given threshold in a large matrix
import h5py
import numpy as np
import sys, getopt
import re
import time
import sys
import time

def get_options(argv):
    input_pval_file = ''
    try:
        opts, args = getopt.getopt(argv, "i:", ["input="])
    except getopt.GetoptError:
        print "Incorrect input for script \n"
        print "-i <input_pval_h5>\n"
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-i", "--input"):
            input_pval_file = arg

    print "Input pval file = ", input_pval_file
    return input_pval_file

def create_idx(n_genotype, chunk_size):
    div_interval = range(0, n_genotype, chunk_size)
    batch_idx = []
    total_batches = len(div_interval)
    for i in range(0, total_batches):
        if (i == (total_batches - 1 )):
            batch_idx.append([div_interval[i], n_genotype])
        else :
            batch_idx.append([div_interval[i], div_interval[i+1]])
    return batch_idx


#input pvalue matrix & threshold
localtime = time.asctime( time.localtime(time.time()) )
print "P-val Median Script Started :", localtime

pval_h5_file = get_options(sys.argv[1:])

#pval_h5_file = "Smith2008_PARTICApy_pvals.h5"

print "Running script pvals_median.py"
print "Loading p-values saved in %s" % pval_h5_file

prefix = re.sub("\.h5", '', pval_h5_file)

open_h5 = h5py.File(pval_h5_file, 'a')

pval_mx = open_h5["/pvals"] 

N_geno = pval_mx.shape[1]

N_pheno = pval_mx.shape[0]

N_total_pval = N_pheno * N_geno


if pval_mx.chunks is None:
    geno_chunk_size = 100
    pheno_chunk_size = 100
else :
    geno_chunk_size = pval_mx.chunks[1]
    pheno_chunk_size = pval_mx.chunks[0]


median_h5_file = prefix + "_median.h5"

print "Chunks in hdf5 accessed by phenotypes"
chunk_idx = create_idx(N_pheno, pheno_chunk_size)
n_chunks = len(chunk_idx)
pval_median = np.empty([N_pheno])
for i in range(0,n_chunks):
    subset_idx = chunk_idx[i]
    chunk_pvals = pval_mx[range(subset_idx[0], subset_idx[1]), : ] 
    chunk_medians = np.median(chunk_pvals, axis = 1)
    pval_median[range(subset_idx[0], subset_idx[1])] = chunk_medians
    sys.stdout.write(".")
    sys.stdout.flush()

h5py.File.close(open_h5)
# declare empty list to append p-values to
localtime = time.asctime( time.localtime(time.time()) )
print "P-value median extraction completed : ", localtime

total_medians = pval_median.shape[0]

print "Total %d medians extracted for %d phenotypes" % (total_medians, N_pheno)

save_median = h5py.File(median_h5_file, 'a')
save_median.create_dataset("median_pvals", data = pval_median, chunks = True)


localtime = time.asctime( time.localtime(time.time()) )
print "P-val median Script completed :", localtime


