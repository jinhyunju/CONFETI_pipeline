# python script for saving values below a given threshold in a large matrix
import h5py
import numpy as np
import sys, getopt
import re
import time
import sys
import time

def get_options(argv):
    inputfile = ''
    try:
        opts, args = getopt.getopt(argv, "i:t:", ["input=", "threshold="])
    except getopt.GetoptError:
        print "Incorrect input for script \n"
        print "-i <input_pval_h5> -t <threshold>\n"
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-i", "--input"):
            input_pval_file = arg
        elif opt in ("-t", "--threshold"):
            threshold = float(arg)

    print "Input pval file = ", input_pval_file
    print "Threshold = ", threshold
    return input_pval_file, threshold


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
testing = "off"
direction = "geno"
localtime = time.asctime( time.localtime(time.time()) )
print "P-val Rank Script with containers Started :", localtime

if testing == "on":
    pval_h5_file = "SmithTest_SCALEICA_pvals.h5"
    threshold = 0.3
else :
    pval_h5_file, threshold = get_options(sys.argv[1:])

print "Running script pval_rank.py"
print "Loading p-values saved in %s" % pval_h5_file

prefix = re.sub("\.h5", '', pval_h5_file)

open_h5 = h5py.File(pval_h5_file, 'r')

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


if direction == "geno":
    print "Chunks in hdf5 accessed by genotypes"
    chunk_idx = create_idx(N_geno, geno_chunk_size)
    n_chunks = len(chunk_idx)
    pval_sl = []
    print "Loading %d batches of p-values based on hdf5 chunk size" % n_chunks

    chunk_pvals = np.empty((N_pheno, geno_chunk_size), dtype = np.float32)

    for i in range(0,n_chunks):
        print "Reading Chunk %d / %d " % ( i + 1, n_chunks)
        subset_idx = chunk_idx[i]
        if i == (n_chunks -1):
            chunk_pvals = pval_mx[:, range(subset_idx[0], subset_idx[1])]
        else :
            pval_mx.read_direct(chunk_pvals, np.s_[:,subset_idx[0]:subset_idx[1]])
        pval_sl.append(np.extract(chunk_pvals < threshold, chunk_pvals))
else :
    print "Chunks in hdf5 accessed by phenotypes"
    chunk_idx = create_idx(N_pheno, pheno_chunk_size)
    n_chunks = len(chunk_idx)
    pval_sl = []
    for i in range(0,n_chunks):
        subset_idx = chunk_idx[i]
        chunk_pvals = pval_mx[range(subset_idx[0], subset_idx[1]), : ] 
        pvals_to_save = np.extract(chunk_pvals < threshold, chunk_pvals)
        pval_sl.append(pvals_to_save)
        sys.stdout.write(".")
        sys.stdout.flush()


# declare empty list to append p-values to
localtime = time.asctime( time.localtime(time.time()) )
print "P-value read in completed starting sorting :", localtime

pval_sl = [item for sublist in pval_sl for item in sublist]
pval_sl.sort()

total_saved = len(pval_sl)

print "Total of %d p-values smaller than threshold = %f" % (total_saved, threshold)

CurrentMin = 1
indices = []
pval_list = []
adjpval_list = []

k = total_saved


while k > 0 :
    pval_index = k - 1
    pval = pval_sl[pval_index]
    adj_pval = pval * N_total_pval / k
    if adj_pval < CurrentMin:
        if adj_pval <= threshold:
            indices.append(pval_index)
            pval_list.append(pval)
            adjpval_list.append(adj_pval)
    CurrentMin = adj_pval
    k -= 1


if len(adjpval_list) > 0:
    print "BH p-values stopped at rank %d with adj pval %f" % (len(pval_list), max(adjpval_list))
else :
    print "No BH p-values were below threshold"
h5py.File.close(open_h5)

# save output into hdf5 format
# parallelize first for loop to make it fast
output_file = prefix + "_bh_" + str(threshold) + ".h5"
output_h5 = h5py.File(output_file, 'w')

output_h5.create_dataset("rank", data = indices)
output_h5.create_dataset("pval", data = pval_list)
output_h5.create_dataset("p.bh", data = adjpval_list)
h5py.File.close(output_h5)


localtime = time.asctime( time.localtime(time.time()) )
print "P-val Rank Script with containers ended :", localtime


