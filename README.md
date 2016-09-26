## Scripts that are needed for running the CONFETI eQTL pipeline


File | Description | Example
------------ | ------------- | ------------
`qsub_python_xqtl.sh` | Submiting eQTL job to cluster | `qsub -N GTExAdiSLinear -t 1-6 qsub_python_xqtl.sh -F "-i=GTExAdiposeSubq.h5" -F "-b=120"`                    
`eqtl_hits_summary.R` | Generate summary figures for eQTL analysis | `Rscript eqtl_hits_summary.R "GTExAdipose"`
`eqtl_Kmx_variance.R` | Calculate LMM similarity matrix for various methods | `Rscript eqtl_Kmx_variance.R <inputfile> <method> <variance to use>` 
`PANAMA_mx_var.R` | Calculate LMM similarity matrix for PANAMA | `python PANAMA_mx_var.py -i $INPUTFILE -m $METHOD -v $VAR` 


### Process outline for running an eQTL analysis

1) Create hdf5 file based on RData objects or in memory matrices.

Need phentype, genotype, covariate information, gene information, SNP information

- Template script location: `/home/jij2009/eqtl_multi_method/create_h5/`


2) submit job script using 

`/zenodotus/dat01/mezeylab_store/jij2009/eqtl_qsub_scripts/qsub_python_xqtl.sh`

Example

`qsub -N NHSAEM00 -t 1-6 -pe smp 16 -l "h_vmem=20G" qsub_python_xqtl.sh -F "-i=NHSAEM00.h5" -F "-b=100" -F "-v=95"`

