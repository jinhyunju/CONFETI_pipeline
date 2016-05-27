## Scripts that are needed for running the CONFETI eQTL pipeline


File | Description | Example
------------ | ------------- | ------------
`qsub_python_xqtl.sh` | Submiting eQTL job to cluster | `qsub -N GTExAdiSLinear -t 1-6 qsub_python_xqtl.sh -F "-i=GTExAdiposeSubq.h5" -F "-b=120"`                    
`eqtl_hits_summary.R` | Generate summary figures for eQTL analysis | `Rscript eqtl_hits_summary.R "GTExAdipose"`
`eqtl_Kmx_variance.R` | Calculate LMM similarity matrix for various methods | `Rscript eqtl_Kmx_variance.R <inputfile> <method> <variance to use>` 
`PANAMA_mx_var.R` | Calculate LMM similarity matrix for PANAMA | `python PANAMA_mx_var.py -i $INPUTFILE -m $METHOD -v $VAR` 