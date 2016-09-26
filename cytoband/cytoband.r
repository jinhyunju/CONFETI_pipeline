
load('chrom_bands.RData')


get.cband <- function(chrom, position) {
  if ( is.na(position) ) {
    return ( chrom )
  } else {
   return ( paste(chrom, chrom.bands$name[which.max(chrom.bands$chrom.num==chrom & position < chrom.bands$chromEnd)], sep='') )
  }
}


test3 = apply(test, 1, function(x) get.cband(x["geno_chr"], x["geno_pos"]))



# get unique genotypes

geno_only <- plyrr::ldply(lapply(eqtl.merged.list, function(x) 
                          plyr::ldply(lapply(x, function(a) 
                          subset(a, select = c(genotype, geno_chr, geno_pos))))))

geno_only <- geno_only[!duplicated(geno_only$genotype),]

geno_only$chrom_band <- apply(geno_only, 1, function(x) get.cband(x["geno_chr"], x["geno_pos"]))

test <- eqtl.merged.list
for(m in names(eqtl.merged.list)){
    for (d in names(eqtl.merged.list[[m]])){

      eqtl.merged.list[[m]][[d]]$cytoband <- geno_only[match(eqtl.merged.list[[m]][[d]]$genotype, geno_only$genotype), "chrom_band"]
      eqtl.merged.list[[m]][[d]]$eqtl_id <- with(eqtl.merged.list[[m]][[d]], paste(pheno_symbol, cis_trans, cytoband)
    }


}

# Load map position data

maps <- list()

for ( i in 1:22 ) {

  load(paste0('genetic_maps/chr_', i, '.RData'))

  maps[[i]] <- chrom.map  

}



# Create functions for chromosomes with map position data

map.funs <- lapply(maps, function (m) {

  approxfun(x=m$new.pos, y=m$cm)

})



# Dummy function for other chromosomes

map.funs$X <- map.funs$Y <- map.funs$MT <- function (x) x/1e6



# Main function which operates on vectors


map.pos <- function(chrom, position) {

  stopifnot(length(chrom)==length(position))

  

  sapply(1:length(position), function (i) {

    cm <- map.funs[[chrom[i]]]

    cm(position)

  })

}