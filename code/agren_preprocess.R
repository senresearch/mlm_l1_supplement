library(qtl) # mapping quantitative trait loci

# Read in the data as a cross object
# Downloaded from https://datadryad.org/resource/doi:10.5061/dryad.77971
agren = read.cross("csvs",dir="../data", 
               genfile="geno.csv", 
               phefile= "RIL_DataForSelectionAnalyses3yrs2.csv",
               genotypes=c("a","b"))
class(agren)[1] = "riself"

# Only keep IDs and phenotypes for fruits per seedling
agren$pheno = agren$pheno[,c(1,5:10)]
agren = subset(agren, ind=!apply(is.na(agren$pheno[,-1]), 1, any))

# Impute missing genotypes
agren = fill.geno(agren)

# Write the genotypes and phenotypes to CSVs
write.cross(agren, "csvs", "../processed/agren")

# Calculate conditional genotype probabilities
agrenGenoprobs = calc.genoprob(agren, step=1) 

# Pull out the probabilities from each chromosome and put them all into the 
# same matrix
genoprobs = do.call(cbind, lapply(agrenGenoprobs$geno, 
                                  function(x){x$prob[,,1]})) 

# Set the probabilities about 0.5 to 1 and below 0.5 to 0
genoprobs[genoprobs > 0.5] = 1
genoprobs[genoprobs <= 0.5] = 0

# Write the genotype probabilities to CSV
write.csv(genoprobs, "../processed/agren_genoprobs.csv", row.names = FALSE)
