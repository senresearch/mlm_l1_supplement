library(qtl) # mapping quantitative trait loci

###############################################################################

# Attempt to pull in data from online 

# Downloaded here: https://datadryad.org/resource/doi:10.5061/dryad.77971
# geno.csv and pheno.csv are formatted for R/qtl: 
# 544 individuals, 348 markers, 17 phenotypes 
# (6 QN, 6 not, 4 from ED paper, 1 ID)

# The CSV file I used originally has 400 individuals, 348 markers, and 7 
# phenotypes (6 phenotypes, 1 ID)

# Problem: 
# 1. phenotypes in CSV file I used are on a different scale compared to 
# pheno.csv (both QN and not QN). 
# 2. Unclear how to subset individuals to match up. While it is clear 
# that you can reduce 544 to around 400 by dropping individuals with 
# missing phenotype information for the non-ED columns, it's not clear 
# how the exactly 400 were chosen. 

a = read.cross("csvs",dir="./processed2",
               genotypes=c("a","b"), 
               genfile="geno.csv", 
               phefile= "pheno.csv", 
               na.strings=c("NA","-"))
class(a)[1] = "riself"

sum(!apply(is.na(a$pheno[,11:17]), 1, any)) # 390, compared to 392 in original CSV
length(intersect(a$pheno[,17], agren$pheno[,7])) # 399 intersect

###############################################################################

# Read in the data as a cross object
agren = read.cross("csv",dir="./processed", genotypes=c("a","b"), 
                   file="agren2013_fullGenoMatrix4qtl_withParents.csv", 
                   na.strings=c("NA","-"))

# Impute missing genotypes
agren = fill.geno(agren)

# Write the genotypes and phenotypes to CSVs
write.cross(agren, "csvs", "./processed/agren")

# Calculate conditional genotype probabilities
agren.genoprob = calc.genoprob(agren, step=1) 

# Pull out the probabilities from each chromosome and put them all into the 
# same matrix
genoprobs = do.call(cbind, lapply(agren.genoprob$geno, 
                                  function(x){x$prob[,,1]})) 

# Set the probabilities about 0.5 to 1 and below 0.5 to 0
genoprobs[genoprobs > 0.5] = 1
genoprobs[genoprobs <= 0.5] = 0

# Write the genotype probabilities to CSV
write.csv(genoprobs, "./processed/agren_genoprobs.csv", row.names = FALSE)
