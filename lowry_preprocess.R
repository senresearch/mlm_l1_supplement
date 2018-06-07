library(qtl)

###############################################################################

# Attempt to pull in data from online 

# Phenotype data seems to be stored in the GEO accession page
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42408
# I downloaded the "Series Matrix File" GSE42408_series_matrix.txt
# It has the 25667 phenotypes for 208 individuals 

# Problem: trying to match this up with the genotype data. Supplemental Data 
# Set 5 from the paper, downloaded as "tpc115352SupplementalDS5.csv" has the 
# genotype matrix. But there are 341 individuals, so it's unclear how they 
# should be matched up to their phenotypes. 

# Idea: read in the raw files, subset by some type of ID to get the correct 
# 208 (104, but each has wet and dry). Write the files out again to CSV to 
# run read.cross. 

a = read.cross("csv", file="./processed2/tpc115352SupplementalDS5.csv", skip=1)
class(a)[1] = "riself"

###############################################################################

# Load cross
load("./processed/cross.ge.rda")

# Impute missing genotypes
cross.ge = fill.geno(cross.ge)
# Calculate conditional genotype probabilities
cross.ge.genoprobs = calc.genoprob(cross.ge, step=0)

# Pull out the probabilities from each chromosome and stick them all into the 
# same matrix
genoprobs = do.call(cbind, lapply(cross.ge.genoprobs$geno, 
                                  function(x){x$prob[,,1]})) 

# Rows from the same individuals aren't totally identical
# We will use the ones from the "dry" environments. 
genoprobs = genoprobs[cross.ge$pheno$treatment == "dry",]

# Set the probabilities about 0.5 to 1 and below 0.5 to 0
genoprobs[genoprobs > 0.5] = 1
genoprobs[genoprobs <= 0.5] = 0

# Check that each individual ID (rows) is consecutively duplicated in the 
# phenotypes
all(cross.ge$pheno$id[seq(1, nrow(cross.ge$pheno), by=2)] == 
      cross.ge$pheno$id[seq(2, nrow(cross.ge$pheno), by=2)])

# Split the phenotype rows into "dry" and "wet" environments and save it in 
# wide format. (Should be "dry, wet, dry, wet, ...")
pheno = as.data.frame(matrix(0, nrow(cross.ge$pheno)/2, 
                             2*ncol(cross.ge$pheno[,-(1:5)])))
pheno[,seq(1, ncol(pheno), by=2)] = cross.ge$pheno[,-(1:5)][
  cross.ge$pheno$treatment == "dry", ]
pheno[,seq(2, ncol(pheno), by=2)] = cross.ge$pheno[,-(1:5)][
  cross.ge$pheno$treatment == "wet", ]

# Rename the phenotype columns to the form "phenotype.environment", e.g. 
# "AT1G01010.dry" and "AT1G01010.wet". 
names(pheno) = paste(rep(names(cross.ge$pheno[,-(1:5)]), each=2), c("dry", "wet"), sep=".")

# Cyto contrast, to be added onto the genoprobs (X) matrix. 
cyto = cross.ge$pheno$cyto[cross.ge$pheno$treatment == "dry"]
cyto.genoprobs = cbind(cyto, genoprobs)

# Write the genotype probabilities and the cyto contrast to CSV
write.csv(cyto.genoprobs, "./processed/lowry_cyto_genoprobs.csv", row.names = FALSE)
# write.csv(genoprobs, "./processed/lowry_genoprobs.csv", row.names = FALSE)
# Write the phenotypes to CSV
write.csv(pheno, "./processed/lowry_pheno.csv", row.names = FALSE)

