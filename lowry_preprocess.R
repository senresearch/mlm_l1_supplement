library(qtl) # mapping quantitative trait loci

###############################################################################

# Attempt to pull in data from online 

# Phenotype data seems to be stored in the GEO accession page
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42408
# I downloaded the "Series Matrix File" GSE42408_series_matrix.txt
# It has the 25662 phenotypes for 208 individuals 

# The files I used has 208 individuals (104, but each has wet and dry) and 
# 51324/2 = 25662 phenotypes. There are 450 markers (plus the cyto contrasts)

# Problem: trying to match this up with the genotype data. Supplemental Data 
# Set 5 from the paper, downloaded as "tpc115352SupplementalDS5.csv" has the 
# genotype matrix. There are 341 individuals, 104 of which match up with the 
# individuals (use "line marker" or "line" to match) in cross.ge. Problem is 
# that there are only 168 markers, 123 are matches for the 450 in cross.ge. 

# Idea: read in the raw files, subset by some type of ID to get the correct 
# 208 (104, but each has wet and dry). Write the files out again to CSV to 
# run read.cross. 

a = read.cross("csv", file="./processed2/tpc115352SupplementalDS5.csv", skip=1)
class(a)[1] = "riself"

library(data.table)
b = fread("./processed2/GSE42408_series_matrix.txt", skip=72, header=TRUE)

# 104 matching RILs 
length(intersect(a$pheno$Line.Marker, cross.ge$pheno$line))

# 123 matching markers
sum(sapply(1:5, function(i){
  length(intersect(colnames(cross.ge$geno[[i]]$data), colnames(a$geno[[i]]$data)))
}))


d = read.cross("csvs", 
               genfile="./processed2/TPC2015-00122-RAR1_Supplemental_Data_set_1b.csv", 
               phefile="./processed2/TPC2015-00122-RAR1_Supplemental_Data_set_1a.csv", 
               genotypes=c("a","b"), skip=1)
b = fread("./processed2/TPC2015-00122-RAR1_Supplemental_Data_set_1b.csv", skip=1)
sum(sapply(1:3, function(i){
  length(intersect(colnames(cross.ge$geno[[i]]$data), colnames(d$geno[[i]]$data)))
}))

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

# Rows from the same individuals are similar but not totally identical
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
names(pheno) = paste(rep(names(cross.ge$pheno[,-(1:5)]), each=2), 
                     c("dry", "wet"), sep=".")

# Cyto contrast, to be added onto the genoprobs (X) matrix. 
cyto = cross.ge$pheno$cyto[cross.ge$pheno$treatment == "dry"]
cyto.genoprobs = cbind(cyto, genoprobs)

# Write the genotype probabilities and the cyto contrast to CSV
write.csv(cyto.genoprobs, "./processed/lowry_cyto_genoprobs.csv", 
          row.names = FALSE)
# Write the phenotypes to CSV
write.csv(pheno, "./processed/lowry_pheno.csv", row.names = FALSE)

