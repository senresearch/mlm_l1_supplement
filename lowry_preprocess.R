library(qtl) # mapping quantitative trait loci
library(data.table) # quickly read in tables

# Table that includes IDs and cytoplasm for markers
ids = read.csv("./processed/dd2014_cytocovar.csv") 
# Create indicator for cyto
ids$cyto = ifelse(ids$cyto.num==1, 0, 1)

# Genotype matrix, downloaded from Supplemental Table 1 in
# http://www.plantcell.org/content/27/4/969/tab-figures-data
# (needs to be updated to include all the markers)
geno = fread("./processed/TKrils_map55.csv")

# Series matrix file, downloaded from 
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42408
# Contains phenotypes and some other information about the experiment
SMF = fread("./processed/GSE42408_series_matrix.txt", header=FALSE, fill=TRUE)
# Pull out the sample title column, which includes the line markers and 
# treatments
sample_title = strsplit(as.character(SMF[37,-1]), " ")
# Build phenotype matrix
pheno = data.frame(sapply(sample_title, function(x){sub(",", "", x[2])}), 
                   # marker line
                   sapply(sample_title, function(x){x[3]}), 
                   # treatment
                   apply(t(SMF[-c(1:73,nrow(SMF)),-1]), 2, as.numeric))
names(pheno) = c("code", "treatment", 
                 as.character(SMF[-c(1:73,nrow(SMF))]$V1))

# Merge phenotype matrix with ID matrix 
pheno = merge(ids[c("id", "code", "cyto")], pheno, by="code", all.y=TRUE)

# Merge the phenotypes and genotypes to make one big matrix
out = merge(pheno, geno, by="id", all.x=TRUE)

# Make the matrix ready to be read into R/qtl
blanks = matrix(NA, 2, ncol(pheno))
colnames(blanks) = colnames(pheno)
out = rbind(cbind(blanks, geno[1:2,-1]), out)
write.csv(out, "./processed/lowry_raw.csv", row.names = FALSE, na="")

# Read in the data as a cross object
cross.ge = read.cross("csv", file="./processed/lowry_raw.csv", 
                      genotypes=c("a","b"))
class(cross.ge)[1] = "riself"
save(cross.ge, file="./processed/cross.ge.rda")

###############################################################################

if(!exists("cross.ge")){
  load("./processed/cross.ge.rda")
}

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
                             2*ncol(cross.ge$pheno[,-(1:4)])))
pheno[,seq(1, ncol(pheno), by=2)] = cross.ge$pheno[,-(1:4)][
  cross.ge$pheno$treatment == "dry", ]
pheno[,seq(2, ncol(pheno), by=2)] = cross.ge$pheno[,-(1:4)][
  cross.ge$pheno$treatment == "wet", ]

# Rename the phenotype columns to the form "phenotype.environment", e.g. 
# "AT1G01010.dry" and "AT1G01010.wet". 
names(pheno) = paste(rep(names(cross.ge$pheno[,-(1:4)]), each=2), 
                     c("dry", "wet"), sep=".")

# Cyto contrast, to be added onto the genoprobs (X) matrix. 
cyto = cross.ge$pheno$cyto[cross.ge$pheno$treatment == "dry"]
cyto.genoprobs = cbind(cyto, genoprobs)

# Write the genotype probabilities and the cyto contrast to CSV
write.csv(cyto.genoprobs, "./processed/lowry_cyto_genoprobs.csv", 
          row.names = FALSE)
# Write the phenotypes to CSV
write.csv(pheno, "./processed/lowry_pheno.csv", row.names = FALSE)
