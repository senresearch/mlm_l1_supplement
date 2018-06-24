library(qtl) # mapping quantitative trait loci

# Load cross
load("./processed/cross.ge.rda")

# Impute missing genotypes
cross.ge = fill.geno(cross.ge)
# Calculate conditional genotype probabilities
cross.ge.genoprobs = calc.genoprob(cross.ge, step=0)

# Subset the individuals with the "wet" treatment 
cross.wet = subset(cross.ge.genoprobs, 
                   ind=which(cross.ge.genoprobs$pheno$treatment=="wet"))
# Drop line, id, treatment from phenotypes
cross.wet$pheno = cross.wet$pheno[,-c(1,2,5)] 

# Calculate LOD penalties from scantwo permutations
set.seed(2)
cross.wet.perms.2dim = scantwo(cross.wet, method="hk", n.perm=100)
cross.wet.penalties = calc.penalties(cross.wet.perms.2dim)

# Function to run multiple QTL (`stepwiseqtl`) on the phenotypes, one at a time
# Only allows additive QTL models (no pairwise interactions)
# cross = a cross object
# pheno.idx = indices of phenotype column indices to analyze
# ... additional arguments to pass into `stepwiseqtl`
loop_stepwiseqtl = function(cross, pheno.idx, ...){
  # Iterate through phenotypes
  for (i in pheno.idx){
    # Run stepwiseqtl on one phenotype column
    stepwiseqtl(cross, pheno.col=i, method="hk", additive.only=TRUE, ...)
  }
}

# Sample 100 random phenotypes
set.seed(100)
pheno.idx = sample(1:ncol(cross.wet$pheno), 100)

# Runtime for 100 phenotype columns, max.qtl = 3 
# (does not count time to do scantwo permutations to get LOD penalties)
times_maxqtl3 = system.time(loop_stepwiseqtl(cross.wet, pheno.idx, max.qtl=3,
                                             penalties = cross.wet.penalties))

# Runtime for 100 phenotype columns, max.qtl = 10 (default value)
# (does not count time to do scantwo permutations to get LOD penalties)
times_maxqtl10 = system.time(loop_stepwiseqtl(cross.wet, pheno.idx, max.qtl=10, 
                                              penalties = cross.wet.penalties))

# Write times to CSV
times_out = rbind("max.qtl.3"=times_maxqtl3, "max.qtl.10"=times_maxqtl10)
write.csv(times_out, file="./processed/rqtl_100pheno_times.csv")