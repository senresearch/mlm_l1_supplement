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
# Drop code, id, cyto, and treatment from phenotypes
cross.wet$pheno = cross.wet$pheno[,-(1:4)]

# Calculate LOD penalties from scantwo permutations
set.seed(2)
cross.wet.perms.2dim = scantwo(cross.wet, method="hk", n.perm=100)
cross.wet.penalties = calc.penalties(cross.wet.perms.2dim)

# Function to run multiple QTL analysis (`stepwiseqtl`) on individual phenotypes
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

reps = 10

# Runtime for 100 phenotype columns, max.qtl = 3 
# (does not count time to do scantwo permutations to get LOD penalties)
times_maxqtl3 = replicate(reps, system.time(
  loop_stepwiseqtl(cross.wet, pheno.idx, max.qtl=3,
                   penalties = cross.wet.penalties))[3])

# Runtime for 100 phenotype columns, max.qtl = 10 (default value)
# (does not count time to do scantwo permutations to get LOD penalties)
times_maxqtl10 = replicate(reps, system.time(
  loop_stepwiseqtl(cross.wet, pheno.idx, max.qtl=10,
                   penalties = cross.wet.penalties))[3])

# Write times to CSV
times_out = rbind("max.qtl.3"=times_maxqtl3, "max.qtl.10"=times_maxqtl10)
times_out = cbind(means=rowMeans(times_out), times_out)
write.csv(times_out, file="./processed/lowry_rqtl_100pheno_times.csv")

# Estimate based on the runtime for 100 phenotypes how much time it would 
# take to run stepwiseqtl on all the phenotypes
ncol(cross.wet$pheno)*2 * times_out$means/100/3600
