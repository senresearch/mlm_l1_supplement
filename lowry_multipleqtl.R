library(qtl) # mapping quantitative trait loci

# Load cross
load("./processed/cross_ge.rda")

# Impute missing genotypes
crossGE = fill.geno(crossGE)
# Calculate conditional genotype probabilities
crossGEGenoprobs = calc.genoprob(crossGE, step=0)

# Subset the individuals with the "wet" treatment 
crossWet = subset(crossGEGenoprobs, 
                  ind=which(crossGEGenoprobs$pheno$treatment=="wet"))
# Drop code, id, cyto, and treatment from phenotypes
crossWet$pheno = crossWet$pheno[,-(1:4)]

# Calculate LOD penalties from scantwo permutations
set.seed(2)
crossWetPerms2Dim = scantwo(crossWet, method="hk", n.perm=100)
crossWetPenalties = calc.penalties(crossWetPerms2Dim)

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
pheno.idx = sample(1:ncol(crossWet$pheno), 100)

reps = 10

# Runtime for 100 phenotype columns, max.qtl = 3 
# (does not count time to do scantwo permutations to get LOD penalties)
timesMaxQTL3 = replicate(reps, system.time(
  loop_stepwiseqtl(crossWet, pheno.idx, max.qtl=3,
                   penalties = crossWetPenalties))[3])

# Runtime for 100 phenotype columns, max.qtl = 10 (default value)
# (does not count time to do scantwo permutations to get LOD penalties)
timesMaxQTL10 = replicate(reps, system.time(
  loop_stepwiseqtl(crossWet, pheno.idx, max.qtl=10,
                   penalties = crossWetPenalties))[3])

# Write times to CSV
times_out = rbind("max.qtl.3"=timesMaxQTL3, "max.qtl.10"=timesMaxQTL10)
times_out = cbind(means=rowMeans(times_out), times_out)
write.csv(times_out, file="./processed/lowry_rqtl_100pheno_times.csv")

# Estimate based on the runtime for 100 phenotypes how much time it would 
# take to run stepwiseqtl on all the phenotypes
ncol(crossWet$pheno)*2 * times_out$means/100/3600
