# Load libraries
library(reshape2) # reshape data frames
library(MESS) # AUC

# Read in unique condition-concentration names
Xnames = lapply(1:6, function(i){
  temp = as.character(read.csv(paste0("../processed/processed_KEIO_data/p", i, 
                                      "_krit_cond_conc_names.csv"), 
                               header=FALSE)[,1])
  out = unique(temp)
  out = out[-length(out)]
  return(out)
})
# Read in unique mutant names
Znames = lapply(1:6, function(i){
  temp = as.character(read.csv(paste0("../processed/processed_KEIO_data/p", i, 
                                      "_krit_mut_names.csv"), 
                               header=FALSE)[,1])
  out = unique(temp)
  out = out[-length(out)]
  return(out)
})


# Read in least squares MLM t-statistics
tStats = lapply(1:6, function(i){
  temp = as.matrix(read.csv(paste("../processed/nichols_p", i, "_tStats.csv", 
                                  sep=""), sep=",", header=FALSE))
  return(temp[-nrow(temp), -ncol(temp)])
})
# Use condition-concentrations and mutants for the interaction names
for (i in 1:6) {
  rownames(tStats[[i]]) = Xnames[[i]]
  colnames(tStats[[i]]) = Znames[[i]]
}


# Pull out coefficients corresponding to lambda of interest
lambda = 10 # MSE is minimized here + 1 SE (lambda = 0.455166)
# Read in L1 coefficient estimates
# Use condition-concentrations and mutants for the interaction names
l1Coeffs = lapply(1:6, function(i){
  temp = read.csv(paste("../processed/nichols_p", i, "_l1_coeffs.csv", sep=""), 
                  sep=",", header=FALSE)
  out = matrix(temp[,lambda], 
               length(Xnames[[i]])+1, length(Znames[[i]])+1)[-1,-1]
  rownames(out) = Xnames[[i]]
  colnames(out) = Znames[[i]]
  return(out)
})


# Indices of minimal media conditions for each plate
minimalIdx = sapply(Xnames, function(x){grep("M9min", x)})
# Auxotrophs are mutants with no growth in minimial media
# "No growth" status for a mutant is determined by whether the the quantile at 
# quantCutoff for the coefficients corresponding to minimal media conditions 
# is negative
# Quantile cutoff for negative coefficients
quantCutoff = 0.95


# Least squares t-statistics
# Range of median cutoffs used to determine auxotroph status
cutoffsl2 = seq(-30, 30, by=0.5)
# Pull out the mutants whose t-statistics corresponding to minimial media 
# conditions that have medians below the cutoffs
mlmAuxol2 = lapply(cutoffsl2, function(cutoff){
  names(do.call(c, sapply(lapply(1:6, function(i){
    apply(tStats[[i]][minimalIdx[[i]],], 2, median)}), function(x){
      which(x < cutoff)})))
})
# Get the TRUE/FALSE labels for auxotrophs (above/below cutoffs)
mlmLabelsl2 = sapply(mlmAuxol2, function(y){
  do.call(c, sapply(1:6, function(i){
    sapply(Znames[[i]], function(x){x %in% y})}))})


# L1-penalized estimates
# Range of median cutoffs used to determine auxotroph status
cutoffsl1 = seq(-3, 3, by=0.05) 
# Pull out the mutants whose coefficients corresponding to minimial media 
# conditions that have medians below the cutoffs
mlmAuxol1 = lapply(cutoffsl1, function(cutoff){
  names(do.call(c, sapply(lapply(1:6, function(i){
    apply(l1Coeffs[[i]][minimalIdx[[i]],], 2, median)}), function(x){
      which(x < cutoff)})))
})
# Get the TRUE/FALSE labels for auxotrophs (above/below cutoffs)
mlmLabelsl1 = sapply(mlmAuxol1, function(y){
  do.call(c, sapply(1:6, function(i){
    sapply(Znames[[i]], function(x){x %in% y})}))})

###############################################################################

# Read in Nichols auxotrophs (Supplemental Table 4)
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3060659/bin/NIHMS261392-supplement-04.xls
nicholsAuxo = read.csv("../data/NIHMS261392-supplement-04.csv", 
                       header=TRUE, skip=1)
# Pull out names of Nichols auxotrophs
nicholsAuxoNames = na.omit(
  sapply(strsplit(as.character(nicholsAuxo$Auxotrophs), "-"), 
         function(x){tolower(x[2])})) 
# Get the indices of the mutants corresponding to the Nichols auxotrophs
nicholsAuxoIdx = sapply(Znames, function(x){na.omit(
  match(nicholsAuxoNames, tolower(x)))})
# Pull out the L1-penalized coefficients of the interactions for minimal 
# media conditions and the Nichols auxotrophs
nicholsAuxoMinl1 = sapply(1:6, function(i){
  l1Coeffs[[i]][minimalIdx[[i]],nicholsAuxoIdx[[i]]]})
# Proportion of Nichols auxotrophs with L1-penalized coefficients that 
# have a negative quantile at quantCutoff
mean(do.call(c, sapply(nicholsAuxoMinl1, function(x){
  apply(x, 2, quantile, quantCutoff)})) <= 0)


# Reshape list of Nichols auxotrophs into a data frame
nicholsAuxoMinl1Df = melt(nicholsAuxoMinl1)
names(nicholsAuxoMinl1Df) = c("Cond_Conc", "Mutant", "l1Coeff", "Plate")


png("../pictures/nichols_auxo_dot.png", 
     width=15, height=10, units="cm", res=300)
par(mar=c(2.1,4.1,1.1,1.1))

# Dot plot of coefficients corresponding to Nichols auxotrophs
plot(l1Coeff~as.numeric(Mutant), nicholsAuxoMinl1Df, 
     xlab="", ylab=expression("L"[1]*"-Penalized Interactions"), 
     xaxt="n", pch=16, cex=0.6)
title(xlab="Auxotroph Mutants", line=1)
# Draw horizontal bars at medians
points(unique(as.numeric(nicholsAuxoMinl1Df$Mutant)), 
       tapply(nicholsAuxoMinl1Df$l1Coeff, nicholsAuxoMinl1Df$Mutant, median), 
       pch="-", cex=2.5)
# Horizontal reference line at 0
abline(h=0, col="grey")
dev.off()


# Get the labels for whether or not each mutant is a Nichols auxotroph
nicholsLabels = do.call(c, sapply(1:6, function(i){
  sapply(tolower(Znames[[i]]), function(x){x %in% nicholsAuxoNames})}))

# TPR and FPR when taking the Nichols auxotrophs as the "truth"
# Least squares t-statistics
nicholsTPRl2 = apply(mlmLabelsl2, 2, function(x){
  sum(x==nicholsLabels & x==TRUE)/sum(nicholsLabels)})
nicholsFPRl2 = apply(mlmLabelsl2, 2, function(x){
  1 - sum(x==nicholsLabels & x==FALSE)/sum(nicholsLabels==FALSE)})
# L1-penalized estimates
nicholsTPRl1 = apply(mlmLabelsl1, 2, function(x){
  sum(x==nicholsLabels & x==TRUE)/sum(nicholsLabels)})
nicholsFPRl1 = apply(mlmLabelsl1, 2, function(x){
  1 - sum(x==nicholsLabels & x==FALSE)/sum(nicholsLabels==FALSE)})

png("../pictures/nichols_auxo_ROC.png", 
     width=10, height=10, units="cm", res=300)
par(mar=c(4.1,4.1,1.1,1.1))

# ROC curve for L1-penalized estimates
plot(nicholsFPRl1, nicholsTPRl1,
     xlab="False Positive Rate", ylab="True Positive Rate", 
     xaxs="i", yaxs="i", type="l", 
     cex.lab=0.8, cex.axis=0.8, cex.main=0.8)
# ROC curve for least squares t-statistics
lines(nicholsFPRl2, nicholsTPRl2, lty=2, col="dimgrey")
# Legend
legend("bottomright", bty="n", 
       legend=c(expression("L"[1]*"-Penalized"), "Least Squares"), 
       lty=c(1,2), col=c("black", "dimgrey"), cex=0.8)
# Reference line
abline(0, 1, col="grey", lty=3)
dev.off()

# AUC
# Least squares t-statistics
auc(nicholsFPRl2, nicholsTPRl2, type="spline")
# L1-penalized estimates
auc(nicholsFPRl1, nicholsTPRl1, type="spline")
