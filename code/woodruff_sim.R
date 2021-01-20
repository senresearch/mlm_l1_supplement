library(MESS) # AUC
library(mutoss) # adaptive Benjamini-Hochberg

# 100 chemicals
nChem = 100
# 10 tissues
nTiss = 10

# Read in simulated X matrix of demographics
X = read.csv("../processed/woodruff_sim_X.csv", header=FALSE)
# Read in simulated Y matrix 
Y = read.csv("../processed/woodruff_sim_Y.csv", header=FALSE)
# Read in simulated interactions
interactions = read.csv("../processed/woodruff_sim_interactions.csv", 
                        header=FALSE)

# Read in L1 coefficient estimates
woodruffSimL1Coeffs = read.csv("../processed/woodruff_sim_l1_coeffs.csv", 
                               header=FALSE)
# Reshape the coefficients so that there is a separate matrix for each lambda
coeffsList = lapply(1:length(woodruffSimL1Coeffs), function(j){
  matrix(woodruffSimL1Coeffs[,j], 
         dim(interactions)[1]+1, dim(interactions)[2]+1)
})


# Run linear regression one-by-one for each chemical x tissue combination 
# and get p-values
lmPvals = sapply(1:ncol(Y), function(j){
  summary(lm(y~., data=data.frame(y=Y[,j], X)))$coefficients[-1,4]
})
# Apply adaptive Benjamini-Hochberg adjustment to p-values
adjlmPvals = matrix(adaptiveBH(as.matrix(lmPvals), alpha=0.05, 
                               silent=TRUE)$adjPValues, ncol=ncol(lmPvals))


# Positives
posMat = interactions[,-(1:(nChem+nTiss))] != 0
# Negatives
negMat = interactions[,-(1:(nChem+nTiss))] == 0
# FDR cutoffs, chosen to get smooth curves
FDR = c(10^(seq(-164, -44, by=20)), 10^(seq(-40, -4, by=4)), 
        seq(0.1,1,length=20))


# Function to map interaction space to chemicals
# interactions = interactions
# nChem = number of chemicals
# nTiss = number of tissues
# fun = function to apply for reduction
inter2chem = function (interactions, nChem, nTiss, fun=mean) {
  # Initialize output matrix
  chems = matrix(nrow=nrow(interactions), ncol=nChem)
  # Reduce interactions by chemical
  for (j in 1:nChem) {
    chems[,j] = apply(interactions[,seq(j, by=nChem, length=nTiss)], 1, fun) 
  }
  return(chems)
}

# Repeated matrix of simulated chemical effects to map to chemicals
interStackChem = repmat(as.matrix(interactions[,1:nChem]), ncol=nTiss) 


# TPR for L1 coefficient estimates
mlmL1TPR = sapply(1:length(coeffsList), function(i){
  sum((coeffsList[[i]][-1, -1][,-(1:(nChem+nTiss))] != 0) & posMat)
}) / sum(posMat)
# FPR for L1 coefficient estimates
mlmL1FPR = sapply(1:length(coeffsList), function(i){
  sum((coeffsList[[i]][-1, -1][,-(1:(nChem+nTiss))] != 0) & negMat)
}) / sum(negMat)


# TPR for univariate linear regression
lmTPR = sapply(1:length(FDR), function(i){
  sum((adjlmPvals <= FDR[i]) & posMat)
}) / sum(posMat)
# FPR for univariate linear regression
lmFPR = sapply(1:length(FDR), function(i){
  sum((adjlmPvals <= FDR[i]) & negMat)
}) / sum(negMat)


# TPR for univariate linear regression when using number of hits to detect 
# chemical significance
lmChemTPR = sapply(1:4, function(k){
  sapply(1:length(FDR), function(i){
    sum(inter2chem(((adjlmPvals <= FDR[i]) & 
                      (interStackChem != 0)), nChem, nTiss) > k/5)
  })
}) / sum(interactions[,1:nChem] != 0)
# FPR for univariate linear regression when using number of hits to detect 
# chemical significance
lmChemFPR = sapply(1:length(FDR), function(i){
  sum(inter2chem(((adjlmPvals <= FDR[i]) & 
                    (interStackChem == 0)), nChem, nTiss) > 1/5)
}) / sum(interactions[,1:nChem] == 0)


# Calculate AUC for each method
AUCs = c(MESS::auc(c(1,mlmL1FPR,0), c(1,mlmL1TPR,0), type="spline"), 
         MESS::auc(c(0,lmFPR,1), c(0,lmTPR,1), type="spline"),
         sapply(1:4, function(j){
           auc(c(0,lmChemFPR,1), c(0,lmChemTPR[,j],1), type="spline")}))


png("../pictures/woodruff_sim_ROC_chem.png", width=380, height=380)
par(mar=c(4.1,4.1,1.1,1.1))

# L1-penalized matrix linear models
plot(c(mlmL1FPR,0), c(mlmL1TPR,0), 
     xlab="False Positive Rate", ylab="True Positive Rate", 
     xaxs="i", yaxs="i", type="l", lwd=1.5) 

# Univariate linear regression for interactions
lines(lmFPR, lmTPR, lty="dotted", lwd=1.5) 

# Univariate linear regression with at least 1/5 chemical hits
lines(lmChemFPR, lmChemTPR[,1], col="grey", lty="dashed")
# Univariate linear regression with at least 2/5 chemical hits
lines(lmChemFPR, lmChemTPR[,2], col="grey", lty="dotdash")
# Univariate linear regression with at least 3/5 chemical hits
lines(lmChemFPR, lmChemTPR[,3], col="grey", lty="solid")
# Univariate linear regression with at least 4/5 chemical hits
lines(lmChemFPR, lmChemTPR[,4], col="grey", lty="dotted")

# Reference line
abline(0, 1, col="grey") 

# Legend for different methods
legChem = paste0(c("Interactions",  "1/5 Hits", "2/5 Hits", 
                    "3/5 Hits", "4/5 Hits"), 
                  " (", round(AUCs[-1], 3), ")")
legend(0.4, 0.4, 
       c(as.expression(bquote(L[1]-penalized ~ (.(round(AUCs[1], 3))))), 
         sapply(legChem, function(x){bquote(.(x))})), 
       col=c("black", "black", rep("grey", 4)), 
       lty=c("solid", "dotted", "dashed", 
             "dotdash", "solid", "dotted"), 
       lwd=c(rep(1.5,2), rep(1,4)), bty="n")
dev.off()
