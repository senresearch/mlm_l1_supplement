library(MESS) # AUC
library(mutoss) # adaptive Benjamini-Hochberg

# Read in interactions
interactions = lapply(1:6, function(i){
  out = read.csv(paste("../processed/nichols_sim_p", i, 
                       "_interactions.csv", sep=""), header=FALSE)
  return(out[-nrow(out), -ncol(out)])
})


# Read in MLM p-values 
pvals = lapply(1:6, function(i){
  out = read.csv(paste("../processed/nichols_sim_p", i, 
                       "_pvals.csv", sep=""), header=FALSE)
  return(out[-nrow(out), -ncol(out)])
})

# Convert MLM p-values to adaptive BH-adjusted p-values
adjPvals = lapply(1:6, function(i){
  matrix(adaptiveBH(as.matrix(pvals[[i]]), 
                    alpha=0.05, silent=TRUE)$adjPValues, 
         ncol=ncol(pvals[[i]]))
})


# Read in L1 coefficient estimates
l1Coeffs = lapply(1:6, function(i){
  read.csv(paste("../processed/nichols_sim_p", i, "_l1_coeffs.csv", sep=""), 
           sep=",", header=FALSE)
})

# Reshape the coefficients so that there is a separate matrix for each lambda
coeffsList = lapply(1:length(l1Coeffs), function(i){
  lapply(1:ncol(l1Coeffs[[i]]), function(j){
    matrix(l1Coeffs[[i]][,j], 
           dim(interactions[[i]])[1]+1, dim(interactions[[i]])[2]+1)[-1,-1]
  })
})


# Positives
posMat = lapply(interactions, function(x){ x != 0})
# Negatives
negMat = lapply(interactions, function(x){ x == 0}) 


# Range of FDR cutoffs for least squares p-values
FDRs = seq(0, 1, length=50)

# TPR least squares p-values
mlmTPR = lapply(1:6, function(i){
  sapply(FDRs, function(FDR){
    sum((adjPvals[[i]] <= FDR) & posMat[[i]]) / 
      sum(posMat[[i]])
  })
})
# FPR least squares p-values
mlmFPR = lapply(1:6, function(i){
  sapply(FDRs, function(FDR){
    sum((adjPvals[[i]] <= FDR) & negMat[[i]]) / 
      sum(negMat[[i]])
  })
})


# TPR for L1 coefficient estimates
mlmL1TPR = lapply(1:length(coeffsList), function(i){
  sapply(1:length(coeffsList[[i]]), function(j){
    sum((coeffsList[[i]][[j]] != 0) & posMat[[i]])
  }) / sum(posMat[[i]])
})
# FPR for L1 coefficient estimates
mlmL1FPR = lapply(1:length(coeffsList), function(i){
  sapply(1:length(coeffsList[[i]]), function(j){
    sum((coeffsList[[i]][[j]] != 0) & negMat[[i]])
  }) / sum(negMat[[i]])
})


png("../pictures/nichols_sim_p%01d_ROC.png", 
    width=10, height=10, units="cm", res=300)
par(mar=c(4.1,4.1,1.1,1.1))

AUCs = sapply(1:6, function(i) {
  # ROC curve for MLM
  plot(c(0, sort(mlmL1FPR[[i]])), c(0, sort(mlmL1TPR[[i]])), 
       xlab="False Positive Rate", ylab="True Positive Rate", 
       main=paste("Plate", i), xaxs="i", yaxs="i", type="l", 
       cex.lab=0.8, cex.axis=0.8, cex.main=0.8)
  # ROC curve for S scores
  lines(c(0, mlmFPR[[i]]), c(0, mlmTPR[[i]]),
        type="l", lty=2, col="dimgrey")
  
  # Legend for different methods
  legend("bottomright", bty="n", 
         legend=c(expression("L"[1]*"-Penalized"), "Least Squares"), 
         lty=c(1,2), col=c("black", "dimgrey"), cex=0.8)
  # Reference line
  abline(0, 1, col="grey", lty=3)
  
  # Return AUCs
  return(c(auc(c(0, mlmL1FPR[[i]]), c(0, mlmL1TPR[[i]]), type="spline"), 
           auc(c(0, mlmFPR[[i]]), c(0, mlmTPR[[i]]), type="spline")))
})
dev.off()

rownames(AUCs) = c("L1-Penalized", "Least Squares")
colnames(AUCs) = paste("Plate", 1:6)
