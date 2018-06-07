library(MESS) # For AUC
library(mutoss) # For adaptive BH

# 100 chemicals
nchem = 100
# 10 tissues
ntiss = 10

# Simulated X matrix of demographics
X = read.csv("./processed/enviro_Xnoint_simX.csv", header=FALSE)
# Simulated Y matrix 
Y = read.csv("./processed/enviro_Y_sim_simX.csv", header=FALSE)
# Simulated interactions
inter = read.csv("./processed/enviro_sim_inter_simX.csv", header=FALSE)

# Lambdas used for L1
lambdas = as.character(read.csv("./processed/enviro_lambdas.csv", header=FALSE)[,1])
# Some annoying manipulation to get lambdas to match file names
lambdas[lambdas=="1"] = "1.0"
lambdas[lambdas=="5e-04"] = "0.0005"
lambdas[lambdas=="1e-04"] = "0.0001"
lambdas[lambdas=="8e-05"] = "8.0e-5"
lambdas[lambdas=="6e-05"] = "6.0e-5"

# Read in coefficient estimates from L1
coeffs_l1 = lapply(lambdas, function(lambda){
  read.csv(paste("./processed/enviro_sim_lambda", lambda, "_foldall_simX.csv", sep=""), header=FALSE)
})

# TPR and FPR for L1
tpr_l1 = sapply(1:length(lambdas), function(i){
  sum((coeffs_l1[[i]][-1, -1][,-(1:(nchem+ntiss))] == 0) & (inter[,-(1:(nchem+ntiss))] == 0)) /sum(inter[,-(1:(nchem+ntiss))] == 0)
})
fpr_l1 = 1 - sapply(1:length(lambdas), function(i){
  sum((coeffs_l1[[i]][-1, -1][,-(1:(nchem+ntiss))] != 0) & (inter[,-(1:(nchem+ntiss))] != 0)) /sum(inter[,-(1:(nchem+ntiss))] != 0)
})


# Run linear regression one by one for each chemical x tissue and get p-values
pvals_linreg = sapply(1:ncol(Y), function(j){
  summary(lm(y~., data=data.frame(y=Y[,j], X)))$coefficients[-1,4]
})
# B-H adjusted p-values
pvalsa_linreg = matrix(adaptiveBH(as.matrix(pvals_linreg), alpha=0.05, silent=TRUE)$adjPValues, ncol=ncol(pvals_linreg))

# Cutoffs
FDR = c(10^(seq(-164, -44, by=20)), 10^(seq(-40, -4, by=4)), seq(0.1,1,length=20))

# TPR and FPR from univariate linear regression
tprall = sapply(1:length(FDR), function(i){
  sum((pvalsa_linreg <= FDR[i]) & (inter[,-(1:(nchem+ntiss))] != 0)) /sum(inter[,-(1:(nchem+ntiss))] != 0)
})
fprall = sapply(1:length(FDR), function(i){
  sum((pvalsa_linreg <= FDR[i]) & (inter[,-(1:(nchem+ntiss))] == 0)) /sum(inter[,-(1:(nchem+ntiss))] == 0)
})


# Map interaction space to chemicals
B_reduce_chem = function (B, nchem, ntiss, fun=mean) {
  B_reduced = matrix(nrow=nrow(B), ncol=nchem)
  for (j in 1:nchem) {
    B_reduced[,j] = apply(B[,seq(j, by=nchem, length=ntiss)], 1, fun) 
  }
  return(B_reduced)
}

# Repeated matrix of simulated chemical effects to map to chemicals
inter_stack_chem = repmat(as.matrix(inter[,1:nchem]), ncol=ntiss) 

# TPR and FPR when using number of hits to detect chemical significance
tpr_chem = sapply(1:4, function(k){
  sapply(1:length(FDR), function(i){
    sum(B_reduce_chem(((pvalsa_linreg <= FDR[i]) & (inter_stack_chem != 0)), nchem, ntiss) > k/5) /sum(inter[,1:nchem] != 0)
  })
})
fpr_chem = sapply(1:4, function(k){
  sapply(1:length(FDR), function(i){
    sum(B_reduce_chem(((pvalsa_linreg <= FDR[i]) & (inter_stack_chem == 0)), nchem, ntiss) > 1/5) /sum(inter[,1:nchem] == 0)
  })
})
sapply(1:4, function(j){auc(c(0,fpr_chem[,j],1), c(0,tpr_chem[,j],1))})



auc_chem = c(MESS::auc(c(1,fpr_l1,0), c(1,tpr_l1,0), type="spline"),
             MESS::auc(c(0,fprall,1), c(0,tprall,1)),
             sapply(1:4, function(j){auc(c(0,fpr_chem[,j],1), c(0,tpr_chem[,j],1))}))
png("./pictures/enviro_sim_ROC_simX_chem.png", width=380, height=380)
par(mar=c(4.1,4.1,1.1,1.1))
plot(c(fpr_l1,0), c(tpr_l1,0), xlab="False Positive Rate", ylab="True Positive Rate", 
     xaxs="i", yaxs="i", type="l") # L1
lines(fprall, tprall, col="firebrick3") # Univariate linear regression for interactions

# At least 1, 2, 3, and 4 hits
lines(fpr_chem[,1], tpr_chem[,1], col="dodgerblue3")
lines(fpr_chem[,2], tpr_chem[,2], col="dodgerblue3", lty="dashed")
lines(fpr_chem[,3], tpr_chem[,3], col="dodgerblue3", lty="dotdash")
lines(fpr_chem[,4], tpr_chem[,4], col="dodgerblue3", lty="dotted")

abline(0,1, col="grey") # Reference

leg_chem = paste0(c("Interactions",  "1/5 Hits", "2/5 Hits", "3/5 Hits", "4/5 Hits"), 
                  " (", round(auc_chem[-1], 3), ")")
legend(0.4, 0.4, 
       c(expression(L[1]-penalized ~ (0.905)), leg_chem), 
       col=c("black", "firebrick3", rep("dodgerblue3", 4)), 
       lty=c(rep("solid",3), "dashed", "dotdash", "dotted"), bty="n")
dev.off()
