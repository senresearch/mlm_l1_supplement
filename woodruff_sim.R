library(MESS) # AUC
library(mutoss) # adaptive Benjamini-Hochberg

# 100 chemicals
nchem = 100
# 10 tissues
ntiss = 10

# Read in simulated X matrix of demographics
X = read.csv("./processed/woodruff_sim_X.csv", header=FALSE)
# Read in simulated Y matrix 
Y = read.csv("./processed/woodruff_sim_Y.csv", header=FALSE)
# Read in simulated interactions
interactions = read.csv("./processed/woodruff_sim_interactions.csv", 
                        header=FALSE)

# Read in L1 coefficient estimates
woodruff_sim_l1_coeffs = read.csv("./processed/woodruff_sim_l1_coeffs.csv", 
                                  header=FALSE)
# Reshape the coefficients so that there is a separate matrix for each lambda
coeffs_list = lapply(1:length(woodruff_sim_l1_coeffs), function(j){
  matrix(woodruff_sim_l1_coeffs[,j], 
         dim(interactions)[1]+1, dim(interactions)[2]+1)
})


# Run linear regression one-by-one for each chemical x tissue combination 
# and get p-values
p_lm = sapply(1:ncol(Y), function(j){
  summary(lm(y~., data=data.frame(y=Y[,j], X)))$coefficients[-1,4]
})
# Apply adaptive Benjamini-Hochberg adjustment to p-values
p.adj_lm = matrix(adaptiveBH(as.matrix(p_lm), alpha=0.05, 
                             silent=TRUE)$adjPValues, ncol=ncol(p_lm))


# Positives
P_mat = interactions[,-(1:(nchem+ntiss))] != 0
# Negatives
N_mat = interactions[,-(1:(nchem+ntiss))] == 0
# FDR cutoffs, chosen to get smooth curves
FDR = c(10^(seq(-164, -44, by=20)), 10^(seq(-40, -4, by=4)), 
        seq(0.1,1,length=20))


# Function to map interaction space to chemicals
# interact = interactions
# nchem = number of chemicals
# ntiss = number of tissues
# FUN = function to apply for reduction
inter2chem = function (interact, nchem, ntiss, FUN=mean) {
  # Initialize output matrix
  chem = matrix(nrow=nrow(interact), ncol=nchem)
  # Reduce interactions by chemical
  for (j in 1:nchem) {
    chem[,j] = apply(interact[,seq(j, by=nchem, length=ntiss)], 1, FUN) 
  }
  return(chem)
}

# Repeated matrix of simulated chemical effects to map to chemicals
inter_chem_stack = repmat(as.matrix(interactions[,1:nchem]), ncol=ntiss) 


# TPR for L1 coefficient estimates
tpr_mlm_l1 = sapply(1:length(coeffs_list), function(i){
  sum((coeffs_list[[i]][-1, -1][,-(1:(nchem+ntiss))] != 0) & P_mat)
}) / sum(P_mat)
# FPR for L1 coefficient estimates
fpr_mlm_l1 = sapply(1:length(coeffs_list), function(i){
  sum((coeffs_list[[i]][-1, -1][,-(1:(nchem+ntiss))] != 0) & N_mat)
}) / sum(N_mat)


# TPR for univariate linear regression
tpr_lm = sapply(1:length(FDR), function(i){
  sum((p.adj_lm <= FDR[i]) & P_mat)
}) / sum(P_mat)
# FPR for univariate linear regression
fpr_lm = sapply(1:length(FDR), function(i){
  sum((p.adj_lm <= FDR[i]) & N_mat)
}) / sum(N_mat)


# TPR for univariate linear regression when using number of hits to detect 
# chemical significance
tpr_lm_chem = sapply(1:4, function(k){
  sapply(1:length(FDR), function(i){
    sum(inter2chem(((p.adj_lm <= FDR[i]) & 
                      (inter_chem_stack != 0)), nchem, ntiss) > k/5)
  })
}) / sum(interactions[,1:nchem] != 0)
# FPR for univariate linear regression when using number of hits to detect 
# chemical significance
fpr_lm_chem = sapply(1:length(FDR), function(i){
  sum(inter2chem(((p.adj_lm <= FDR[i]) & 
                    (inter_chem_stack == 0)), nchem, ntiss) > 1/5)
}) / sum(interactions[,1:nchem] == 0)


# Calculate AUC for each method
AUCs = c(MESS::auc(c(1,fpr_mlm_l1,0), c(1,tpr_mlm_l1,0), type="spline"), 
         MESS::auc(c(0,fpr_lm,1), c(0,tpr_lm,1)),
         sapply(1:4, function(j){
           auc(c(0,fpr_lm_chem,1), c(0,tpr_lm_chem[,j],1))}))


png("./pictures/woodruff_sim_ROC_chem.png", width=380, height=380)
par(mar=c(4.1,4.1,1.1,1.1))

# L1-penalized matrix linear models
plot(c(fpr_mlm_l1,0), c(tpr_mlm_l1,0), 
     xlab="False Positive Rate", ylab="True Positive Rate", 
     xaxs="i", yaxs="i", type="l") 

# Univariate linear regression for interactions
lines(fpr_lm, tpr_lm, col="firebrick3") 

# Univariate linear regression with at least 1/5 chemical hits
lines(fpr_lm_chem, tpr_lm_chem[,1], col="dodgerblue3")
# Univariate linear regression with at least 2/5 chemical hits
lines(fpr_lm_chem, tpr_lm_chem[,2], col="dodgerblue3", lty="dashed")
# Univariate linear regression with at least 3/5 chemical hits
lines(fpr_lm_chem, tpr_lm_chem[,3], col="dodgerblue3", lty="dotdash")
# Univariate linear regression with at least 4/5 chemical hits
lines(fpr_lm_chem, tpr_lm_chem[,4], col="dodgerblue3", lty="dotted")

# Reference line
abline(0,1, col="grey") 

# Legend for different methods
leg_chem = paste0(c("Interactions",  "1/5 Hits", "2/5 Hits", 
                    "3/5 Hits", "4/5 Hits"), 
                  " (", round(AUCs[-1], 3), ")")
legend(0.4, 0.4, 
       c(bquote(L[1]-penalized~(.(round(AUCs[1], 3)))), 
         sapply(leg_chem, function(x){bquote(.(x))})), 
       col=c("black", "firebrick3", rep("dodgerblue3", 4)), 
       lty=c(rep("solid",3), "dashed", "dotdash", "dotted"), bty="n")
dev.off()
