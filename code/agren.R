# Read in L1 coefficient estimates
agrenL1Coeffs = as.matrix(read.csv("../processed/agren_l1_coeffs.csv", header=FALSE))
# Pull out coefficients corresponding to lambda of interest
lambda = 8 # MSE is minimized here (lambda = 6.19174)
coeffs = matrix(agrenL1Coeffs[,lambda], 700, 2) 


png("../pictures/agren_l1_coeffs_lambda_6.2.png", width=380, height=150)
par(mar=c(1.1,4.1,1.1,1.1))

# Main effects
plot(1:700, abs(coeffs[,1]), type="l", 
     xlab="", ylab="QTL Effects", xaxt="n", 
     ylim=c(min(-abs(coeffs[,2])), max(abs(coeffs[,1]))))

# Site interactions
lines(1:700, -abs(coeffs[,2]))

# Reference lines to delineate chromosomes
invisible(sapply(cumsum(c(173, 118, 137, 113)), 
                 function(v){abline(v=v, col="grey", lty=3)}))
dev.off()
