# Read in L1 coefficient estimates
agrenL1Coeffs = as.matrix(read.csv("../processed/agren_l1_coeffs.csv", header=FALSE))
# Pull out coefficients corresponding to lambda of interest
lambda = 8 # MSE is minimized here (lambda = 6.19174)
coeffs = matrix(agrenL1Coeffs[,lambda], 700, 2) 


png("../pictures/agren_l1_coeffs_lambda_6.2.png", width=380, height=150)
par(mar=c(1.1,4.1,1.1,1.1))

# Main effects
plot(1:700, abs(coeffs[,1]), type="l", 
     xlab="", xaxt="n", yaxt="n", 
     ylab="Absolute QTL Effects", 
     ylim=max(abs(coeffs))*c(-1,1))

# Site interactions
lines(1:700, -abs(coeffs[,2]))

# Reference lines to delineate chromosomes
invisible(sapply(cumsum(c(173, 118, 137, 113)), 
                 function(v){abline(v=v, col="grey", lty=3)}))

# Customized y-axis without negative labels
ticks = seq(-0.15, 0.15, by=0.05)
axis(2, at = ticks, labels = abs(ticks))

# Text labels for main and interaction effects
text(700, 0.13, "Main effects", adj=1)
text(700, -0.13, "Interactions", adj=1)
dev.off()
