# Read in L1 coefficient estimates
agrenL1Coeffs = as.matrix(read.csv("./processed/agren_l1_coeffs.csv", header=FALSE))
# Pull out coefficients corresponding to lambda of interest
lambda = 7 # 7th lambda has value 7.4
coeffs = matrix(agrenL1Coeffs[,7], 700, 8) 


png("./pictures/agren_l1_coeffs_lambda7.4.png", width=380, height=150)
par(mar=c(1.1,4.1,1.1,1.1))

# Site interactions
plot(1:700, abs(coeffs[,2]), type="l", 
     xlab="", ylab="Interactions", xaxt="n", ylim=c(0,0.1), 
     col="firebrick3")

# All other interactions
lines(1:700, abs(coeffs[,3]), lty=2)
lines(1:700, abs(coeffs[,4]), lty=2)
lines(1:700, abs(coeffs[,5]), lty=2)
lines(1:700, abs(coeffs[,6]), lty=2)
lines(1:700, abs(coeffs[,7]), lty=2)
lines(1:700, abs(coeffs[,8]), lty=2)

# Reference lines to delineate chromosomes
invisible(sapply(cumsum(c(173, 118, 137, 113)), 
                 function(v){abline(v=v, col="grey", lty=3)}))
dev.off()
