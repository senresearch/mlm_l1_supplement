# Read in L1 estimates for Agren data
coeffs = as.matrix(read.csv("./processed/agren_l1_coeffs.csv", header=FALSE))

coeffs_lambda7.4 = matrix(coeffs[,7], 700, 8)

# Lambda = 7.4 (7th lambda), all interactions
# location interactions colored red
png("./pictures/agren_l1_coeffs_lambda7.4.png", width=380, height=150)
par(mar=c(1.1,4.1,1.1,1.1))

# Site interactions
plot(1:700, abs(coeffs_lambda7.4[,2]), type="l", 
     xlab="", ylab="Interactions", xaxt="n", ylim=c(0,0.2), 
     col="firebrick3")

# All other interactions
lines(1:700, abs(coeffs_lambda7.4[,3]), lty=2)
lines(1:700, abs(coeffs_lambda7.4[,4]), lty=2)
lines(1:700, abs(coeffs_lambda7.4[,5]), lty=2)
lines(1:700, abs(coeffs_lambda7.4[,6]), lty=2)
lines(1:700, abs(coeffs_lambda7.4[,7]), lty=2)
lines(1:700, abs(coeffs_lambda7.4[,8]), lty=2)

# Reference lines
invisible(sapply(cumsum(c(173, 118, 137, 113)), 
                 function(v){abline(v=v, col="grey", lty=3)}))
dev.off()
