# Read in L1 estimates for Agren data
agren_out = as.matrix(read.csv("./processed/agren_lambdas_optim.csv", header=FALSE))

# Lambda = 7.4 (7th lambda), all interactions
# location interactions colored red
png("./pictures/Agren_site_allQTL_lambda7.4.png", width=380, height=150)
par(mar=c(1.1,4.1,1.1,1.1))

# Site interactions
plot(1:700, abs(agren_out[(700+1):(2*700),7]), type="l", 
     xlab="", ylab="Interactions", xaxt="n", ylim=c(0,0.055), 
     col="firebrick3")

# All other interactions
lines(1:700, abs(agren_out[(2*700+1):(3*700),7]), lty=2)
lines(1:700, abs(agren_out[(3*700+1):(4*700),7]), lty=2)
lines(1:700, abs(agren_out[(4*700+1):(5*700),7]), lty=2)
lines(1:700, abs(agren_out[(5*700+1):(6*700),7]), lty=2)
lines(1:700, abs(agren_out[(6*700+1):(7*700),7]), lty=2)
lines(1:700, abs(agren_out[(7*700+1):(8*700),7]), lty=2)

# Reference lines
invisible(sapply(cumsum(c(173, 118, 137, 113)), 
                 function(v){abline(v=v, col="grey", lty=3)}))
dev.off()
