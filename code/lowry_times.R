# Number of random genes that were sampled
randSamp = c(25, 50, 100, 200, 400, 800, 1600)

# Runtimes for as the the number of genes increases
lowryTimes = lapply(randSamp, function(r) {
  read.csv(paste0("../processed/lowry_times_", r, ".csv"), row.names=1)
})

# Extract mean runtimes across replicates
meanTimes = t(do.call(cbind, lapply(lowryTimes, function(x) { x["mean"]})))
rownames(meanTimes) = randSamp


# Colors and point shapes for plotting
myCols = c("darkred", "firebrick2", 
           "darkblue", "dodgerblue3", "deepskyblue1", 
           "darkolivegreen4")
myPch = c(1, 0, 16, 15, 17, 4)


png("../pictures/lowry_times.png",  width=450, height=350)
par(mar=c(4.1,4.1,1.1,1.1))

# Plot runtimes against the number of genes
matplot(log(randSamp, 2), log(meanTimes, 2), type="b", xaxt="n", yaxt="n", 
        pch=myPch, lty=1:6, col=myCols, 
        xlab="Number of Genes", ylab="Runtime (seconds)")

# Manually add x-axis on the log2 scale
axis(1, log(randSamp, 2), randSamp)
# Manually add y-axis on the log2 scale
yTicks = c(5, 20, 80, 320, 1280, 5120, 20480)
axis(2, log(yTicks, 2), yTicks)

# Legend for different estimation methods
legend("topleft", 
       legend=c("Coordinate descent. (cyclic)", "Coordinate descent (random)", 
                "ISTA (fixed step size)", "FISTA (fixed step size)", 
                "FISTA (backtracking)", "ADMM"), 
       bty="n", pch=myPch, lty=1:6, col=myCols)
dev.off()
