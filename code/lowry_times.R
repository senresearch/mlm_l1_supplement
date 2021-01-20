
randSamp = c(25, 50, 100, 200, 400, 800, 1600)
for (r in randSamp) {
  old_times = read.csv(paste0("../processed/lowry_times_", r, "_old.csv"), 
                       row.names=1)
  new_times = read.csv(paste0("../processed/lowry_times_", r, "_new.csv"), 
                       row.names=1)
  
  agg_times = cbind(old_times[,-1], new_times[,-1])
  agg_times = data.frame(method=rownames(agg_times), 
                         mean=rowMeans(agg_times), agg_times)
  colnames(agg_times)[-(1:2)] = 1:(ncol(agg_times)-2)
  write.csv(agg_times, paste0("../processed/lowry_times_", r, ".csv"), 
            row.names=FALSE)
}


# Number of random genes that were sampled
randSamp = c(25, 50, 100, 200, 400, 800, 1600)

# Runtimes for as the the number of genes increases
lowryTimes = lapply(randSamp, function(r) {
  read.csv(paste0("../processed/lowry_times_", r, ".csv"), row.names=1)
})

# Extract mean runtimes across replicates
meanTimes = t(do.call(cbind, lapply(lowryTimes, function(x) { x["mean"]})))
rownames(meanTimes) = randSamp


# Point shapes for plotting
myPch = c(1, 0, 16, 15, 17, 4)


png("../pictures/lowry_times.png",  width=450, height=350)
par(mar=c(4.1,4.1,1.1,1.1))

# Plot runtimes against the number of genes
matplot(log(randSamp, 2), log(meanTimes, 2), type="b", xaxt="n", yaxt="n", 
        col="black", pch=myPch, lty=1:6, 
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
       bty="n", pch=myPch, lty=1:6)
dev.off()
