# Runtimes for FISTA with backtracking as the dimensions scale up
FISTApqTimes = read.csv("../processed/FISTA_pq_times.csv")
FISTApqMeans = matrix(FISTApqTimes$mean, ncol=5)
FISTAnmTimes = read.csv("../processed/FISTA_nm_times.csv")
FISTAnmMeans = matrix(FISTAnmTimes$mean, ncol=5)

# Runtimes for ADMM as the dimensions scale up
ADMMpqTimes = read.csv("../processed/ADMM_pq_times.csv")
ADMMpqMeans = matrix(ADMMpqTimes$mean, ncol=5)
ADMMnmTimes = read.csv("../processed/ADMM_nm_times.csv")
ADMMnmMeans = matrix(ADMMnmTimes$mean, ncol=5)

# Values of p & q and n & m
p = q = c(200, 400, 600, 800, 1000)
n = m = c(400, 800, 1200, 1600, 2000)

# Function used with `image` to get the plot to match the matrix.
orient = function(mat){
  return(t(mat)[,nrow(mat):1])
}


# Breaks and colors for heatmaps of the difference between runtimes
diffBreaks = c(-205, -120, seq(-100, 100, length=8), 120, 205)
diffCols = c(colorRampPalette(colors = c("blue","white"))(6), 
             colorRampPalette(colors = c("white","red"))(6)[-1])

# Heatmap of the difference between FISTA and ADMM runtimes, varying p and q
image(q, p, orient(FISTApqMeans - ADMMpqMeans), 
      col=diffCols, breaks=diffBreaks, yaxt="n", 
      main="FISTA - ADMM runtimes \n n=m=1200")
axis(2, p, rev(p))
text(matrix(rep(p, 5), ncol=5), matrix(rep(q, each=5), ncol=5), 
     round(orient(FISTApqMeans - ADMMpqMeans), 2))

# Heatmap of the difference between FISTA and ADMM runtimes, varying n and m
image(m, n, orient(FISTAnmMeans - ADMMnmMeans), 
      col=diffCols, breaks=diffBreaks, yaxt="n", 
      main="FISTA - ADMM runtimes \n p=q=400")
axis(2, n, rev(n))
text(matrix(rep(n, 5), ncol=5), matrix(rep(m, each=5), ncol=5), 
     round(orient(FISTAnmMeans - ADMMnmMeans), 2))


# Breaks and colors for heatmaps of the ratio of runtimes
ratioBreaks = seq(0.3, 1.7, length=12)
ratioCols = c(colorRampPalette(colors = c("black","white"))(6), 
              colorRampPalette(colors = c("white","black"))(6)[-1])

# Heatmap of the ratio of FISTA and ADMM runtimes, varying p and q
image(q, p, orient(FISTApqMeans/ADMMpqMeans), 
      col=ratioCols, breaks=ratioBreaks, yaxt="n", 
      main="FISTA / ADMM runtimes \n n=m=1200")
axis(2, p, rev(p))
text(matrix(rep(p, 5), ncol=5), matrix(rep(q, each=5), ncol=5), 
     round(orient(FISTApqMeans/ADMMpqMeans), 2))

# Heatmap of the ratio of FISTA and ADMM runtimes, varying n and m
image(m, n, orient(FISTAnmMeans/ADMMnmMeans), 
      col=ratioCols, breaks=ratioBreaks, yaxt="n", 
      main="FISTA / ADMM runtimes \n p=q=400")
axis(2, n, rev(n))
text(matrix(rep(n, 5), ncol=5), matrix(rep(m, each=5), ncol=5), 
     round(orient(FISTAnmMeans/ADMMnmMeans), 2))
