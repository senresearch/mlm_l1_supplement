using Distributed
using DataFrames
using LinearAlgebra
import LinearAlgebra.I
using Random
using CSV

# L1-penalized matrix linear models
using matrixLMnet


# Read in X (genotype probabilities) with cytoplasm contrast. 
# The first row is a header. 
X = convert(Array{Float64, 2}, CSV.read("../processed/lowry_cyto_genoprobs.csv", 
                                        delim=',', header=true))

# Read in Y (phenotypes). The first row is a header. 
Y = convert(Array{Float64, 2}, CSV.read("../processed/lowry_pheno.csv", 
                                        delim=',', header=true))

# Number of phenotypes 
nPheno = convert(Int64, size(Y,2)/2) 
# Create Z matrix. 
Z = kron(Matrix{Float64}(I, nPheno, nPheno), vcat([1 1], [1 -1]))

# Put together RawData object for MLM
MLMData = RawData(Response(Y), Predictors(X, Z))

# Array of 16 lambdas
lambdas = reverse(1.2.^(-9:6))


# Randomly sample a small number of phenotypes for a dry run
randSamp = 500 # Number of random phenotypes to sample
Random.seed!(10)
randPhenos = sort(Random.shuffle(1:nPheno)[1:randSamp])
idx = Array{Int64}(undef, 2*randSamp)
for i in 1:randSamp
	  idx[2*i-1] = randPhenos[i]*2-1
	  idx[2*i] = randPhenos[i]*2
end

# Subset the Y matrix to keep only these phenotypes
YSub = Y[:,idx]

# Create the Z matrix for this subset of phenotypes
nPhenoSub = convert(Int64, size(YSub,2)/2) 
ZSub = kron(Matrix{Float64}(I, nPhenoSub, nPhenoSub), vcat([1 1], [1 -1]))

# Put together RawData object for MLM using subsetted data
MLMDataSub = RawData(Response(YSub), Predictors(X, ZSub))

# Dry run of L1-penalized matrix linear model using subset of phenotypes
resultsSub = mlmnet(fista_bt!, MLMDataSub, lambdas; isZIntercept=false) 


# Run L1-penalized matrix linear model while timing it
start = time()
results = mlmnet(fista_bt!, MLMData, lambdas; isZIntercept=false) 
elapsedTime = time() - start

# Flatten coefficients and write results to CSV
flat_coeffs = coef_2d(results)
CSV.write("../processed/lowry_l1_coeffs.csv", 
          DataFrame(flat_coeffs), writeheader=false)

# Print and write time to CSV
println(elapsedTime)
CSV.write("../processed/lowry_time.csv", 
          DataFrame(time=elapsedTime), writeheader=false)
