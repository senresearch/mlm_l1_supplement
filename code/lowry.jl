using Distributed
using DataFrames
using LinearAlgebra
using Random
using CSV

# L1-penalized matrix linear models
@everywhere include("../../mlm_packages/matrixLMnet/src/matrixLMnet.jl")
@everywhere using Main.matrixLMnet


# Read in X (genotype probabilities) with cytoplasm contrast. 
# The first row is a header. 
X = convert(Array{Float64, 2}, CSV.read("../processed/lowry_cyto_genoprobs.csv", 
                                        delim=',', header=true))

# Read in Y (phenotypes). The first row is a header. 
Y = convert(Array{Float64, 2}, CSV.read("../processed/lowry_pheno.csv", 
                                        delim=',', header=true))

# Number of phenotypes 
nPheno = convert(Int64, size(Y,2)/2) 
# Create Z matrix. The first column indicates treatment (wet/dry)
Z = hcat(repeat([1, -1], outer=nPheno), 
         kron(Matrix{Float64}(LinearAlgebra.I, nPheno, nPheno), 
              vcat([1 1], [1 -1])))

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
ZSub = hcat(repeat([1, -1], outer=nPhenoSub), 
            kron(Matrix{Float64}(LinearAlgebra.I, nPhenoSub, nPhenoSub), 
                 vcat([1 1], [1 -1])))

# Put together RawData object for MLM using subsetted data
MLMDataSub = RawData(Response(YSub), Predictors(X, ZSub))

# Dry run of L1-penalized matrix linear model using subset of phenotypes
resultsSub = mlmnet(fista_bt!, MLMDataSub, lambdas; isZInterceptReg=true) 


# Run L1-penalized matrix linear model while timing it
tic()
results = mlmnet(fista_bt!, MLMData, lambdas; isZInterceptReg=true) 
elapsedTime = toc()

# Flatten coefficients and write results to CSV
flat_coeffs = coef_2d(results)
CSV.write("../processed/lowry_l1_coeffs.csv", DataFrame(flat_coeffs))

# Print and write time to CSV
println(elapsedTime)
CSV.write("../processed/lowry_time.csv", DataFrame(elapsedTime))
