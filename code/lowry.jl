# L1-penalized matrix linear models
include("../../mlm_packages/matrixLMnet/src/matrixLMnet.jl")
using matrixLMnet

# DataFrames 
using DataFrames


# Read in X (genotype probabilities) with cytoplasm contrast. 
# The first row is a header. 
X = convert(Array{Float64}, readtable("../processed/lowry_cyto_genoprobs.csv", 
                                      separator=',', header=true))

# Read in Y (phenotypes). The first row is a header. 
Y = convert(Array{Float64}, readtable("../processed/lowry_pheno.csv", 
                                      separator=',', header=true))

# Number of phenotypes 
npheno = convert(Int64, size(Y,2)/2) 
# Create Z matrix. The first column indicates treatment (wet/dry)
Z = hcat(repeat([1, -1], outer=npheno), 
         kron(eye(npheno), vcat([1 1], [1 -1])))

# Put together RawData object for MLM
MLMData = RawData(Response(Y), Predictors(X, Z))

# Array of 16 lambdas
lambdas = reverse(1.2.^(-9:6))


# Randomly sample a small number of phenotypes for a dry run
randSamp = 500 # Number of random phenotypes to sample
srand(10)
randPhenos = sort(shuffle(1:npheno)[1:randSamp])
idx = Array(Int64, 2*randSamp)
for i in 1:randSamp
	idx[2*i-1] = randPhenos[i]*2-1
	idx[2*i] = randPhenos[i]*2
end

# Subset the Y matrix to keep only these phenotypes
YSub = Y[:,idx]

# Create the Z matrix for this subset of phenotypes
nPhenoSub = convert(Int64, size(YSub,2)/2) 
ZSub = hcat(repeat([1, -1], outer=nPhenoSub), 
            kron(eye(nPhenoSub), vcat([1 1], [1 -1])))

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
writecsv("../processed/lowry_l1_coeffs.csv", flat_coeffs)

# Print and write time to CSV
println(elapsedTime)
writecsv("../processed/lowry_time.csv", elapsedTime)
