# L1-penalized matrix linear models
include("../mlm_packages/matrixLMnet/src/matrixLMnet.jl")
using matrixLMnet

# DataFrames 
using DataFrames

# Functions for creating dummy contrasts for categorical variables
include("dummy_fun.jl")


# Read in X (genotype probabilities) with cytoplasm contrast. 
# The first row is a header. 
X = convert(Array{Float64}, readtable("./processed/lowry_cyto_genoprobs.csv", 
                                      separator = ',', header=true))

# Read in Y (phenotypes). The first row is a header. 
Y = convert(Array, readtable("./processed/lowry_pheno.csv", 
                             separator = ',', header=true))

# Number of phenotypes 
npheno = convert(Int64, size(Y,2)/2) 
# Create Z matrix 
Z = hcat(repeat([1, -1], outer=npheno), # wet/dry
         kron(eye(npheno), vcat([1 1], [1 -1])))

# Put together RawData object for MLM
MLM_data = RawData(Response(Y), Predictors(X, Z))

# Array of 16 lambdas
lambdas = reverse(1.2.^(-9:6))


# Randomly sample a small number of phenotypes for a dry run
randsamp = 500 # Number of random phenotypes to sample
srand(10)
rand_phenos = sort(shuffle(1:npheno)[1:randsamp])
idx = Array(Int64, 2*randsamp)
for i in 1:randsamp
	idx[2*i-1] = rand_phenos[i]*2-1
	idx[2*i] = rand_phenos[i]*2
end

# Subset the Y matrix to keep only these phenotypes
Ysub = Y[:,idx]

# Create the Z matrix for this subset of phenotypes
npheno_sub = convert(Int64, size(Ysub,2)/2) 
Zsub = hcat(repeat([1, -1], outer=npheno_sub), 
            kron(eye(npheno_sub), vcat([1 1], [1 -1])))

# Put together RawData object for MLM using subsetted data
MLM_data_sub = RawData(Response(Ysub), Predictors(X, Zsub))

# Dry run of L1-penalized matrix linear model using subset of phenotypes
results_sub = mlmnet(fista_bt!, MLM_data_sub, lambdas; 
                     isZInterceptReg=true, stepsize=0.01)


# Run L1-penalized matrix linear model while timing it
tic()
results = mlmnet(fista_bt!, MLM_data, lambdas; 
                 isZInterceptReg=true, stepsize=0.01)
elapsed_time = toc()

# Flatten coefficients and write results to CSV
flat_coeffs = coef_2d(results)
writecsv("./processed/lowry_l1_coeffs.csv", flat_coeffs)

# Print and write time to CSV
println(elapsed_time)
writecsv("./processed/lowry_time.csv", elapsed_time)
