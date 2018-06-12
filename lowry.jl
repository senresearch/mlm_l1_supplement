# Source files and dependencies
include("../mlm_packages/matrixLMnet/src/matrixLMnet.jl")
using matrixLMnet

using DataFrames
include("dummyfun.jl")

#include("l1_streamlined.jl")
#include("FISTA_backtrack.jl")


# Read in X (genotype probabilities) with cytoplasm contrast. The first row is a header. 
X = convert(Array{Float64}, readtable("./processed/juenger_cyto_genoprobs.csv", separator = ',', header=true))

# Read in Y (phenotypes). The first row is a header. 
Y = convert(Array, readtable("./processed/juenger_pheno.csv", separator = ',', header=true))

# Create the Z contrast. 
npheno = convert(Int64, size(Y,2)/2) 
Z = hcat(repeat([1, -1], outer=npheno), kron(eye(npheno), vcat([1 1], [1 -1])))

MLM_data = RawData(Response(Y), Predictors(X, Z))

# 16 lambda values
lambdas = reverse(1.2.^(-9:6))

# Real data
println("Starting")
#@time reals = run_l1(fista!, X, Y, Z, lambdas; filename="./processed/juenger_l1_real.csv", isZInterceptReg=true)

# Dry run 
results = mlmnet(fista_bt!, MLM_data, lambdas; isZInterceptReg=true; maxiter=5)







tic()
results = mlmnet(fista_bt!, MLM_data, lambdas; isZInterceptReg=true)
elapsed_time = toc()


flat_coeffs = coef_2d(results)
writecsv("./processed/lowry_l1_coeffs.csv", flat_coeffs)

writecsv("./processed/lowry_time.csv", elapsed_time)
