using Distributed
using DataFrames
import Statistics.mean, Statistics.std
using Random
using CSV
using JLD2

# L1-penalized matrix linear models
@everywhere using matrixLMnet


# Read in Y (phenotypes). The first row is a header. The first column is IDs. 
Y = convert(Array{Float64, 2}, CSV.read("../processed/agren_phe.csv", 
			                            delim=',', header=true)[:,2:end])
# Take the log of Y
Y = log.(Y)
# Standardize Y 
Y = (Y.-mean(Y, dims=1)) ./ std(Y, dims=1) 

# Read in X (genotype probabilities). The first row is a header. 
X = convert(Array{Float64, 2}, CSV.read("../processed/agren_genoprobs.csv", 
                                        delim=',', header=true))

# Create Z matrix, indicating country (Italy/Sweden). 
Z = reshape([1, -1, 1, -1, 1, -1], 6, 1)

# Put together RawData object for MLM 
MLMData = RawData(Response(Y), Predictors(X, Z))

# Array of 50 lambdas
lambdas = reverse(1.2.^(-32:17))

# Run L1-penalized matrix linear model
results = mlmnet(fista_bt!, MLMData, lambdas, isZInterceptReg=true)

# Flatten coefficients and write results to CSV
flatCoeffs = coef_2d(results)
CSV.write("../processed/agren_l1_coeffs.csv", 
          DataFrame(flatCoeffs), writeheader=false)


# Run 10-fold cross-validation (on the rows)
Random.seed!(120)
mlmnetCVObjs = mlmnet_cv(fista_bt!, MLMData, lambdas, 10, 1; 
                         isZInterceptReg=true)
# Look at summary information from cross-validation
println(mlmnet_cv_summary(mlmnetCVObjs))

# Save Mlmnet_cv object
@save "../processed/agren_l1_cv.jld2" mlmnetCVObjs
