# L1-penalized matrix linear models
@everywhere include("../mlm_packages/matrixLMnet/src/matrixLMnet.jl")
@everywhere using matrixLMnet

# DataFrames 
using DataFrames

# JLD for saving variables
using JLD


# Read in Y (phenotypes). The first row is a header. The first column is IDs. 
Y = convert(Array{Float64}, readtable("./processed/agren_phe.csv", 
			                                separator = ',', header=true)[:,2:7])
# Take the log of Y
Y = log.(Y)
# Standardize Y 
Y = (Y.-mean(Y,1))./std(Y,1) 

# Read in X (genotype probabilities). The first row is a header. 
X = convert(Array{Float64}, readtable("./processed/agren_genoprobs.csv", 
                                      separator = ',', header=true))

# Create Z matrix. The first column indicates country (Italy/Sweden). 
Z = hcat([1, -1, 1, -1, 1, -1], eye(6))

# Put together RawData object for MLM 
MLM_data = RawData(Response(Y), Predictors(X, Z))

# Array of 50 lambdas
lambdas = reverse(1.2.^(-32:17))


# Run L1-penalized matrix linear model
results = mlmnet(fista_bt!, MLM_data, lambdas, isZInterceptReg=true)

# Flatten coefficients and write results to CSV
flat_coeffs = coef_2d(results)
writecsv("./processed/agren_l1_coeffs.csv", flat_coeffs)


# Run 10-fold cross-validation (on the rows)
srand(120)
mlmnet_cv_objs = mlmnet_cv(fista_bt!, MLM_data, lambdas, 10, 1; 
                           isZInterceptReg=true)
# Look at summary information from cross-validation
println(mlmnet_cv_summary(mlmnet_cv_objs))

# Save Mlmnet_cv object
save("./processed/agren_l1_cv.jld", "mlmnet_cv_objs", mlmnet_cv_objs)
