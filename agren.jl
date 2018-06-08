# Source files and dependencies
@everywhere include("../mlm_packages/matrixLMnet/src/matrixLMnet.jl")
@everywhere using matrixLMnet

using DataFrames
using CSV
#using StatsBase
using MLBase
#using Gadfly
using JLD

#@everywhere include("FISTA_backtrack.jl")
#@everywhere include("FISTA.jl")


# Read in Y (pheno).
# Y is only the first 6 columns (the last one is ids). The first row is a header. 
Y = convert(Array, readtable("./processed/agren_phe.csv", separator = ',', header=true)[:,1:6])
# Drop missing rows of Y from X and Y. 
dropidx = vec(!any(Y.=="-",2))
Ystd = Y[dropidx, :]
for i in 1:length(Ystd)
  if (typeof(Ystd[i]) == String)
    Ystd[i] = parse(Float64, Ystd[i])
  end
end
# Take the log of Y
Ystd = log(convert(Array{Float64}, Ystd))
# Standardize Y using means/standard deviation
Ystd = (Ystd.-mean(Ystd,1))./std(Ystd,1) 

# Read in X (genotype probabilities)
X = readtable("./processed/agren_genoprobs.csv", separator = ',', header=true)
Xnoint = convert(Array{Float64}, X[dropidx, :])


# Create the Z contrast
# First column indicates country
# Last six columns form an identity matrix 
Znoint = hcat([1, 1, 1, -1, -1, -1], eye(6))




# @everywhere function fista_backtrack_pathwise(X::Array{Float64,2}, Y::Array{Float64,2}, Z::Array{Float64,2}, lambdas::Array{Float64,1}, fold::Int64, 
#                    isXIntercept::Bool=true, isZIntercept::Bool=true, isXReg::BitArray{1}=trues(size(X,2)), isZReg::BitArray{1}=trues(size(Z,2)),
#                    isStandardize::Bool=true, thresh::Float64=10.0^(-7), maxiter::Int=10^10)
#   # Check that the lambdas are in descending order. 
#   if any(lambdas .!= sort(lambdas, rev=true))
#     println("Sorting lambdas into descending order.")
#     lambdas = sort(lambdas, rev=true)
#   end

#   # Pre-allocate array for coefficients
#   coeffs = Array(Float64, (size(X,2)+isXIntercept)*(size(Z,2)+isZIntercept), length(lambdas)) 

#   # Start with coefficients initalized at zero for the largest lambda value
#   startB = zeros((size(X,2)+isXIntercept),(size(Z,2)+isZIntercept))

#   # Run FISTA for largest lambda value using the start coefficients of all zeros
#   tic()
#   coeffs[:,1] = vec(fista_backtrack(X, Y, Z, lambdas[1], startB, isZInterceptReg=true))
#   writecsv(string("./processed/agren_lambda",lambdas[1],"_fold", fold,".csv"), reshape(coeffs[:,1], size(X,2)+isXIntercept, size(Z,2)+isZIntercept))
#   toc()

#   # Iterate through the remaining lambdas
#   for i=2:length(lambdas) 
#     tic()
#     coeffs[:,i] = vec(fista_backtrack(X, Y, Z, lambdas[i], reshape(coeffs[:,i-1], size(X,2)+isXIntercept, size(Z,2)+isZIntercept), isZInterceptReg=true))
#     writecsv(string("./processed/agren_lambda",lambdas[i],"_fold", fold,".csv"), reshape(coeffs[:,i], size(X,2)+isXIntercept, size(Z,2)+isZIntercept))
#     toc()
#   end

#   return coeffs
# end


# Run L1 with this stuff. Penalize the Z intercept
srand(120)
nfolds = 8

# Provide headers where you just pass in the number of row/column folds and then they get generated.
mlmnet_cv(fun::Function, data::RawData, lambdas::Array{Float64,1}, rowFolds::Array{Array,1}, colFolds::Array{Array,1}; 
    isXIntercept::Bool=true, isZIntercept::Bool=true, 
    isXReg::BitArray{1}=trues(size(X,2)), isZReg::BitArray{1}=trues(size(Z,2)), 
    isXInterceptReg::Bool=false, isZInterceptReg::Bool=false, 
    isStandardize::Bool=true, stepsize::Float64=1.0, funArgs...)

save("./processed/agren_l1_cv.jld", "my_mlm_cv_obj", my_mlm_cv_obj)



#### general routine for generating folds ####

# 0. Base function that gets called by everything else 
# Pass in something like collect(Kfold(size(Xnoint,1), 8)) for rowFolds, colFolds
# mlmnet_cv(fun::Function, data::RawData, lambdas::Array{Float64,1}, 
#              rowFolds::Array{Array,1}, colFolds::Array{Array,1}; ...)

# 1. mlmnet_cv(fun::Function, data::RawData, lambdas::Array{Float64,1}, 
#              nRowFolds::Int64=10, nColFolds::Int64=nRowFolds; ...)
# Assume that if nFolds = 1, that means all the data will be used

# 2. mlmnet_cv(fun::Function, data::RawData, lambdas::Array{Float64,1}, 
#              rowFolds::Array{Array,1}, nColFolds::Int64=10; ...)
# Generate only column folds

# 3. mlmnet_cv(fun::Function, data::RawData, lambdas::Array{Float64,1}, 
#              nRowFolds::Int64=10, colFolds::Array{Array,1}; ...)
# Generate only row folds

# add get_rest function (returns array of booleans) to helpers 
# if there was only 1 fold along one dimension, the whole array needs to be passed into the MSE function
# incorporate into get_rest function so that it returns an array of all trues if the fold=everything


#### 


X_folds_idx = collect(Kfold(size(Xnoint,1), 8))


# Write the X, Y, and Z arrays from each fold to CSV.
# Generate the X, Y, and Z matrices for each fold
X_folds = []
Y_folds = []
for i = 1:nfolds
  push!(X_folds, Xnoint[X_folds_idx[i],:]) 
  push!(Y_folds, Ystd[X_folds_idx[i],:])
end 
Z_folds = repeat([Znoint], outer=nfolds)

# Do 8-fold CV for the 50 lambda values using overlapping folds.
lambdas = reverse(1.2.^(-32:17))[1:10]
# lambdas = [1.0, 0.5, 0.1]


# #@time pmap(l1_pathwise, repeat([ista], outer=nfolds), X_folds, Y_folds, Z_folds, repeat([lambdas], outer=nfolds))
# #a = @elapsed pmap(fista_backtrack_pathwise, X_folds, Y_folds, Z_folds, repeat([lambdas], outer=nfolds), collect(1:1:nfolds))
# #println(a)
# @time pmap(fista_backtrack_pathwise, X_folds, Y_folds, Z_folds, repeat([lambdas], outer=nfolds), collect(1:1:nfolds))


# function getrest(fold, allidx)
#   out = Array(Bool, length(allidx))
#   for i = 1:length(out)
#     out[i] = allidx[i] in fold
#   end
#   return(!out)
# end


# """
#   calc_mse(lambda, X_folds_idx, Z_folds_idx, X, Y, Z, isStandardize)

# Calculates three types of MSE for each of the CV folds for a given lambda

# # Arguments 

# - lambda = lambda penalty, a floating scalar
# - X_folds_idx = 2d boolean array generated by make_folds for the rows of X 
# - Z_folds_idx = 2d boolean array generated by make_folds for the rows of Z
# - X, Y, Z = 2d arrays X, Y, Z standardized and coded the way they were when being passed into the call to FISTA
# - isXIntercept = boolean flag indicating whether or not to include an X intercept (row main effects)
#     Defaults to true
# - isZIntercept = boolean flag indicating whether or not to include a Z intercept (column main effects)
#     Defaults to true

# # Value

# 2d array of floats with dimensions 3 by the number of folds 

# # Notes 

# The three types of MSE are calculated based on the residuals corresponding to the observations with 
# 1. X covariates not kept in the fold and Z covariates kept in the fold 
# 2. X covariates kept in the fold and Z covariates not kept in the fold 
# 3. X covariates not kept in the fold and Z covariates not kept in the fold 

# """
# function calc_mse(lambda, X_folds_idx, X, Y, Z, isXIntercept=true, isZIntercept=true)
#   # Include X and Z intercepts if necessary
#   if isXIntercept == true
#     X = hcat(ones(size(X,1)), X)
#   end
#   if isZIntercept == true
#     Z = hcat(ones(size(Z,1)), Z)
#   end

#   # Number of folds 
#   nfolds = length(X_folds_idx)

#   # Initialize array to store MSE values 
#   MSE = Array(Float64, 3, nfolds)

#   # For each fold 
#   for j in 1:nfolds
#     coeffs = readcsv(string("./processed/agren_lambda",lambda,"_fold", j,".csv")) # Read in the coefficients for that lambda
#     resid = Y - X*coeffs*transpose(Z) # Find the residuals 

#     # Calculate each of the MSE values 
#     MSE[1,j] = mean(resid[getrest(X_folds_idx[j], 1:size(X,1)), :].^2)
#     MSE[2,j] = mean(coeffs[2:end,2:end].==0)
#     MSE[3,j] = lambda
#   end

#   return MSE
# end

# MSE_arr = Array(Float64, length(lambdas), 3, nfolds)

# for k in 1:length(lambdas)
#   MSE_arr[k,:,:] = calc_mse(lambdas[k], X_folds_idx, Xnoint, Ystd, Znoint)
# end

# MSE_df = DataFrame(mean(MSE_arr, 3)[:,:,1])
# names!(MSE_df, map(parse, ["AvgMSE", "PercentZero", "Lambda"]))
# println(MSE_df)


####### used for plotting ###########
# out = fista_backtrack_pathwise(Xnoint, Ystd, Znoint, lambdas) 
#writecsv("./processed/agren_lambdas_optim.csv", out)

MLM_data = RawData(Response(Ystd), Predictors(Xnoint, Znoint))
results = mlmnet(fista_bt!, MLM_data, lambdas, isZInterceptReg=true)

flat_coeffs = Array(Float64, size(results.B,2)*size(results.B,3), length(lambdas))
for i=1:length(lambdas)
  flat_coeffs[:,i] = vec(results.B[i,:,:])
end
writecsv("./processed/agren_l1_coeffs.csv", flat_coeffs)












# @everywhere include("../mlm_packages/matrixLMnet/src/matrixLMnet.jl")
# @everywhere using matrixLMnet

# using DataFrames

# # Read in Y (pheno).
# # Y is only the first 6 columns (the last one is ids). The first row is a header. 
# Y = convert(Array, readtable("./processed/agren_phe.csv", separator = ',', header=true)[:,1:6])
# # Drop missing rows of Y from X and Y. 
# dropidx = vec(!any(Y.=="-",2))
# Ystd = Y[dropidx, :]
# for i in 1:length(Ystd)
#   if (typeof(Ystd[i]) == String)
#     Ystd[i] = parse(Float64, Ystd[i])
#   end
# end
# # Take the log of Y
# Ystd = log(convert(Array{Float64}, Ystd))
# # Standardize Y using means/standard deviation
# Ystd = (Ystd.-mean(Ystd,1))./std(Ystd,1) 

# # Read in X (genotype probabilities)
# X = readtable("./processed/agren_genoprobs.csv", separator = ',', header=true)
# Xnoint = convert(Array{Float64}, X[dropidx, :])


# # Create the Z contrast
# # First column indicates country
# # Last six columns form an identity matrix 
# Znoint = hcat([1, 1, 1, -1, -1, -1], eye(6))


# MLM_data = RawData(Response(Ystd), Predictors(Xnoint, Znoint))

# # 50 lambda values
# lambdas = reverse(1.2.^(-32:17))

reps = 10
agren_times = SharedArray{Float64}(5, reps)

println("Starting")

# Dry run
mlmnet(fista_bt!, MLM_data, lambdas, isZInterceptReg=true)
# Get times from running FISTA with backtracking
@sync @parallel for j in 1:reps
  agren_times[5,j] = @elapsed mlmnet(fista_bt!, MLM_data, lambdas, isZInterceptReg=true)
end
println(agren_times[5,:])

# Cluster gives memory error when trying to calculate step size, so manually passing in
# Dry run
mlmnet(fista!, MLM_data, lambdas, isZInterceptReg=true)#, stepsize=0.00389033154159019)
# Get times from running FISTA with fixed step size
@sync @parallel for j in 1:reps
  agren_times[4,j] = @elapsed mlmnet(fista!, MLM_data, lambdas, isZInterceptReg=true)#, stepsize=0.00389033154159019)
end
println(agren_times[4,:])

# Cluster gives memory error when trying to calculate step size, so manually passing in
# Dry run
mlmnet(ista!, MLM_data, lambdas, isZInterceptReg=true)#, stepsize=0.00389033154159019)
# Get times from running ISTA with fixed step size
@sync @parallel for j in 1:reps
  agren_times[3,j] = @elapsed mlmnet(ista!, MLM_data, lambdas, isZInterceptReg=true)#, stepsize=0.00389033154159019)
end
println(agren_times[3,:])

# Dry run
mlmnet(cd_active!, MLM_data, lambdas, isZInterceptReg=true)
# Get times from running active coordinate descent
@sync @parallel for j in 1:reps
  agren_times[2,j] = @elapsed mlmnet(cd_active!, MLM_data, lambdas, isZInterceptReg=true)
end
println(agren_times[2,:])

# Dry run
mlmnet(cd!, MLM_data, lambdas, isZInterceptReg=true)
# Get times from running cyclic coordinate descent
@sync @parallel for j in 1:reps
  agren_times[1,j] = @elapsed mlmnet(cd!, MLM_data, lambdas, isZInterceptReg=true)
end
println(agren_times[1,:])

println(mean(agren_times, 2))

writecsv("./processed/agren_times.csv", agren_times)
writecsv("./processed/agren_times.csv", hcat(mean(agren_times, 2), agren_times))