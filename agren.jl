# L1-penalized matrix linear models
@everywhere include("../mlm_packages/matrixLMnet/src/matrixLMnet.jl")
@everywhere using matrixLMnet

# DataFrames 
using DataFrames
# JLD for saving variables
using JLD


# Read in Y (phenotypes). The first row is a header. The last column is IDs. 
# Y = convert(Array, readtable("./processed/agren_phe.csv", 
			# separator = ',', header=true)[:,1:6])
# first column is IDs
Y = convert(Array, readtable("./processed/agren_phe.csv", 
			separator = ',', header=true)[:,2:7])
# # Drop missing rows of Y.  
# dropidx = vec(.!any(Y.=="-",2))
# Ystd = Y[dropidx, :]
# for i in 1:length(Ystd)
  # if (typeof(Ystd[i]) == String)
    # Ystd[i] = parse(Float64, Ystd[i])
  # end
# end
# Take the log of Y
Ystd = log.(convert(Array{Float64}, Ystd))
# Standardize Y 
Ystd = (Ystd.-mean(Ystd,1))./std(Ystd,1) 

# Read in X (genotype probabilities). The first row is a header. 
X = readtable("./processed/agren_genoprobs.csv", 
              separator = ',', header=true)
# # Drop missing rows of Y from X. 
# Xnoint = convert(Array{Float64}, X[dropidx, :])
Xnoint = convert(Array{Float64}, X)

# Create Z matrix. The first column indicates country (Italy/Sweden). 
# Znoint = hcat([1, 1, 1, -1, -1, -1], eye(6))
Znoint = hcat([1, -1, 1, -1, 1, -1], eye(6))

# Put together RawData object for MLM 
MLM_data = RawData(Response(Ystd), Predictors(Xnoint, Znoint))

# Array of 50 lambdas
lambdas = reverse(1.2.^(-32:17))


# Run L1-penalized matrix linear model
results = mlmnet(fista_bt!, MLM_data, lambdas, isZInterceptReg=true, 
				 stepsize=0.01)

# Flatten coefficients and write results to CSV
flat_coeffs = coef_2d(results)
writecsv("./processed/agren_l1_coeffs.csv", flat_coeffs)


# Run 8-fold cross-validation (on the rows)
srand(120)
mlmnet_cv_objs = mlmnet_cv(fista_bt!, MLM_data, lambdas, 8, 1; 
                           isZInterceptReg=true, stepsize=0.01)
# Look at summary information from cross-validation
println(mlmnet_cv_summary(mlmnet_cv_objs))

# Save Mlmnet_cv object
save("./processed/agren_l1_cv.jld", "mlmnet_cv_objs", mlmnet_cv_objs)


##############################################################################


# Number of replicates
reps = 15
# Initialize array for storing times
agren_times = SharedArray{Float64}(5, reps)


# Get times from running FISTA with backtracking
@sync @parallel for j in 1:reps
  agren_times[5,j] = @elapsed mlmnet(fista_bt!, MLM_data, lambdas, isZInterceptReg=true, stepsize=0.01)
end

# Get times from running FISTA with fixed step size
@sync @parallel for j in 1:reps
  agren_times[4,j] = @elapsed mlmnet(fista!, MLM_data, lambdas, isZInterceptReg=true)
end

# Get times from running ISTA with fixed step size
@sync @parallel for j in 1:reps
  agren_times[3,j] = @elapsed mlmnet(ista!, MLM_data, lambdas, isZInterceptReg=true)
end

# Get times from running active coordinate descent
@sync @parallel for j in 1:reps
  agren_times[2,j] = @elapsed mlmnet(cd_active!, MLM_data, lambdas, isZInterceptReg=true)
end

# Get times from running cyclic coordinate descent
@sync @parallel for j in 1:reps
  agren_times[1,j] = @elapsed mlmnet(cd!, MLM_data, lambdas, isZInterceptReg=true)
end


# Print and write times to CSV
println(mean(agren_times, 2))
writecsv("./processed/agren_times.csv",  
          vcat(["method" "mean" transpose(collect(1:reps))], 
                hcat(["cd", "cd_active", "ista", "fista", "fista_bt"], 
                      mean(agren_times, 2), agren_times)))