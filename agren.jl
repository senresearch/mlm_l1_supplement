# Source files and dependencies
@everywhere include("../mlm_packages/matrixLMnet/src/matrixLMnet.jl")
@everywhere using matrixLMnet

using DataFrames
using JLD

# Read in Y (pheno).
# Y is only the first 6 columns (the last one is ids). The first row is a header. 
Y = convert(Array, readtable("./processed/agren_phe.csv", separator = ',', header=true)[:,1:6])
# Drop missing rows of Y from X and Y. 
dropidx = vec(.!any(Y.=="-",2))
Ystd = Y[dropidx, :]
for i in 1:length(Ystd)
  if (typeof(Ystd[i]) == String)
    Ystd[i] = parse(Float64, Ystd[i])
  end
end
# Take the log of Y
Ystd = log.(convert(Array{Float64}, Ystd))
# Standardize Y using means/standard deviation
Ystd = (Ystd.-mean(Ystd,1))./std(Ystd,1) 

# Read in X (genotype probabilities)
X = readtable("./processed/agren_genoprobs.csv", separator = ',', header=true)
Xnoint = convert(Array{Float64}, X[dropidx, :])


# Create the Z contrast
# First column indicates country
# Last six columns form an identity matrix 
Znoint = hcat([1, 1, 1, -1, -1, -1], eye(6))






# Run L1 with this stuff. Penalize the Z intercept
MLM_data = RawData(Response(Ystd), Predictors(Xnoint, Znoint))
lambdas = reverse(1.2.^(-32:17))

srand(120)
mlmnet_cv_objs = mlmnet_cv(fista_bt!, MLM_data, lambdas, 8, 1; isZInterceptReg=true)
mlmnet_cv_summary(mlmnet_cv_objs)

save("./processed/agren_l1_cv.jld", "mlmnet_cv_objs", mlmnet_cv_objs)


results = mlmnet(fista_bt!, MLM_data, lambdas, isZInterceptReg=true)

flat_coeffs = coef_2d(results)
writecsv("./processed/agren_l1_coeffs.csv", flat_coeffs)












reps = 10
agren_times = SharedArray{Float64}(5, reps)

println("Starting")

# Dry run
mlmnet(fista_bt!, MLM_data, lambdas, isZInterceptReg=true)
# Get times from running FISTA with backtracking
@sync @parallel for j in 1:reps
  agren_times[5,j] = @elapsed mlmnet(fista_bt!, MLM_data, lambdas, isZInterceptReg=true)
end

# Cluster gives memory error when trying to calculate step size, so manually passing in
# Dry run
mlmnet(fista!, MLM_data, lambdas, isZInterceptReg=true)#, stepsize=0.00389033154159019)
# Get times from running FISTA with fixed step size
@sync @parallel for j in 1:reps
  agren_times[4,j] = @elapsed mlmnet(fista!, MLM_data, lambdas, isZInterceptReg=true)#, stepsize=0.00389033154159019)
end

# Cluster gives memory error when trying to calculate step size, so manually passing in
# Dry run
mlmnet(ista!, MLM_data, lambdas, isZInterceptReg=true)#, stepsize=0.00389033154159019)
# Get times from running ISTA with fixed step size
@sync @parallel for j in 1:reps
  agren_times[3,j] = @elapsed mlmnet(ista!, MLM_data, lambdas, isZInterceptReg=true)#, stepsize=0.00389033154159019)
end

# Dry run
mlmnet(cd_active!, MLM_data, lambdas, isZInterceptReg=true)
# Get times from running active coordinate descent
@sync @parallel for j in 1:reps
  agren_times[2,j] = @elapsed mlmnet(cd_active!, MLM_data, lambdas, isZInterceptReg=true)
end

# Dry run
mlmnet(cd!, MLM_data, lambdas, isZInterceptReg=true)
# Get times from running cyclic coordinate descent
@sync @parallel for j in 1:reps
  agren_times[1,j] = @elapsed mlmnet(cd!, MLM_data, lambdas, isZInterceptReg=true)
end

println(mean(agren_times, 2))

writecsv("./processed/agren_times.csv",  
          vcat(["method" "mean" transpose(collect(1:reps))], 
                hcat(["cd", "cd_active", "ista", "fista", "fista_bt"], 
                      mean(agren_times, 2), agren_times)))