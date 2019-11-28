using Distributed
using DataFrames
import Statistics.mean, Statistics.std
using CSV
using SharedArrays

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


# Number of replicates
reps = 100
# Initialize array for storing times
agrenTimes = SharedArrays.SharedArray{Float64}(6, reps)


# Get times from running ADMM
@sync @distributed for j in 1:reps
    agrenTimes[6,j] = @elapsed mlmnet(admm!, MLMData, lambdas, 
                                      isZInterceptReg=true, rho=0.005) 
end

# Get times from running FISTA with backtracking
@sync @distributed for j in 1:reps
    agrenTimes[5,j] = @elapsed mlmnet(fista_bt!, MLMData, lambdas, 
                                      isZInterceptReg=true) 
end

# Get times from running FISTA with fixed step size
@sync @distributed for j in 1:reps
    agrenTimes[4,j] = @elapsed mlmnet(fista!, MLMData, lambdas, 
                                      isZInterceptReg=true)
end

# Get times from running ISTA with fixed step size
@sync @distributed for j in 1:reps
    agrenTimes[3,j] = @elapsed mlmnet(ista!, MLMData, lambdas, 
                                      isZInterceptReg=true)
end

# Get times from running random coordinate descent
@sync @distributed for j in 1:reps
    agrenTimes[2,j] = @elapsed mlmnet(cd_active!, MLMData, lambdas, 
                                      isZInterceptReg=true)
end

# Get times from running cyclic coordinate descent
@sync @distributed for j in 1:reps
    agrenTimes[1,j] = @elapsed mlmnet(cd_active!, MLMData, lambdas, 
                                      isZInterceptReg=true, isRandom=false)
end


# Print and write times to CSV
println(mean(agrenTimes, dims=2))
CSV.write("../processed/agren_times.csv",  
          DataFrame(vcat(["method" "mean" transpose(collect(1:reps))], 
                    hcat(["cd_cyclic", "cd_random", "ista", 
                          "fista", "fista_bt", "admm"], 
                         mean(agrenTimes, dims=2), agrenTimes))), 
          writeheader=false)
