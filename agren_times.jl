# L1-penalized matrix linear models
@everywhere include("../mlm_packages/matrixLMnet/src/matrixLMnet.jl")
@everywhere using matrixLMnet

# DataFrames 
using DataFrames


# Read in Y (phenotypes). The first row is a header. The first column is IDs. 
Y = convert(Array{Float64}, readtable("./processed/agren_phe.csv", 
			                          separator=',', header=true)[:,2:end])
# Take the log of Y
Y = log.(Y)
# Standardize Y 
Y = (Y.-mean(Y,1))./std(Y,1) 

# Read in X (genotype probabilities). The first row is a header. 
X = convert(Array{Float64}, readtable("./processed/agren_genoprobs.csv", 
                                      separator=',', header=true))

# Create Z matrix. The first column indicates country (Italy/Sweden). 
Z = hcat([1, -1, 1, -1, 1, -1], eye(6))

# Put together RawData object for MLM 
MLMData = RawData(Response(Y), Predictors(X, Z))

# Array of 50 lambdas
lambdas = reverse(1.2.^(-32:17))


# Number of replicates
reps = 10
# Initialize array for storing times
agrenTimes = SharedArray{Float64}(5, reps)


# Get times from running FISTA with backtracking
@sync @parallel for j in 1:reps
    agrenTimes[5,j] = @elapsed mlmnet(fista_bt!, MLMData, lambdas, 
                                      isZInterceptReg=true) 
end

# Get times from running FISTA with fixed step size
@sync @parallel for j in 1:reps
    agrenTimes[4,j] = @elapsed mlmnet(fista!, MLMData, lambdas, 
                                      isZInterceptReg=true)
end

# Get times from running ISTA with fixed step size
@sync @parallel for j in 1:reps
    agrenTimes[3,j] = @elapsed mlmnet(ista!, MLMData, lambdas, 
                                      isZInterceptReg=true)
end

# Get times from running active coordinate descent
@sync @parallel for j in 1:reps
    agrenTimes[2,j] = @elapsed mlmnet(cd_active!, MLMData, lambdas, 
                                      isZInterceptReg=true)
end

# Get times from running cyclic coordinate descent
@sync @parallel for j in 1:reps
    agrenTimes[1,j] = @elapsed mlmnet(cd_active!, MLMData, lambdas, 
                                      isZInterceptReg=true, isRandom=false)
end


# Print and write times to CSV
println(mean(agrenTimes, 2))
writecsv("./processed/agren_times.csv",  
         vcat(["method" "mean" transpose(collect(1:reps))], 
              hcat(["cd", "cd_active", "ista", "fista", "fista_bt"], 
                   mean(agrenTimes, 2), agrenTimes)))