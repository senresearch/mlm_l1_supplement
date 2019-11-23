using Distributed
using DataFrames
using Statistics
using Random
using CSV
using JLD2

# Matrix linear models for genetic screening data
@everywhere using GeneticScreen

# L1-penalized matrix linear models
@everywhere using matrixLMnet


# Number of permutations 
nPerms = 1000

# Array of 40 lambdas
lambdas = reverse(1.2.^(-30:9))

# Iterate through the six plates
for i in 1:6
    
    # Read in data for each plate
    # Colony opacity
    Y = CSV.read(string("../processed/processed_KEIO_data/p", i, 
                        "_krit_dat.csv"), delim=',', header=true) 
    
    # Conditions
    X = CSV.read(string("../processed/processed_KEIO_data/p", i, 
                        "_krit_cond.csv"), delim=',', header=true) 
    
    # Mutant keys
    Z = CSV.read(string("../data/raw_KEIO_data/KEIO", i, 
                        "_KEY.csv"), delim='\t', header=true) 
    
    # Put together RawData object for matrix linear models 
    MLMData = read_plate(X[[:Cond_Conc]], Y, Z[[:name]]; 
                            XCVar=:Cond_Conc, ZCVar=:name,
                            XCType="sum", ZCType="sum", isYstd=true)
    
    # Run L1-penalized matrix linear model
    results = mlmnet(fista_bt!, MLMData, lambdas)
    
    # Flatten coefficients and write results to CSV
    flat_coeffs = coef_2d(results)
    CSV.write(string("../processed/nichols_p", i, "_l1_coeffs.csv"), 
              DataFrame(flat_coeffs), writeheader=false)

    # Create 10 folds (on the rows) that ensure that at least one condition 
	# replicate is included in each fold
    Random.seed!(120)
	mlmnetCVfolds = make_folds_conds(convert(Array{String}, X[:Cond_Conc]))
	# Run 10-fold cross-validation (on the rows)
    mlmnetCVObjs = mlmnet_cv(fista_bt!, MLMData, lambdas, mlmnetCVfolds, 1)
    # Look at summary information from cross-validation
    println(mlmnet_cv_summary(mlmnetCVObjs))
    
    # Save Mlmnet_cv object
	@save string("../processed/nichols_p", i, "_l1_cv.jld2") mlmnetCVObjs
    
end

function valid_reduce2(A::Array{Float64,2}, fun::Function=mean)
	out = Array{Float64}(undef, size(A, 1))
	for i in 1:size(A, 1)
		out[i] = fun(A[i, A[i, :] .< Inf])
	end
	return out
end


# Search for the "best" lambda (1 SE more than the lowest MSE)
for i in 1:6
    # Load CV output
    @load string("../processed/nichols_p", i, "_l1_cv.jld2") mlmnetCVObjs
    # Find the index of the lambda with the lowest MSE
    lambdaMinIdx = argmin(mean(mlmnetCVObjs.mse, dims=2))
    # Compute standard error across folds for the lowest MSE
    mse1se = std(mlmnetCVObjs.mse, dims=2)[lambdaMinIdx] + 
             mean(mlmnetCVObjs.mse, dims=2)[lambdaMinIdx]
    # Find the index of the lambda that is closest to being 1 SE greater than 
    # the lowest lambda, in the direction of the bigger lambdas
    lambdaMin1seIdx = argmin(abs.(mean(mlmnetCVObjs.mse, 
                                  dims=2)[1:lambdaMinIdx[1]].-mse1se))
    # Print the summary info for that lambda
    println(mlmnet_cv_summary(mlmnetCVObjs)[lambdaMin1seIdx,:])
end