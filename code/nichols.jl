using Distributed
using DataFrames
using Random
using CSV
using JLD2

# Matrix linear models for genetic screening data
@everywhere using GeneticScreen

# L1-penalized matrix linear models
@everywhere using matrixLMnet


# Number of permutations 
nPerms = 1000

# Array of 25 lambdas
lambdas = reverse(1.3.^(-18:6))


# Iterate through the six plates
for i in 1:6
    
    # Read in data for each plate
    # Colony opacity
    Y = CSV.read(string("../processed/processed_KEIO_data/p", i, 
                        "_krit_dat.csv"), DataFrame, delim=',', header=true) 
    
    # Conditions
    X = CSV.read(string("../processed/processed_KEIO_data/p", i, 
                        "_krit_cond.csv"), DataFrame, delim=',', header=true) 
    
    # Mutant keys
    Z = CSV.read(string("../data/raw_KEIO_data/KEIO", i, 
                        "_KEY.csv"), DataFrame, delim='\t', header=true) 
    
    # Put together RawData object for matrix linear models 
    MLMData = read_plate(X[:,[:Cond_Conc]], Y, Z[:,[:name]]; 
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
	mlmnetCVfolds = make_folds_conds(convert(Array{String}, X[:,:Cond_Conc]))
	# Run 10-fold cross-validation (on the rows)
    mlmnetCVObjs = mlmnet_cv(fista_bt!, MLMData, lambdas, mlmnetCVfolds, 1)
# Look at summary information from cross-validation for best lambdas
    println(lambda_min(mlmnetCVObjs))
    
    # Save Mlmnet_cv object
	@save string("../processed/nichols_p", i, "_l1_cv.jld2") mlmnetCVObjs
    
end
