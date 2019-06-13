using Distributed
using DataFrames
using Statistics
using LinearAlgebra
using Distributions
using Random
using CSV
using SharedArrays

# L1-penalized matrix linear models
@everywhere include("../../mlm_packages/matrixLMnet/src/matrixLMnet.jl")
@everywhere using Main.matrixLMnet


@everywhere begin

    """
        sim_effect(n, propNonzero; eDist)
    
    Simulate effects with a given proportion of nonzero effects drawn from some 
    distribution. The remaining effects will be set to zero.
    
    # Arguments 
    
    - n = length of 1d effect array. 
    - propNonzero = proportion of nonzero effects. Defaults to `1/2`.
    
    # Keyword arguments
    
    - eDist = distribution from which the nonzero effects should be simulated. 
      Defaults to Normal(0,2). 
    
    # Value
    
    1d array of floats 
    
    """
    function sim_effect(n::Int64, propNonzero::Float64=1/2, 
                        eDist::Distribution=Normal(0,2))
        
        # Initialize vector for storing effects 
        effect = zeros(n)
        
        # Randomly sample indices of nonzero effects
        idx = sample(1:n, convert(Integer, round(n*propNonzero)); 
                     replace=false)
        
        # Simulate and assign nonzero effects  
        effect[idx] = rand(eDist, convert(Integer, round(n*propNonzero)))
        
        return effect
    end
end


@everywhere begin

    """
        make_Y(n, m, fixed, edist)
    
    Simulate a 2d response array for simulations, with user-specified fixed 
    effects. 

    # Arguments 
    
    - n = number of rows
    - m = number of columns
    - fixed = 2d array of fixed effects (should be n by m)
    - edist = distribution from which the non-fixed effects should be randomly 
      sampled. Defaults to Normal(0,3). 
    
    # Value
    
    2d response array
    
    """
    function make_Y(n::Int64, m::Int64, fixed::Array{Float64,2}, 
                    edist::Distribution=Normal(0,3))
        
        return fixed + rand(edist, n, m)
    end
end 


@everywhere begin
    
    """
        sim_raw_data(n, m, p, q, seed)

    Simulate RawData object
    
    # Arguments 
    
    - n = number of rows of X = number of rows of Y
    - m = number of rows of Z = number of columns of Y
    - p = number of columns of X
    - q = number of columns of Z
    - seed = random seed
    
    # Value
    
    RawData object with simulated data
    
    # Additional notes
    
    Row and column main effects are 1/2 nonzero drawn from Normal(0,2). 
    Interactions are 1/8 nonzero drawn from Normal(0,2). Errors are drawn 
    from Normal(0,3). 
    
    """
    function sim_raw_data(n::Int64, m::Int64, p::Int64, q::Int64, 
                          seed::Int64=10)
        
        # Set random seed
        Random.seed!(seed)
        
        # Simulate row main effects, 1/2 nonzero from Normal(0,2). 
        d = sim_effect(p) 
        # Simulate column main effects, 1/2 nonzero from Normal(0,2). 
        c = sim_effect(q) 
        # Simulate interaction effects, 1/8 nonzero from Normal(0,2). 
        inter = reshape(sim_effect((p)*(q), 1/8), p, q) 
        
        # Simulate X and Z
        X = repeat(Matrix{Float64}(LinearAlgebra.I, p, p), 
                   outer=[Int(ceil(n/p)),1])[1:n,:]
        Z = repeat(Matrix{Float64}(LinearAlgebra.I, q, q), 
                   outer=[Int(ceil(m/q)),1])[1:m,:]
        
        # Generate the fixed effects
        fixed = X*d .+ transpose(Z*c) .+ X*inter*transpose(Z)
        # Simulate Y using fixed effects 
        Y = make_Y(n, m, fixed)
        # Standardize Y
        Y = (Y.-Statistics.mean(Y, dims=1)) ./ Statistics.std(Y, dims=1) 
        	
        # Put together RawData object for MLM
        return RawData(Response(Y), Predictors(X, Z))
    end
end


@everywhere begin
    
    """
        run_sim(lambdas, fun; n, m, p, q, seed)
    
    Simulate RawData object and run L1-penalized matrix linear model

    # Arguments 
    
    - lambdas = 1d array of floats consisting of lambda penalties in 
	descending order
    - fun = function that applies an L1-penalty estimate method. Defaults to 
    fista_bt!
    - n = number of rows of X = number of rows of Y
    - m = number of rows of Z = number of columns of Y
    - p = number of columns of X
    - q = number of columns of Z
    - seed = random seed
    - funArgs = variable keyword arguments to be passed into `fun`

    # Value
    
    Mlmnet object
    
    # Additional notes
    
    Row and column main effects are 1/2 nonzero drawn from Normal(0,2). 
    Interactions are 1/8 nonzero drawn from Normal(0,2). Errors are drawn 
    from Normal(0,3). Default parameters are used for the call to mlmnet. 

    """
    function run_sim(lambdas::Array{Float64,1}, fun::Function=fista_bt!; 
                     n::Int64=600, m::Int64=600, 
                     p::Int64=200, q::Int64=200, 
                     seed::Int64=10, funArgs...)
        
        # Simulate data
        MLMData = sim_raw_data(n, m, p, q, seed)
        
        # Run L1-penalized matrix linear model
        return mlmnet(fun, MLMData, lambdas; funArgs...)
    end
end


# Array of 50 lambdas
lambdas = reverse(1.2.^(-32:17))
# Number of replicates
reps = 10

# Range of p and q values to try
pqVals = collect(200:200:1000)
# Range of n and m values to try 
nmVals = collect(400:400:2000)


# Generate grid of p and q values
pqGrid = vec(collect(Base.product(pqVals, pqVals)))
# Initialize array for storing times
pqTimes = SharedArrays.SharedArray{Float64}(length(pqGrid), reps)
# Hold n and m fixed at 1200 and vary p and q over a grid
@sync @distributed for j in 1:reps
    for i in 1:length(pqGrid)
        pqTimes[i,j] = @elapsed run_sim(lambdas; n=Int64(mean(nmVals)), 
                                        m=Int64(mean(nmVals)), 
                                        p=pqGrid[i][1], 
                                        q=pqGrid[i][2], stepsize=1.0)
    end
end

# Print and write times to CSV
println(reshape(Statistics.mean(pqTimes, dims=2), 
                length(pqVals), length(pqVals)))
CSV.write("../processed/pq_times.csv",  
          DataFrame(vcat(["p" "q" "mean" transpose(collect(1:reps))], 
                         hcat([x[1] for x in pqGrid], [x[2] for x in pqGrid], 
                              Statistics.mean(pqTimes, dims=2), pqTimes))), 
          writeheader=false)


# Generate grid of n and m values
nmGrid = vec(collect(Base.product(nmVals, nmVals)))
# Initialize array for storing times
nmTimes = SharedArrays.SharedArray{Float64}(length(nmGrid), reps)
# Hold p and q fixed at 600 and vary n and m over a grid
@sync @distributed for j in 1:reps
    for i in 1:length(nmGrid)
        nmTimes[i,j] = @elapsed run_sim(lambdas; n=nmGrid[i][1], 
                                        m=nmGrid[i][2], 
                                        p=Int64(mean(pqVals)), 
                                        q=Int64(mean(pqVals)), stepsize=1.0)
    end
end

# Print and write times to CSV
println(reshape(Statistics.mean(nmTimes, dims=2), 
                length(nmVals), length(nmVals)))
CSV.write("../processed/nm_times.csv",  
          DataFrame(vcat(["n" "m" "mean" transpose(collect(1:reps))], 
                         hcat([x[1] for x in nmGrid], [x[2] for x in nmGrid], 
                              Statistics.mean(nmTimes, dims=2), nmTimes))), 
          writeheader=false)
