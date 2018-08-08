# L1-penalized matrix linear models
@everywhere include("../../mlm_packages/matrixLMnet/src/matrixLMnet.jl")
@everywhere using matrixLMnet


@everywhere begin
    
    """
        simRawData(n, m, p, q, seed)

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
    
    function simRawData(n::Int64, m::Int64, p::Int64, q::Int64, 
                        seed::Int64=10)
        # Set random seed
        srand(seed)
        
        # Simulate row main effects, 1/2 nonzero from Normal(0,2). 
        d = make_effect(p) 
        # Simulate column main effects, 1/2 nonzero from Normal(0,2). 
        c = make_effect(q) 
        # Simulate interaction effects, 1/8 nonzero from Normal(0,2). 
        inter = reshape(make_effect((p)*(q), 1/8), p, q) 
        
        # Simulate X and Z
        X = repeat(eye(p), outer=[Int(ceil(n/p)),1])[1:n,:]
        Z = repeat(eye(q), outer=[Int(ceil(m/q)),1])[1:m,:]
        
        # Generate the fixed effects
        fixed = X*d .+ transpose(Z*c) .+ A_mul_Bt(X*inter, Z)
        # Simulate Y using fixed effects 
        Y = make_Y(n, m, fixed)
        # Standardize Y
        Y = (Y.-mean(Y,1))./std(Y,1) 
        	
        # Put together RawData object for MLM
        return RawData(Response(Y), Predictors(X, Z))
    end
end


@everywhere begin
    
    """
        runSim(lambdas, fun; n, m, p, q, seed)
    
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
    
    function runSim(lambdas::Array{Float64,1}, fun::Function=fista_bt!; 
                    n::Int64=600, m::Int64=600, 
                    p::Int64=200, q::Int64=200, seed::Int64=10, funArgs...)
        # Simulate data
        MLMData = simRawData(n, m, p, q, seed)
        
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
pqTimes =  SharedArray{Float64}(length(pqGrid), reps)
# Hold n and m fixed at 1200 and vary p and q over a grid
@sync @parallel for j in 1:reps
    for i in 1:length(pqGrid)
        pqTimes[i,j] = @elapsed runSim(lambdas; n=Int64(mean(nmVals)), 
                                       m=Int64(mean(nmVals)), 
                                       p=pqGrid[i][1], 
                                       q=pqGrid[i][2], stepsize=1.0)
    end
end

# Print and write times to CSV
println(reshape(mean(pqTimes, 2), length(pqVals), length(pqVals)))
writecsv("../processed/pq_times.csv",  
          vcat(["p" "q" "mean" transpose(collect(1:reps))], 
               hcat([x[1] for x in pqGrid], [x[2] for x in pqGrid],
                    mean(pqTimes, 2), pqTimes)))


# Generate grid of n and m values
nmGrid = vec(collect(Base.product(nmVals, nmVals)))
# Initialize array for storing times
nmTimes =  SharedArray{Float64}(length(nmGrid), reps)
# Hold p and q fixed at 600 and vary n and m over a grid
@sync @parallel for j in 1:reps
    for i in 1:length(nmGrid)
        nmTimes[i,j] = @elapsed runSim(lambdas; n=nmGrid[i][1], 
                                       m=nmGrid[i][2], 
                                       p=Int64(mean(pqVals)), 
                                       q=Int64(mean(pqVals)), stepsize=1.0)
    end
end

# Print and write times to CSV
println(reshape(mean(nmTimes, 2), length(nmVals), length(nmVals)))
writecsv("../processed/nm_times.csv",  
         vcat(["n" "m" "mean" transpose(collect(1:reps))], 
              hcat([x[1] for x in nmGrid], [x[2] for x in nmGrid],
                   mean(nmTimes, 2), nmTimes)))
