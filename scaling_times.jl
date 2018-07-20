# L1-penalized matrix linear models
@everywhere include("../mlm_packages/matrixLMnet/src/matrixLMnet.jl")
@everywhere using matrixLMnet

# Functions for simulating data
@everywhere include("sim_funs.jl")


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
        MLMdata = simRawData(n, m, p, q, seed)
        
        # Run L1-penalized matrix linear model
        return mlmnet(fun, MLMdata, lambdas; funArgs...)
    end
end


# Array of 50 lambdas
lambdas = reverse(1.2.^(-32:17))
# Number of replicates
reps = 10

# Range of p and q values to try
pq_vals = collect(200:200:1000)
# Range of n and m values to try 
nm_vals = collect(400:400:2000)


# Generate grid of p and q values
pq_grid = vec(collect(Base.product(pq_vals, pq_vals)))
# Initialize array for storing times
pq_times =  SharedArray{Float64}(length(pq_grid), reps)
# Hold n and m fixed at 1200 and vary p and q over a grid
@sync @parallel for j in 1:reps
    for i in 1:length(pq_grid)
        pq_times[i,j] = @elapsed runSim(lambdas; n=Int64(mean(nm_vals)), 
                                        m=Int64(mean(nm_vals)), 
                                        p=pq_grid[i][1], 
                                        q=pq_grid[i][2], stepsize=1.0)
    end
end

# Print and write times to CSV
println(reshape(mean(pq_times, 2), length(pq_vals), length(pq_vals)))
writecsv("./processed/pq_times.csv",  
          vcat(["p" "q" "mean" transpose(collect(1:reps))], 
               hcat([x[1] for x in pq_grid], [x[2] for x in pq_grid],
                    mean(pq_times, 2), pq_times)))


# Generate grid of n and m values
nm_grid = vec(collect(Base.product(nm_vals, nm_vals)))
# Initialize array for storing times
nm_times =  SharedArray{Float64}(length(nm_grid), reps)
# Hold p and q fixed at 600 and vary n and m over a grid
@sync @parallel for j in 1:reps
    for i in 1:length(nm_grid)
        nm_times[i,j] = @elapsed runSim(lambdas; n=nm_grid[i][1], 
                                        m=nm_grid[i][2], 
                                        p=Int64(mean(pq_vals)), 
                                        q=Int64(mean(pq_vals)), stepsize=1.0)
    end
end

# Print and write times to CSV
println(reshape(mean(nm_times, 2), length(nm_vals), length(nm_vals)))
writecsv("./processed/nm_times.csv",  
         vcat(["n" "m" "mean" transpose(collect(1:reps))], 
              hcat([x[1] for x in nm_grid], [x[2] for x in nm_grid],
                   mean(nm_times, 2), nm_times)))
