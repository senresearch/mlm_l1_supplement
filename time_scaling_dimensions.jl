@everywhere include("../mlm_packages/matrixLMnet/src/matrixLMnet.jl")
@everywhere using matrixLMnet

@everywhere using Distributions

@everywhere include("sim_funs.jl")



@everywhere function simRawData(n, m, p, q, seed=10)
	srand(seed)
	
	# Simulate effects 
	d = makeEffect(p) # Row main effects. 1/2 nonzero with SD 2.
	c = makeEffect(q) # Column main effects. 1/2 nonzero with SD 2.
	inter = reshape(makeEffect((p)*(q), 1/8), p, q) # Interaction effects. 1/8 nonzero with SD 2.

	# Simulate X and Z
	X = repeat(eye(p), outer=[Int(ceil(n/p)),1])[1:n,:]
	Z = repeat(eye(q), outer=[Int(ceil(m/q)),1])[1:m,:]
	
	# Get the fixed effects and generate Y
	# Calculate the fixed effects
	fixed = X*d .+ transpose(Z*c) .+ A_mul_Bt(X*inter, Z)
	Y = makeY(n, m, fixed)
	# Standardize Y
	Y = (Y.-mean(Y,1))./std(Y,1) 
	
	# Put together RawData object for MLM
	return RawData(Response(Y), Predictors(X, Z))
end


@everywhere function runSim(lambdas, fun=fista_bt!; n=600, m=600, p=200, q=200, seed=10)
	MLM_data = simRawData(n, m, p, q, seed)
	
	return mlmnet(fun, MLM_data, lambdas)
end

lambdas = reverse(1.2.^(-32:17))
reps = 10

println("Starting")

# Dry run 
runSim(lambdas)

# Hold n and m fixed at 600 and vary p and q over a grid
pq_grid = vec(collect(Base.product(collect(100:100:600), collect(100:100:600))))
pq_times =  SharedArray{Float64}(length(pq_grid), reps)
@sync @parallel for j in 1:reps
	for i in 1:length(pq_grid)
		pq_times[i,j] = @elapsed runSim(lambdas; n=600, m=600, p=pq_grid[i][1], q=pq_grid[i][2])
	end
end

println(reshape(mean(pq_times, 2), 6, 6))
writecsv("./processed/pq_times.csv", hcat(mean(pq_times, 2), pq_times))


# Hold p and q fixed at 200 and vary n and m over a grid
nm_grid = vec(collect(Base.product(collect(200:200:1200), collect(200:200:1200))))
nm_times =  SharedArray{Float64}(length(nm_grid), reps)
@sync @parallel for j in 1:reps
	for i in 1:length(nm_grid)
		nm_times[i,j] = @elapsed runSim(lambdas; n=nm_grid[i][1], m=nm_grid[i][2], p=200, q=200)
	end
end

println(reshape(mean(nm_times, 2), 6, 6))
writecsv("./processed/nm_times.csv", hcat(mean(nm_times, 2), nm_times))
