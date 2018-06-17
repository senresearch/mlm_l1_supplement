@everywhere include("../mlm_packages/matrixLMnet/src/matrixLMnet.jl")
@everywhere using matrixLMnet
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
	Y = make_Y(n, m, fixed)
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
reps = 15

println("Starting")

# Dry run 
runSim(lambdas)

pq_vals = collect(200:200:1000)
nm_vals = collect(400:400:2000)

# Hold n and m fixed at 600 and vary p and q over a grid
pq_grid = vec(collect(Base.product(pq_vals, pq_vals)))
pq_times =  SharedArray{Float64}(length(pq_grid), reps)
@sync @parallel for j in 1:reps
	for i in 1:length(pq_grid)
		pq_times[i,j] = @elapsed runSim(lambdas; n=Int64(mean(nm_vals)), m=Int64(mean(nm_vals)), p=pq_grid[i][1], q=pq_grid[i][2])
	end
end

println(reshape(mean(pq_times, 2), length(pq_vals), length(pq_vals)))
writecsv("./processed/pq_times.csv",  
          vcat(["p" "q" "mean" transpose(collect(1:reps))], 
                hcat([x[1] for x in pq_grid], [x[2] for x in pq_grid],
                    	mean(pq_times, 2), pq_times)))


# Hold p and q fixed at 200 and vary n and m over a grid
nm_grid = vec(collect(Base.product(nm_vals, nm_vals)))
nm_times =  SharedArray{Float64}(length(nm_grid), reps)
@sync @parallel for j in 1:reps
	for i in 1:length(nm_grid)
		nm_times[i,j] = @elapsed runSim(lambdas; n=nm_grid[i][1], m=nm_grid[i][2], p=Int64(mean(pq_vals)), q=Int64(mean(pq_vals)))
	end
end

println(reshape(mean(nm_times, 2), length(nm_vals), length(nm_vals)))
writecsv("./processed/nm_times.csv",  
          vcat(["n" "m" "mean" transpose(collect(1:reps))], 
                hcat([x[1] for x in nm_grid], [x[2] for x in nm_grid],
                    	mean(nm_times, 2), nm_times)))
