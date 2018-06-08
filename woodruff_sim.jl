# Load libraries and dependencies
@everywhere include("../mlm_packages/matrixLMnet/src/matrixLMnet.jl")
#@everywhere using matrixLMnet

using DataFrames
using Distributions
using MLBase
include("dummyfun.jl")
#@everywhere include("FISTA_backtrack.jl")
#@everywhere include("l1_streamlined.jl")
#include("collins.jl")
include("sim_funs.jl")


# Read in X (demographics)
X = readtable("./processed/enviro_X.csv", separator = ',', header=true)
# Overparameterized treatment contrasts for the race variable
Xnoint = convert(Array{Float64}, contr(X, [:race], ["noint"]))
Xnoint = randn(size(Xnoint))

# Set number of chemicals and tissues
nchem = 100
ntiss = 10 

# Contrasts for chemicals and tissues
chem = repmat(eye(nchem), ntiss, 1)
tiss = zeros(nchem*ntiss, ntiss)
for j in 1:ntiss
	tiss[(nchem*(j-1)+1):(nchem*j),j] = 1
end

# Contrasts for chemicals, tissues, and combinations
Znoint = hcat(chem, tiss, eye(nchem*ntiss))

# Dimensions of model
n = size(Xnoint,1)
m = size(Znoint,1)
p = size(Xnoint,2)
q = size(Znoint,2)

# Generate fixed effects.
srand(40)
d = makeEffect(p) # Demographic main effects. 1/2 nonzero with SD 2.
c = repeat(makeEffect(nchem, 1/4), outer=ntiss) # Chemical main effects. 1/4 nonzero with SD 2.
t = repeat(makeEffect(ntiss, 1), inner=nchem) # Tissue main effects for each tissue. SD 2
inter = reshape(makeEffect((p)*(q), 1/8), p, q) # Interaction effects. 1/8 nonzero with SD 2.

# Calculate the fixed effects
fixed = Xnoint*d .+ transpose(c .+ t) .+ Xnoint*inter*transpose(Znoint) 
Y_sim = makeY(n,m, fixed)
# Standardize Y
Y_sim = (Y_sim.-mean(Y_sim,1))./std(Y_sim,1) 

# Modified l1_pathwise function that writes the results to CSV
@everywhere function l1_pathwise(fun::Function, X::Array{Float64,2}, Y::Array{Float64,2}, Z::Array{Float64,2}, lambdas::Array{Float64,1}, fold,  
				   regXidx::Array{Int64,1}, regZidx::Array{Int64,1}, reg::BitArray{2}, norms; stepsize=1.0, funArgs...)
  # Check that the lambdas are in descending order. 
  if any(lambdas .!= sort(lambdas, rev=true))
    println("Sorting lambdas into descending order.")
    lambdas = sort(lambdas, rev=true)
  end

  # Pre-allocate array for coefficients
  coeffs = Array(Float64, length(lambdas), size(X,2), size(Z,2)) 

  # Start with coefficients initalized at zero for the largest lambda value
  startB = zeros(size(X,2), size(Z,2))

  # Iterate through the path of lambdas
  for i=1:length(lambdas) 
  	# Get L1-penalty estimates by updating the coefficients from previous iteration in place
  	fun(X, Y, Z, lambdas[i], startB, regXidx, regZidx, reg, norms; stepsize=stepsize, funArgs...)
    # Assign a slice of coeffs to the current coefficient estimates
  	coeffs[i,:,:] = startB 
	writecsv(string("./processed/enviro_sim_lambda",round(lambdas[i], 5),"_fold", fold,".csv"), reshape(coeffs[i,:,:], size(X,2), size(Z,2)))
  end

  return coeffs
end

"""
  run_l1(fun, X, Y, Z, lambdas, fold; 
        filename, isXIntercept, isZIntercept, isXReg, isZReg, isXInterceptReg, isZInterceptReg, isStandardize, stepsize, funArgs...)

Standardizes X and Z, performs the supplied method on a descending list of lambdas using ``warm starts'', 
and backtransforms resulting coefficients. 

# Arguments

- fun = function that applies an L1-penalty estimate method
- X = 2d array of floats consisting of the row covariates, with all categorical variables coded in appropriate contrasts
- Y = 2d array of floats consisting of the multivariate response
- Z = 2d array of floats consisting of the column covariates, with all categorical variables coded in appropriate contrasts
- lambdas = 1d array of floats consisting of lambda penalties in descending order
    If they are not in descending order, they will be sorted

# Keyword arguments

- filename = optional string representing a file path to write out the coefficients 
    Defaults to `nothing`, in which case nothing will be written out
    Should be a `.csv` file (currently)
- isXIntercept = boolean flag indicating whether or not to include an X intercept (row main effects)
    Defaults to true
- isZIntercept = boolean flag indicating whether or not to include a Z intercept (column main effects)
    Defaults to true
- isXReg = 2d array of bit flags indicating whether or not to regularize each of the X (row main) effects
    Defaults to 2d array of trues with length equal to the number of X effects
- isZReg = 2d array of bit flags indicating whether or not to regularize each of the Z (column main) effects
    Defaults to 2d array of trues with length equal to the number of Z effects
- isXInterceptReg = boolean flag indicating whether or not to regularize the X intercept
    Defaults to false 
- isZInterceptReg = boolean flag indicating whether or not to regularize the Z intercept
    Defaults to false
- isStandardize = boolean flag indicating if X and Z should be standardized (mean 0, standard deviation 1) 
    Defaults to true.
- stepsize = float; step size for updates
    Defaults to 1.0
- funArgs = variable keyword arguments to be passed into fun

# Value

A 2d array consisting of the coefficient estimates. Each column contains the vectorized estimates for a given lambda.

"""

function run_l1(X::Array{Float64,2}, Y::Array{Float64,2}, Z::Array{Float64,2}, lambdas::Array{Float64,1}, fold; 
  filename=nothing, isXIntercept::Bool=true, isZIntercept::Bool=true, 
	isXReg::BitArray{1}=trues(size(X,2)), isZReg::BitArray{1}=trues(size(Z,2)), 
  isXInterceptReg::Bool=false, isZInterceptReg::Bool=false, 
  isStandardize::Bool=true, stepsize=1.0, funArgs...)

  fun = fista_bt!
  # Add on X and Z intercepts, if necessary. 
  if isXIntercept==true # Include X intercept
    X = addIntercept(X)
    isXReg = vcat(isXInterceptReg, isXReg)
  end
  if isZIntercept==true # Include Z intercept
    Z = addIntercept(Z)
    isZReg = vcat(isZInterceptReg, isZReg)
  end
  reg = isXReg.*transpose(isZReg) # Matrix to keep track of which coefficients to regularize.
  regXidx = find(isXReg) # Indices corresponding to regularized X covariates. 
  regZidx = find(isZReg) # Indices corresponding to regularized Z covariates. 

  # If intercepts are included and predictors will be standardized, save copies of X and Z for back-transforming later
  if (isStandardize==true) && (isXIntercept==true) && (isZIntercept==true) 
    Xold = copy(X)
    Zold = copy(Z)
  end
  # Standardize predictors if necessary.
  if isStandardize == true 
    meansX, normsX, = standardize!(X, isXIntercept) 
    meansZ, normsZ, = standardize!(Z, isZIntercept)
    norms = nothing # If X and Z are standardized, set the norm to nothing
  else # If not standardizing, calculate the norm matrix
    norms = transpose(sum(X.^2,1)).*sum(Z.^2,1) # 2d array of norms corresponding to each coefficient
  end

  # If chosen method is ista/fista with fixed step size, compute the step size. 
  if (((string(fun) == "ista!") | (string(fun) == "fista!")) & (stepsize == 1.0))
  	# Step size calculate as reciprocal of maximum eigenvalue of kron(Z, X)
  	XTX = transpose(X)*X
  	ZTZ = transpose(Z)*Z
  	stepsize = 1/(maximum(abs(eigvals(XTX))) * maximum(abs(eigvals(ZTZ))))
  	if isStandardize==true
    	# There are some issues if doing this for standardized X and Z (complex eigenvalues)
    	# Hack is to add diagonal matrix where the diagonal is random normal noise
    	stepsize = 1/max(eigmax(transpose(X)*X + diagm(1.0 + randn(size(X,2))/1000)) * eigmax(transpose(Z)*Z + diagm(1.0 + randn(size(Z,2))/1000)), 
    	eigmin(transpose(X)*X + diagm(1.0 + randn(size(X,2))/10000)) * eigmin(transpose(Z)*Z + diagm(1.0 + randn(size(Z,2))/1000)))
  	else 
    	stepsize = 1/max(eigmax(transpose(X)*X) * eigmax(transpose(Z)*Z), 
    	eigmin(transpose(X)*X) * eigmin(transpose(Z)*Z))
  	end
  end

  # Run the specified L1-penalty method on the supplied inputs. 
  coeffs = l1_pathwise(fun, X, Y, Z, lambdas, fold, regXidx, regZidx, reg, norms; stepsize=stepsize/100, funArgs...)
  
  # Back-transform coefficient estimates, if necessary. 
  if isStandardize == true  && (isXIntercept==true) && (isZIntercept==true)  # Case if including both X and Z intercepts. 
    backtransform!(coeffs, meansX, meansZ, normsX, normsZ, Y, Xold, Zold)
  elseif isStandardize == true # Otherwise
    backtransform!(coeffs, isXIntercept, isZIntercept, meansX, meansZ, normsX, normsZ)
  end

  # If a file path is given, flatten the coefficients into column vectors and write them to the path.
  if (filename != nothing) 
  	flat_coeffs = Array(Float64, size(X,2)*size(Z,2), length(lambdas))
  	for i=1:length(lambdas)
  		flat_coeffs[:,i] = vec(coeffs[i,:,:])
  	end
  	writecsv(filename, flat_coeffs)
  end

  return coeffs
end


lambdas = reverse(1.3.^(-37:12))
full_out = run_l1(Xnoint, Y_sim, Znoint, lambdas, "all_simX")

writecsv("./processed/enviro_Y_sim_simX.csv", Y_sim)
writecsv("./processed/enviro_Xnoint_simX.csv", Xnoint)

# writecsv("./processed/enviro_lambdas_simX.csv", round(lambdas, 5))
writecsv("./processed/enviro_sim_inter_simX.csv", inter)