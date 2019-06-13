using DataFrames
using LinearAlgebra
using Distributions
using Random
using CSV

# L1-penalized matrix linear models
include("../../mlm_packages/matrixLMnet/src/matrixLMnet.jl")
using Main.matrixLMnet


"""
    sim_effect(n, propNonzero; eDist)

Simulate effects with a given proportion of nonzero effects drawn from some 
distribution. The remaining effects will be set to zero.

# Arguments 

- n = length of 1d effect array. 
- propNonzero = proportion of nonzero effects. Defaults to `1/2`.

# Keyword arguments

- eDist = distribution from which the nonzero effects should be simulated. 
  Defaults to `Normal(0,2)`. 

# Value

1d array of floats 

"""
function sim_effect(n::Int64, propNonzero::Float64=1/2, 
                    eDist::Distribution=Normal(0,2))
    
    # Initialize vector for storing effects 
    effect = zeros(n)
    
    # Randomly sample indices of nonzero effects
    idx = sample(1:n, convert(Integer, round(n*propNonzero)); replace=false)
    
    # Simulate and assign nonzero effects  
    effect[idx] = rand(eDist, convert(Integer, round(n*propNonzero)))
    
    return effect
end


# Number of chemicals 
nChem = 100
# Number of tissues
nTiss = 10 

# Dimensions of data
n = 108
m = nChem*nTiss
p = 19
q = nChem*nTiss + nChem + nTiss


# Simulate X (demographics)
Random.seed!(100)
X = randn(n, p) 

# Create contrasts for chemicals and tissues
chem = repeat(Matrix{Float64}(LinearAlgebra.I, nChem, nChem), nTiss, 1)
tiss = zeros(nChem*nTiss, nTiss)
for j in 1:nTiss
    tiss[(nChem*(j-1)+1):(nChem*j),j] .= 1
end
# Create Z matrix
Z = hcat(chem, tiss, Matrix{Float64}(LinearAlgebra.I, nChem*nTiss, nChem*nTiss))


Random.seed!(40)
# Simulate demographic main effects, 1/2 nonzero from Normal(0,2). 
demEff = sim_effect(p) 
# Simulate chemical main effects, 1/4 nonzero from Normal(0,2). 
chemEff = repeat(sim_effect(nChem, 1/4), outer=nTiss) 
# Simulate tissue main effects, from Normal(0,2). 
tissEff = repeat(sim_effect(nTiss, 1.0), inner=nChem) 
# Simulate interaction effects, 1/8 nonzero from Normal(0,2). 
interactions = reshape(sim_effect(p*q, 1/8), p, q) 

# Generate the fixed effects
fixedEff = X*demEff .+ transpose(chemEff .+ tissEff) .+ 
             X*interactions*transpose(Z) 
# Simulate Y using fixed effects 
YSim = fixedEff + rand(Normal(0,3), n, m) 
# Standardize Y
YSim = (YSim.-mean(YSim, dims=1))./std(YSim, dims=1) 


# Put together RawData object for MLM 
MLMSimData = RawData(Response(YSim), Predictors(X, Z))

# Array of 50 lambdas
lambdas = reverse(1.3.^(-37:12))


# Run L1-penalized matrix linear model
results = mlmnet(fista_bt!, MLMSimData, lambdas)

# Flatten coefficients and write results to CSV
flat_coeffs = coef_2d(results)
CSV.write("../processed/woodruff_sim_l1_coeffs.csv", 
          DataFrame(flat_coeffs), writeheader=false)

# Write simualted X, Y, and interactions to CSV
CSV.write("../processed/woodruff_sim_Y.csv", 
          DataFrame(YSim), writeheader=false)
CSV.write("../processed/woodruff_sim_X.csv", 
          DataFrame(X), writeheader=false)
CSV.write("../processed/woodruff_sim_interactions.csv", 
          DataFrame(interactions), writeheader=false)
