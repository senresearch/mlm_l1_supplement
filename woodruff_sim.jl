# L1-penalized matrix linear models
include("../mlm_packages/matrixLMnet/src/matrixLMnet.jl")
using matrixLMnet

# DataFrames 
using DataFrames

# Functions for simulating data
include("sim_funs.jl")


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
srand(100)
X = randn(n, p) 

# Create contrasts for chemicals and tissues
chem = repmat(eye(nChem), nTiss, 1)
tiss = zeros(nChem*nTiss, nTiss)
for j in 1:nTiss
	tiss[(nChem*(j-1)+1):(nChem*j),j] = 1
end
# Create Z matrix
Z = hcat(chem, tiss, eye(nChem*nTiss))


srand(40)
# Simulate demographic main effects, 1/2 nonzero from Normal(0,2). 
demEff = make_effect(p) 
# Simulate chemical main effects, 1/4 nonzero from Normal(0,2). 
chemEff = repeat(make_effect(nChem, 1/4), outer=nTiss) 
# Simulate tissue main effects, from Normal(0,2). 
tissEff = repeat(make_effect(nTiss, 1.0), inner=nChem) 
# Simulate interaction effects, 1/8 nonzero from Normal(0,2). 
interactions = reshape(make_effect((p)*(q), 1/8), p, q) 

# Generate the fixed effects
fixedEff = X*demEff .+ transpose(chemEff .+ tissEff) .+ 
            X*interactions*transpose(Z) 
# Simulate Y using fixed effects 
YSim = make_Y(n, m, fixedEff)
# Standardize Y
YSim = (YSim.-mean(YSim,1))./std(YSim,1) 


# Put together RawData object for MLM 
MLMData = RawData(Response(YSim), Predictors(X, Z))

# Array of 50 lambdas
lambdas = reverse(1.3.^(-37:12))


# Run L1-penalized matrix linear model
results = mlmnet(fista_bt!, MLMData, lambdas)

# Flatten coefficients and write results to CSV
flat_coeffs = coef_2d(results)
writecsv("./processed/woodruff_sim_l1_coeffs.csv", flat_coeffs)

# Write simualted X, Y, and interactions to CSV
writecsv("./processed/woodruff_sim_Y.csv", YSim)
writecsv("./processed/woodruff_sim_X.csv", X)
writecsv("./processed/woodruff_sim_interactions.csv", interactions)
