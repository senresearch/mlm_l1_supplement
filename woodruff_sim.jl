# L1-penalized matrix linear models
include("../mlm_packages/matrixLMnet/src/matrixLMnet.jl")
using matrixLMnet

# DataFrames 
using DataFrames

# Functions for simulating data
include("sim_funs.jl")


# Number of chemicals 
nchem = 100
# Number of tissues
ntiss = 10 

# Dimensions of data
n = 108
m = nchem*ntiss
p = 19
q = nchem*ntiss + nchem + ntiss


# Simulate X (demographics)
srand(100)
X = randn(n, p) 

# Create contrasts for chemicals and tissues
chem = repmat(eye(nchem), ntiss, 1)
tiss = zeros(nchem*ntiss, ntiss)
for j in 1:ntiss
	tiss[(nchem*(j-1)+1):(nchem*j),j] = 1
end
# Create Z matrix
Z = hcat(chem, tiss, eye(nchem*ntiss))


srand(40)
# Simulate demographic main effects, 1/2 nonzero from Normal(0,2). 
dem_eff = make_effect(p) 
# Simulate chemical main effects, 1/4 nonzero from Normal(0,2). 
chem_eff = repeat(make_effect(nchem, 1/4), outer=ntiss) 
# Simulate tissue main effects, from Normal(0,2). 
tiss_eff = repeat(make_effect(ntiss, 1.0), inner=nchem) 
# Simulate interaction effects, 1/8 nonzero from Normal(0,2). 
interactions = reshape(make_effect((p)*(q), 1/8), p, q) 

# Generate the fixed effects
fixed_eff = X*dem_eff .+ transpose(chem_eff .+ tiss_eff) .+ 
            X*interactions*transpose(Z) 
# Simulate Y using fixed effects 
Ysim = make_Y(n, m, fixed_eff)
# Standardize Y
Ysim = (Ysim.-mean(Ysim,1))./std(Ysim,1) 


# Put together RawData object for MLM 
MLMdata = RawData(Response(Ysim), Predictors(X, Z))

# Array of 50 lambdas
lambdas = reverse(1.3.^(-37:12))


# Run L1-penalized matrix linear model
results = mlmnet(fista_bt!, MLMdata, lambdas)

# Flatten coefficients and write results to CSV
flat_coeffs = coef_2d(results)
writecsv("./processed/woodruff_sim_l1_coeffs.csv", flat_coeffs)

# Write simualted X, Y, and interactions to CSV
writecsv("./processed/woodruff_sim_Y.csv", Ysim)
writecsv("./processed/woodruff_sim_X.csv", X)
writecsv("./processed/woodruff_sim_interactions.csv", interactions)
