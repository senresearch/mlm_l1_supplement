# Load libraries and dependencies
@everywhere include("../mlm_packages/matrixLMnet/src/matrixLMnet.jl")
#@everywhere using matrixLMnet

using DataFrames

#using MLBase
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


lambdas = reverse(1.3.^(-37:12))

MLM_data = RawData(Response(Y_sim), Predictors(Xnoint, Znoint))
results = mlmnet(fista_bt!, MLM_data, lambdas)
flat_coeffs = coef_2d(results)
writecsv("./processed/woodruff_sim_l1_coeffs.csv", flat_coeffs)

# full_out = run_l1(Xnoint, Y_sim, Znoint, lambdas, "all_simX")

writecsv("./processed/enviro_Y_sim_simX.csv", Y_sim)
writecsv("./processed/enviro_Xnoint_simX.csv", Xnoint)

writecsv("./processed/enviro_lambdas_simX.csv", round(lambdas, 5))
writecsv("./processed/enviro_sim_inter_simX.csv", inter)