# Load libraries and dependencies
include("../mlm_packages/matrixLMnet/src/matrixLMnet.jl")
using matrixLMnet

using DataFrames

#using MLBase
include("dummyfun.jl")
#@everywhere include("FISTA_backtrack.jl")
#@everywhere include("l1_streamlined.jl")
#include("collins.jl")
include("sim_funs.jl")

# Set number of chemicals and tissues
nchem = 100
ntiss = 10 

# Dimensions of model
n = 108
m = nchem*ntiss
p = 19
q = nchem*ntiss + nchem + ntiss


# Read in X (demographics)
#X = readtable("./processed/woodruff_X.csv", separator = ',', header=true)
# Overparameterized treatment contrasts for the race variable
#Xnoint = convert(Array{Float64}, contr(X, [:race], ["noint"]))
srand(100)
Xnoint = randn(n, p) #size(Xnoint) = (108, 19)


# Contrasts for chemicals and tissues
chem = repmat(eye(nchem), ntiss, 1)
tiss = zeros(nchem*ntiss, ntiss)
for j in 1:ntiss
	tiss[(nchem*(j-1)+1):(nchem*j),j] = 1
end

# Contrasts for chemicals, tissues, and combinations
Znoint = hcat(chem, tiss, eye(nchem*ntiss))



# Generate fixed effects.
srand(40)
dem_eff = makeEffect(p) # Demographic main effects. 1/2 nonzero with SD 2.
chem_eff = repeat(makeEffect(nchem, 1/4), outer=ntiss) # Chemical main effects. 1/4 nonzero with SD 2.
tiss_eff = repeat(makeEffect(ntiss, 1), inner=nchem) # Tissue main effects for each tissue. SD 2
interactions = reshape(makeEffect((p)*(q), 1/8), p, q) # Interaction effects. 1/8 nonzero with SD 2.

# Calculate the fixed effects
fixed_eff = Xnoint*dem_eff .+ transpose(chem_eff .+ tiss_eff) .+ Xnoint*interactions*transpose(Znoint) 
Ysim = makeY(n, m, fixed_eff)
# Standardize Y
Ysim = (Ysim.-mean(Ysim,1))./std(Ysim,1) 


lambdas = reverse(1.3.^(-37:12))

MLM_data = RawData(Response(Ysim), Predictors(Xnoint, Znoint))
results = mlmnet(fista_bt!, MLM_data, lambdas)
flat_coeffs = coef_2d(results)
writecsv("./processed/woodruff_sim_l1_coeffs.csv", flat_coeffs)

writecsv("./processed/woodruff_sim_Y.csv", Ysim)
writecsv("./processed/woodruff_sim_X.csv", Xnoint)

writecsv("./processed/woodruff_sim_lambdas.csv", lambdas)
writecsv("./processed/woodruff_sim_interactions.csv", interactions)

# full_out = run_l1(Xnoint, Y_sim, Znoint, lambdas, "all_simX")

# writecsv("./processed/enviro_Y_sim_simX.csv", Y_sim)
# writecsv("./processed/enviro_Xnoint_simX.csv", Xnoint)

# writecsv("./processed/enviro_lambdas_simX.csv", round(lambdas, 5))
# writecsv("./processed/enviro_sim_inter_simX.csv", inter)