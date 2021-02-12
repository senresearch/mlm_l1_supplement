import LinearAlgebra.I, LinearAlgebra.eigvals
using Random
using Distributed
using DataFrames
import Statistics.mean, Statistics.std
using CSV
using SharedArrays

# L1-penalized matrix linear models
@everywhere using MatrixLMnet


# Read in X (genotype probabilities) with cytoplasm contrast. 
# The first row is a header. 
X = convert(Array{Float64, 2}, CSV.read("../processed/lowry_cyto_genoprobs.csv", 
                                        DataFrame, delim=',', header=true))

# Read in Y (phenotypes). The first row is a header. 
Y = convert(Array{Float64, 2}, CSV.read("../processed/lowry_pheno.csv", 
                                        DataFrame, delim=',', header=true))

# Number of phenotypes 
nPheno = convert(Int64, round(size(Y,2)/2)) 

# Array of 16 lambdas
lambdas = reverse(1.2.^(-9:6))

# Number of random phenotypes to sample
randSamp = [25, 50, 100, 200, 400, 800, 1600] 

# Number of replicates
reps = 100
# Initialize array for storing times
lowryTimes = SharedArrays.SharedArray{Float64}(6, reps, length(randSamp))


for r in 1:length(randSamp)
    # Randomly choose a subset of the phenotypes
    Random.seed!(randSamp[r])
    randPhenos = sort(Random.shuffle(1:nPheno)[1:randSamp[r]])
    
    # Identify corresponding indices to subset from the Y matrix
    idx = Array{Int64}(undef, 2*randSamp[r])
    for i in 1:randSamp[r]
        idx[2*i-1] = randPhenos[i]*2-1
        idx[2*i] = randPhenos[i]*2
    end
    
    # Subset the Y matrix to keep only these phenotypes
    YSub = Y[:,idx]
    
    # Create the Z matrix for this subset of phenotypes
    nPhenoSub = convert(Int64, round(size(YSub,2)/2)) 
    ZSub = kron(Matrix{Float64}(I, nPhenoSub, nPhenoSub), vcat([1 1], [1 -1]))
    
    # Put together RawData object for MLM using subsetted data
    MLMDataSub = RawData(Response(YSub), Predictors(X, ZSub))
    
    # Get times from running ADMM
    @sync @distributed for j in 1:reps
        lowryTimes[6,j,r] = @elapsed mlmnet(admm!, MLMDataSub, lambdas, 
                                            isZInterceptReg=false) #, rho=rhoSub) 
    end
    
    # Get times from running FISTA with backtracking
    @sync @distributed for j in 1:reps
        lowryTimes[5,j,r] = @elapsed mlmnet(fista_bt!, MLMDataSub, lambdas, 
                                            isZInterceptReg=false) 
    end
    
    # Get times from running FISTA with fixed step size
    @sync @distributed for j in 1:reps
        lowryTimes[4,j,r] = @elapsed mlmnet(fista!, MLMDataSub, lambdas, 
                                            isZInterceptReg=false)
    end
    
    # Get times from running ISTA with fixed step size
    @sync @distributed for j in 1:reps
        lowryTimes[3,j,r] = @elapsed mlmnet(ista!, MLMDataSub, lambdas, 
                                            isZInterceptReg=false)
    end
    
    # Get times from running random coordinate descent
    @sync @distributed for j in 1:reps
        lowryTimes[2,j,r] = @elapsed mlmnet(cd_active!, MLMDataSub, lambdas, 
                                            isZInterceptReg=false)
    end
    
    # Get times from running cyclic coordinate descent
    @sync @distributed for j in 1:reps
        lowryTimes[1,j,r] = @elapsed mlmnet(cd_active!, MLMDataSub, lambdas, 
                                            isZInterceptReg=false, isRandom=false)
    end
    
    
    # Print and write times to CSV
    println(mean(lowryTimes[:,:,r], dims=2))
    CSV.write(string("../processed/lowry_times_", randSamp[r], ".csv"),  
            DataFrame(vcat(["method" "mean" transpose(collect(1:reps))], 
                        hcat(["cd_cyclic", "cd_random", "ista", 
                            "fista", "fista_bt", "admm"], 
                            mean(lowryTimes[:,:,r], dims=2), lowryTimes[:,:,r]))), 
            writeheader=false)
    
end