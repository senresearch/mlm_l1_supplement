using Distributions


"""
    make_effect(n, prop_nonzero, dist)

Simulate effects with a given proportion of nonzero effects drawn from some 
distribution. The remaining effects will be set to zero.

# Arguments 

- n = length of 1d effect array. 
- prop_nonzero = proportion of nonzero effects. Defaults to 0.5.
- dist = distribution from which the nonzero effects should be simulated. 
Defaults to Normal(0,2). 

# Value

1d array of effects 

"""

function make_effect(n::Int64, prop_nonzero::Float64=0.5, 
                     dist::Distribution=Normal(0,2))
    # Initialize vector for storing effects 
    effect = zeros(n)
    
    # Randomly sample indices of nonzero effects
    idx = sample(1:n, convert(Integer,round(n*prop_nonzero)); replace=false)
    
    # Simulate and assign nonzero effects  
    effect[idx] = rand(dist, convert(Integer,round(n*prop_nonzero)))
    
    return effect
end


"""
    make_Y(n, m, fixed, edist)

Simulate a 2d response array for simulations, with user-specified fixed 
effects. 

# Arguments 

- n = number of rows
- m = number of columns
- fixed = 2d array of fixed effects (should be n by m)
- edist = distribution from which the non-fixed effects should be randomly 
sampled. Defaults to Normal(0,3). 

# Value

2d response array

"""

function make_Y(n::Int64, m::Int64, fixed::Array{Float64,2}, 
                edist::Distribution=Normal(0,3))
    return fixed + rand(edist, n, m)
end
