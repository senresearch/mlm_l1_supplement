"""
    make_effect(n, propNonzero, dist)

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

function make_effect(n::Int64, propNonzero::Float64=0.5, 
                     dist::Distribution=Normal(0,2))
    # Initialize vector for storing effects 
    effect = zeros(n)
    
    # Randomly sample indices of nonzero effects
    idx = sample(1:n, convert(Integer, round(n*propNonzero)); replace=false)
    
    # Simulate and assign nonzero effects  
    effect[idx] = rand(dist, convert(Integer, round(n*propNonzero)))
    
    return effect
end
