using Distributions

"""
  makeEffect(length, nonzero, dist)

Function to simulate effects with a given proportion of nonzero and the rest drawn from some random distribution

# Arguments 

- length = 
- nonzero = 
- dist = 

# Value

effect vector

"""
function makeEffect(length, nonzero=0.5, dist=Normal(0,2))
  effect = zeros(length)
  effect[sample(1:length, convert(Integer,round(length*nonzero)); replace=false)] = rand(dist, convert(Integer,round(length*nonzero)))
  return effect
end

# Function to set up the Yijs in a simulation.
# nm = dimensions of Y as a tuple
# fixed = fixed effects for Y, should have same dimensions as nm
# rdist, cdist, and edist = distributions or ranges from which the non-fixed effects should be randomly sampled.
# Default to Normal(0,1), the standard normal.

function makeY(n, m, fixed, rdist=Normal(0,1), cdist=Normal(0,1), edist=Normal(0,1))
  return fixed + rand(rdist,n,m) + rand(cdist,n,m) + rand(edist,n,m)
end
