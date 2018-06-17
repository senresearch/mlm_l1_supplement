using DataFrames

# Convert categorical variables to dummy indicators using specified contrast type.
# df = a data frame.
# cvar = a symbol for the categorical variable in df to be converted.
# ctype = character string indicating the type of contrast.
# trtref = missing. Not used in this version of get_dummy. 
function get_dummy(df::DataFrames.DataFrame, cvar::Symbol, ctype::String, trtref::Missing)
  darr = df[cvar]
  if ctype=="treat"
    vals = levels(darr)[2:end]
  elseif (ctype=="sum")
    vals = levels(darr)[(2:end)-1]
  elseif ctype=="noint" || (ctype=="sumnoint")
    vals = levels(darr)
  else
    error(string("Did not recognize contrast type for ", cvar))
  end
  namedict = Dict(zip(vals, 1:length(vals)))
  arr = zeros(length(darr), length(namedict))
  for i=1:length(darr)
    if haskey(namedict, darr[i])
      arr[i, namedict[darr[i]]] = 1
    end
  end
  if (ctype=="sum")
    arr[darr.==(levels(darr)[end]),:] = -1
  end
  if ctype=="sumnoint"
    arr = arr - (1/length(levels(darr)))
  end
  newdf = convert(DataFrame, arr)
  names!(newdf, [Symbol("$(cvar)_$k") for k in vals])
  return newdf
end

# Convert categorical variables to dummy indicators using specified contrast type.
# Special case of the get_dummy implementation above. For treatment contrasts with a specified
# reference level. 
# df = a data frame.
# cvar = a symbol for the categorical variable in df to be converted.
# ctype = character string indicating the type of contrast. In this case, it must be "treat". 
# trtref = a string specifying the level in cvar to use as the reference. 
function get_dummy(df::DataFrames.DataFrame, cvar::Symbol, ctype::String, trtref::String)
  darr = df[cvar]
  if ctype=="treat"
    vals = levels(darr)[levels(darr) .!= trtref]
  else
    error("Can only specify trtref for treatment contrasts.")
  end
  namedict = Dict(zip(vals, 1:length(vals)))
  arr = zeros(length(darr), length(namedict))
  for i=1:length(darr)
    if haskey(namedict, darr[i])
      arr[i, namedict[darr[i]]] = 1
    end
  end
  newdf = convert(DataFrame, arr)
  names!(newdf, [Symbol("$(cvar)_$k") for k in vals])
  return newdf
end

# Function that sets contrasts by converting categorical values to dummy variables.
# df = a data frame.
# cvars = a symbols vector of categorical variable names in df to be converted.
# ctypes = a vector of the same length as cvars, indicating the types of contrasts.
# Defaults to treatment contrasts ("treat") for all variables in cnames.
# Other options include "sum" for sum contrasts and "noint" for treatment contrasts with no intercept.
# For "treat" ctypes, you can also specify the level to use as the reference treatment. 
# Do this by supplying a vector trtrefs of missing or string values of the same length as ctypes. 
# trtrefs defaults to all missings. 
function contr(df::DataFrames.DataFrame, cvars::Array{Symbol,1}, 
  ctypes::Array{String,1}=repeat(["treat"], inner=length(cvars)), 
  trtrefs::Array=repeat([missing], inner=length(cvars)))
  newdf = DataFrame()
  for var in names(df)
    if !in(var, cvars)
      newdf[var] = df[var]
    else
      dummydf = get_dummy(df, var, ctypes[var.==cvars][1], trtrefs[var.==cvars][1])
      for dummyname in names(dummydf)
        newdf[dummyname] = dummydf[dummyname]
      end
    end
  end
  return newdf
end
