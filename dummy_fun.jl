using DataFrames


"""
    get_dummy(df, cvar, ctype)

Convert categorical variable to dummy indicators using specified contrast 
type.

# Arguments 

- df = DataFrame of variables
- cvar = a symbol for the categorical variable in df to be converted.
- ctypes = character string indicating the type of contrast to use for cvar

# Value

DataFrame of dummy variables for the specified categorical variable

"""

function get_dummy(df::DataFrames.DataFrame, cvar::Symbol, ctype::String)
    # Obtain the levels to use for the dummy indicators, depending on 
    # contrast type
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
    
    # Iterate through levels to make dummy indicators
    namedict = Dict(zip(vals, 1:length(vals)))
    arr = zeros(length(darr), length(namedict))
    for i=1:length(darr)
        if haskey(namedict, darr[i])
            arr[i, namedict[darr[i]]] = 1
        end
    end
    
    # Some additional modifications for sum contrasts
    if (ctype=="sum")
        arr[darr.==(levels(darr)[end]),:] = -1
    end
    if ctype=="sumnoint"
        arr = arr - (1/length(levels(darr)))
    end
    
    # Convert results to a DataFrame and rename columns 
    newdf = convert(DataFrame, arr)
    names!(newdf, [Symbol("$(cvar)_$k") for k in vals])
    return newdf
end



"""
    contr(df, cvars, ctypes)

Converts categorical variables in a DataFrame to specified contrast types. 
All other variables are left as-is. 

# Arguments 

- df = DataFrame of variables
- cvars = 1d array of symbols corresponding to categorical variable names 
in df to be converted
- ctypes = 1d array of the same length as cvars, indicating the types of 
contrasts to use. Defaults to treatment contrasts ("treat") for all 
variables in cvars. Other options include "sum" for sum contrasts and "noint" 
for treatment contrasts with no intercept.

# Value

DataFrame with same variables as the original DataFrame, with categorical 
variables converted to dummy contrasts. 

"""

function contr(df::DataFrames.DataFrame, cvars::Array{Symbol,1}, 
    ctypes::Array{String,1}=repeat(["treat"], inner=length(cvars)))  
    # Initialize new DataFrame
    newdf = DataFrame()
    
    # Iterate through variables in df
    for var in names(df)
        if !in(var, cvars)
            # Add non-categorical variables to the new DataFrame
            newdf[var] = df[var]
        else
            # Convert categorical variables to specified dummy contrasts
            dummydf = get_dummy(df, var, ctypes[var.==cvars][1])
            for dummyname in names(dummydf)
                newdf[dummyname] = dummydf[dummyname]
            end
        end
    end
    
    return newdf
end
