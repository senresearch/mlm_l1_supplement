library(data.table) # quickly read in tables

# Read in keys (with mutant and spatial information) for each of the 6 plates. 
KEIOKeys = lapply(1:6, function(i){
  read.csv(paste0("../data/raw_KEIO_data/KEIO", i, "_KEY.csv"), sep="\t")
})

# Read in the Kritikos condition information. 
kritCond = lapply(1:6, function(i){
  read.csv(paste0("../data/raw_KEIO_data/krit_condition_files/p", i, 
                  "_4120krit.csv"), sep="\t")
})

###############################################################################

# Function to read in the iris/dat files and pull out the opacity column.
# Assumes that rows and columns in the data file are in the same order as in 
# the mutant key. 
# datFile = character string naming the .iris or .dat file to read
# mutantKey = mutant key data frame to use
# directory = directory name holding the data files
# ... additional arguments to be passed into read.csv
read_dat = function(datFile, mutantKey, directory, ...) {
  dat_temp = read.csv(paste(directory, datFile, sep=""), sep="\t", ...)
  if(all(mutantKey$row == dat_temp$row) && 
     all(mutantKey$column == dat_temp$column)) {
    return(dat_temp$opacity)
  } else {
    stop("Check your rows and columns.")
  }
}

# Get the Y response matrix of colony opacities for the Kritikos conditions
kritDat = lapply(1:6, function(i){
  t(sapply(kritCond[[i]]$KRIT.FILE, read_dat, KEIOKeys[[i]], 
           "../data/raw_KEIO_data/krit_dat/"))
})

# Plate 5: Drop rows with conditions "novobiocin" and concentration "null"
kritDat[[5]] = kritDat[[5]][!(kritCond[[5]]$Condition=="novobiocin" & 
                                  kritCond[[5]]$Concentration=="null"),]
kritCond[[5]] = kritCond[[5]][!(kritCond[[5]]$Condition=="novobiocin" & 
                                    kritCond[[5]]$Concentration=="null"), ]

###############################################################################

# Function to clean the concentrations. # Replaces the concentration "20ug/mL" 
# with "20 ug/mL" for the condition "vancomycin". Replaces 
#   - multi-spaces with a single space
#   - ";" with " +"
#   - " sec" with "sec"
# concCol = column of concentrations
# condCol = column of conditions
clean_concentrations = function(concCol, condCol) {
  concCol[condCol=="vancomycin" & concCol=="20ug/mL"] = "20 ug/mL"
  return(gsub("\\s+", " ", 
              gsub(";", " +", 
                   gsub(" sec", "sec", concCol))))
}

# Iterate through Kritikos conditions for each plate 
for (i in 1:6) {
  # Clean concentrations
  kritCond[[i]]$Concentration = 
    clean_concentrations(kritCond[[i]]$Concentration, 
                         kritCond[[i]]$Condition)
  # Create new variable by concatenating conditions and concentrations
  kritCond[[i]]$Cond_Conc = paste(kritCond[[i]]$Concentration, 
                                   kritCond[[i]]$Condition)
}

###############################################################################

# Pull out names of condition-concentrations and mutants
kritCondConcNames = lapply(kritCond, function(x){x$Cond_Conc})
kritMutNames = lapply(KEIOKeys, function(x){as.character(x$name)})

# Total number of condition-concentrations/unique condition-concentrations
length(do.call(c, kritCondConcNames))
length(unique(do.call(c, kritCondConcNames)))

# Total number of mutants/unique mutants
length(do.call(c, kritMutNames))
length(unique(do.call(c, kritMutNames)))

# Number of condition-concentrations/unique condition-concentrations in each 
# plate
sapply(kritCondConcNames, length)
sapply(kritCondConcNames, function(x){length(unique(x))})

# Number of mutants/unique mutants in each plate
sapply(kritMutNames, length)
sapply(kritMutNames, function(x){length(unique(x))})

###############################################################################

# Write processed data to CSV
for (i in 1:6) {
  # X covariate matrix coded with condition-concentrations
  write.csv(kritCond[[i]], 
            file=paste0("../processed/processed_KEIO_data/p", i, 
                        "_krit_cond.csv"), row.names=FALSE)
  # Y response matrix of colony opacities
  write.csv(kritDat[[i]], 
            file=paste0("../processed/processed_KEIO_data/p", i, 
                        "_krit_dat.csv"), row.names=FALSE)
  
  # Condition-concentration names
  write.table(kritCondConcNames[[i]], 
              file=paste0("../processed/processed_KEIO_data/p", i, 
                          "_krit_cond_conc_names.csv"), 
              sep=",", row.names=FALSE, col.names=FALSE)
  # Mutation names
  write.table(kritMutNames[[i]], 
              file=paste0("../processed/processed_KEIO_data/p", i, 
                          "_krit_mut_names.csv"), 
              sep=",", row.names=FALSE, col.names=FALSE)
}
