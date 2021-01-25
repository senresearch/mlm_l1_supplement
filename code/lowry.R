library(qtl) # mapping quantitative trait loci

# Load cross
load("../processed/cross_ge.rda")

# Impute missing genotypes
crossGE = fill.geno(crossGE)
# Calculate conditional genotype probabilities
crossGEGenoprobs = calc.genoprob(crossGE, step=0)

###############################################################################

## Process genetic positions

# Load the Arabidopsis annotation file TAIR10_GFF3_genes.gff 
# https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff
genePos = read.delim("../data/TAIR10_GFF3_genes.gff", 
                     header=F, comment.char="#")

# Keep only gene rows and relevant columns
genePos = genePos[genePos$V3 %in% c("gene","transposable_element_gene",
                                    "pseudogene"),c(1,4,5,9)] 
# Give the columns useful names
names(genePos) = c("sequence", "start", "end", "attr") 

# Pull out the gene IDs
genePos$ids = sapply(as.character(genePos$attr), 
                     function(x){substr(x, 4, 12)})

# Adding obsolete names so that they'll merge correctly with the data set. 
# Create a mapping between old (obsolete) and new names
obsolete = data.frame(
  old = c("AT1G06150", "AT1G21580", "AT1G29720", "AT1G55540", "AT1G55790", 
          "AT1G58230", "AT1G65450", "AT1G77210", "AT3G10820", "AT3G18535", 
          "AT3G22430", "AT3G62900", "AT3G63340", "AT4G00060", "AT4G16180", 
          "AT4G19520", "AT4G38760", "AT4G39420", "AT4G39420", "AT5G22720",
          "AT5G40450", "AT5G47400"),
  new = c("AT1G06145", "AT1G21570", "AT1G29724", "AT1G55545", "AT1G55780", 
          "AT1G58227", "AT1G65445", "AT1G77215", "AT3G10830", "AT3G18540", 
          "AT3G22435", "AT3G62895", "AT3G63330", "AT4G00058", "AT4G16170", 
          "AT4G19515", "AT4G38750", "AT4G39440", "AT4G39450", "AT5G22710",
          "AT5G40451", "AT5G47410"))
# Subset the obsolete genes from the list of gene names
obsoletePos = genePos[sapply(obsolete$old, 
                             function(x){which(genePos$ids == x)}),]
# Replace obsolete with new gene names
obsoletePos$ids = obsolete$new
# Concatenate the new names to the entire list of genes
genePos = rbind(genePos, obsoletePos)

# Get the gene IDs from the cross
phenoIDs = data.frame(ids=names(crossGEGenoprobs$pheno[,-(1:4)]))
# Merge them so that the positions are in the same order as the IDs in 
# the cross
genePos = merge(phenoIDs, genePos, by="ids", all.x=TRUE, all.y=FALSE)
# Only keep relevant columns
genePos = genePos[,1:4]

# Write the list of genetic positions to CSV
write.csv(genePos, "../processed/lowry_gene_pos.csv", row.names = FALSE)

###############################################################################

## Create a map between physical and genetic positions

# Read in list of gene positions if it hasn't already been loaded
if (!exists("genePos")) {
  genePos = read.csv("../processed/lowry_gene_pos.csv", header=TRUE)
}

# TKrils_Marker_PhysPos.csv contains physical positions (provided)
TKrils = read.csv("../data/TKrils_Marker_PhysPos.csv", header=TRUE)
# Split by chromosome
TKrilsChr = lapply(1:5, function(i){TKrils[TKrils$Chr==i, ]}) 

# Get the positions of the markers (from the genoprobs) in centimorgans
markerscM = lapply(crossGEGenoprobs$geno, function(x){x$map})

# Get the positions of the phenotype transcripts in base pairs
transbp = lapply(1:5, function(i){
  genePos$start[genePos$sequence==paste("Chr",i,sep="")]})
# Get the positions in centimorgans
transcM = lapply(1:5, function(i){
  approx(TKrilsChr[[i]]$Col.Phys.Pos, 
         TKrilsChr[[i]]$Gen..Pos..LargeMap., 
         transbp[[i]], rule=2)})

# Re-center the centimorgan positions as though the 5 chromosomes were lined 
# up adjacent to each other on a single continuous scale
# Endpoints for the 5 chromosomes
chrcM = c(0, cumsum(sapply(crossGEGenoprobs$geno, function(x){
  max(x$map)})))
# Marker positions, shifted based on chromosome
cumSumMarkerscM = do.call(c, lapply(1:5, function(i){
  markerscM[[i]] + chrcM[i]}))
# Transcript positions, shifted based on chromosme
cumSumTranscM = do.call(c, lapply(1:5, function(i){
  transcM[[i]]$y + chrcM[i]}))
# Double each transcript position (for dry/wet)
cumSumTranscM2 = rep(cumSumTranscM, each=2)

# Save cM positions
save(chrcM, cumSumMarkerscM, cumSumTranscM, cumSumTranscM2, 
     file="../processed/lowry_plot_data.rData")

###############################################################################

## Plot results

# Load physical distances if they haven't already been loaded
if (!exists("chrcM") || !exists("cumSumMarkerscM") 
    || !exists("cumSumTranscM2")) {
  load("../processed/lowry_plot_data.rData")
}

# Read in L1 coefficient estimates
lowryL1Coeffs = read.csv("../processed/lowry_l1_coeffs.csv", header=FALSE)
# Pull out coefficients corresponding to lambda of interest
lambda = 4 # 4th lambda has value 1.728
coeffs = matrix(lowryL1Coeffs[,lambda], 452, 51324)

# Remove the X intercept and the cyto contrast
coeffs = coeffs[-(1:2), ]
# Find the indices of the non-zero interactions
idx = which(coeffs != 0, arr.ind=TRUE)

# Indices of the main effects
idxMains = idx[idx[,2] %% 2 != 0,] 
# Indices of the interactions 
idxInter = idx[idx[,2] %% 2 == 0,] 

# Chromosome names used for labeling
chrNames = paste(paste("Chr", 1:5, sep=""), 
                 collapse = paste(replicate(14, " "), collapse = ""))

# Map the coefficients to their physical positions
# Main effects
mains = cbind(cumSumMarkerscM[idxMains[,1]], 
              cumSumTranscM2[idxMains[,2]])
# Interactions
inter = cbind(cumSumMarkerscM[idxInter[,1]], 
              cumSumTranscM2[idxInter[,2]])


# Color
png("../pictures/lowry_gene_vs_qtl_pos_lambda_1.7.png", 
    width=380, height=380)
par(mar=c(4.1,4.1,1.1,1.1))

# Main effects
plot(mains[,1], mains[,2], cex=0.8, col="royalblue4", 
     xlab="QTL Position (cM)", ylab="Gene Position (cM)") 

# Interactions
points(inter[,1], inter[,2], pch=15, col="firebrick3")

# Reference lines to delineate chromosomes
invisible(sapply(chrcM[2:5], function(x){abline(v=x)}))
invisible(sapply(chrcM[2:5], function(x){abline(h=x)}))

# Label chromosomes
mtext(side=1, text=chrNames, line=2)
mtext(side=2, text=chrNames, line=2) 

dev.off()


# Black and white
png("../pictures/lowry_gene_vs_qtl_pos_lambda_1.7_bw.png", 
    width=380, height=380)
par(mar=c(4.1,4.1,1.1,1.1))

# Main effects
plot(mains[,1], mains[,2], cex=0.8, 
     xlab="QTL Position (cM)", ylab="Gene Position (cM)") 

# Interactions
points(inter[,1], inter[,2], pch=15)

# Reference lines to delineate chromosomes
invisible(sapply(chrcM[2:5], function(x){abline(v=x)}))
invisible(sapply(chrcM[2:5], function(x){abline(h=x)}))

# Label chromosomes
mtext(side=1, text=chrNames, line=2)
mtext(side=2, text=chrNames, line=2) 

dev.off()
