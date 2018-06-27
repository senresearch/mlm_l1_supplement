library(qtl) # mapping quantitative trait loci

# Load cross
load("./processed/cross.ge.rda")

# Impute missing genotypes
cross.ge = fill.geno(cross.ge)
# Calculate conditional genotype probabilities
cross.ge.genoprobs = calc.genoprob(cross.ge, step=0)

###############################################################################

## Process genetic positions

# Load the Arabidopsis annotation file TAIR10_GFF3_genes.gff 
# https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff
gene.pos = read.delim("./processed/TAIR10_GFF3_genes.gff", 
                      header=F, comment.char="#")

# Keep only gene rows and relevant columns
gene.pos = gene.pos[gene.pos$V3 %in% c("gene","transposable_element_gene",
                                       "pseudogene"),c(1,4,5,9)] 
# Give the columns useful names
names(gene.pos) = c("sequence", "start", "end", "attr") 

# Pull out the gene IDs
gene.pos$ids = sapply(as.character(gene.pos$attr), 
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
obsolete_pos = gene.pos[sapply(obsolete$old, 
                               function(x){which(gene.pos$ids == x)}),]
# Replace obsolete with new gene names
obsolete_pos$ids = obsolete$new
# Concatenate the new names to the entire list of genes
gene.pos = rbind(gene.pos, obsolete_pos)

# Get the gene IDs from the cross
pheno_ids = data.frame(ids=names(cross.ge.genoprobs$pheno[,-(1:4)]))
# Merge them so that the positions are in the same order as the IDs in 
# the cross
gene.pos = merge(pheno_ids, gene.pos, by="ids", all.x=TRUE, all.y=FALSE)
# Only keep relevant columns
gene.pos = gene.pos[,1:4]

# Write the list of genetic positions to CSV
write.csv(gene.pos, "./processed/lowry_gene_pos.csv", row.names = FALSE)

###############################################################################

## Create a map between physical and genetic positions

# Read in list of gene positions if it hasn't already been loaded
if (!exists("gene.pos")) {
  gene.pos = read.csv("./processed/lowry_gene_pos.csv", header=TRUE)
}

# TKrils_Marker_PhysPos.csv contains physical positions (provided)
TKrils = read.csv("./processed/TKrils_Marker_PhysPos.csv", header=TRUE)
# Split by chromosome
TKrils_chr = lapply(1:5, function(i){TKrils[TKrils$Chr==i, ]}) 

# Get the positions of the markers (from the genoprobs) in centimorgans
markers_cM = lapply(cross.ge.genoprobs$geno, function(x){x$map})

# Get the positions of the phenotype transcripts in base pairs
trans_bp = lapply(1:5, function(i){
  gene.pos$start[gene.pos$sequence==paste("Chr",i,sep="")]})
# Get the positions in centimorgans
trans_cM = lapply(1:5, function(i){
  approx(TKrils_chr[[i]]$Col.Phys.Pos, 
         TKrils_chr[[i]]$Gen..Pos..LargeMap., 
         trans_bp[[i]], rule=2)})

# Re-center the centimorgan positions as though the 5 chromosomes were lined 
# up adjacent to each other on a single continuous scale
# Endpoints for the 5 chromosomes
chr_cM = c(0, cumsum(sapply(cross.ge.genoprobs$geno, function(x){
  max(x$map)})))
# Marker positions, shifted based on chromosome
cumsum_markers_cM = do.call(c, lapply(1:5, function(i){
  markers_cM[[i]] + chr_cM[i]}))
# Transcript positions, shifted based on chromosme
cumsum_trans_cM = do.call(c, lapply(1:5, function(i){
  trans_cM[[i]]$y + chr_cM[i]}))
# Double each transcript position (for dry/wet)
cumsum_trans_cM2 = rep(cumsum_trans_cM, each=2)

# Save cM positions
save(chr_cM, cumsum_markers_cM, cumsum_trans_cM, cumsum_trans_cM2, 
     file="./processed/lowry_plot_data.rData")

###############################################################################

## Plot results

# Load physical distances if they haven't already been loaded
if (!exists("chr_cM") || !exists("cumsum_markers") 
    || !exists("cumsum_trans_cM2")) {
  load("./processed/lowry_plot_data.rData")
}

# Read in L1 coefficient estimates
lowry_l1_coeffs = read.csv("./processed/lowry_l1_coeffs.csv", header=FALSE)
# Pull out coefficients corresponding to lambda of interest
lambda = 4 # 4th lambda has value 1.728
coeffs = matrix(lowry_l1_coeffs[,lambda], 452, 51326)

# Remove the intercepts, the cyto (X) contrast, and dry/wet environment 
# (Z) contrast
coeffs = coeffs[-(1:2), -(1:2)]
# Find the indices of the non-zero interactions
idx = which(coeffs != 0, arr.ind=TRUE)

# Indices of the main effects
idx_mains = idx[idx[,2] %% 2 != 0,] 
# Indices of the interactions 
idx_inter = idx[idx[,2] %% 2 == 0,] 

# Chromosome names used for labeling
chr_names = paste(paste("Chr", 1:5, sep=""), 
                  collapse = paste(replicate(14, " "), collapse = ""))

# Map the coefficients to their physical positions
# Main effects
mains = cbind(cumsum_markers_cM[idx_mains[,1]], 
              cumsum_trans_cM2[idx_mains[,2]])
# Interactions
inter = cbind(cumsum_markers_cM[idx_inter[,1]], 
              cumsum_trans_cM2[idx_inter[,2]])


png("./pictures/lowry_gene_vs_qtl_pos_lambda1.7.png", width=380, height=380)
par(mar=c(4.1,4.1,1.1,1.1))

# Main effects
plot(mains[,1], mains[,2], col="royalblue4", 
     xlab="QTL Position (cM)", ylab="Gene Position (cM)") 

# Interactions
points(inter[,1], inter[,2], pch=15, col="firebrick3")

# Reference lines to delineate chromosomes
invisible(sapply(chr_cM[2:5], function(x){abline(v=x)}))
invisible(sapply(chr_cM[2:5], function(x){abline(h=x)}))

# Label chromosomes
mtext(side = 1, text = chr_names, line = 2)
mtext(side = 2, text = chr_names, line = 2) 
dev.off()
