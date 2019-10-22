# Sparse matrix linear models for structured high-throughput data

This repository contains code and data that can be used to reproduce
the results presented in the paper.

## Simulations examining dependence of runtimes on data size

- [`scaling_times.jl`](code/scaling_times.jl): examine how runtime increases 
with dimension size in simulated data for FISTA with backtracking and ADMM
- [`scaling_times.R`](code/scaling_times.R): visually compare runtimes from 
FISTA with backtracking and ADMM

## Simulations inspired by environmental screening data (Woodruff)

---
Woodruff, T. J., Zota, A. R. and Schwartz, J. M. (2011),
'Environmental chemicals in pregnant women in the United States:
NHANES 2003-2004', Environmental health perspectives 119(6), 878.
---

Simulated data inspired by the structure of the above paper.

- [`woodruff_sim.jl`](code/woodruff_sim.jl): run L1-penalized matrix linear model
- [`woodruff_sim.R`](code/woodruff_sim.R): reproduce plot used in manuscript

## Arabidopsis fitness adaptation QTL data (Ågren)

---
Ågren, J., Oakley, C. G., McKay, J. K., Lovell, J. T., & Schemske, D. W. (2013). Genetic mapping of adaptation reveals fitness tradeoffs in Arabidopsis thaliana. *Proceedings of the National Academy of Sciences*, 110(52), 21077-21082.
---

Download `RIL_DataForSelectionAnalyses3yrs.xls` and `geno.csv` from 
[https://datadryad.org/resource/doi:10.5061/dryad.77971](https://datadryad.org/resource/doi:10.5061/dryad.77971) [^fn1] [^fn2]. 

[^fn1]: Ågren, J., Oakley, C. G., Lundemo, S., & Schemske, D. W. (2017). Adaptive divergence in flowering time among natural populations of Arabidopsis thaliana: estimates of selection and QTL mapping. *Evolution*, 71(3), 550-564.

[^fn2]: Ågren, J., Oakley, C. G., Lundemo, S., & Schemske, D. W. (2016), Adaptive divergence in flowering time among natural populations of Arabidopsis thaliana: estimates of selection and QTL mapping. Data from: Dryad Digital Repository. https://doi.org/10.5061/dryad.77971.

- [`agren_preprocess.R`](code/agren_preprocess.R): preprocess data
- [`agren.jl`](code/agren.jl): run L1-penalized matrix linear model
- [`agren.R`](code/agren.R): reproduce plot used in manuscript
- [`agren_times.jl`](code/agren_times.jl): compare runtimes for different
  L1-penalized algorithms

## Arabidopsis eQTL experiment data (Lowry)

---
Lowry, D. B., Logan, T. L., Santuari, L., Hardtke, C. S., Richards, J. H., DeRose-Wilson, L. J., ... & Juenger, T. E. (2013). Expression quantitative trait locus mapping across water availability environments reveals contrasting associations with genomic features in Arabidopsis. *The Plant Cell*, 25(9), 3266-3279.
---

Download the series matrix file from the paper's GEO site
[https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42408] 
and Supplemental Dataset 1b from 
[http://www.plantcell.org/content/27/4/969/tab-figures-data] [^fn3]. 

[^fn3]: Lovell, J. T., Mullen, J. L., Lowry, D. B., Awole, K., Richards, J. H., Sen, S., ... & McKay, J. K. (2015). Exploiting differential gene expression and epistasis to discover candidate genes for drought-associated QTLs in Arabidopsis thaliana. *The Plant Cell*, 27(4), 969-983.

Download annotation file here:  
[https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff](https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff) .

The files for the marker positions, 
[`TKrils_Marker_PhysPos.csv`](data/TKrils_Marker_PhysPos.csv), and the key 
for the TxK RIL IDs, [`dd2014_cytocovar.csv`](data/dd2014_cytocovar.csv), are 
provided.

- [`lowry_preprocess.R`](code/lowry_preprocess.R): preprocess data
- [`lowry.jl`](code/lowry.jl): run L1-penalized matrix linear model
- [`lowry.R`](code/lowry.R): reproduce plot used in manuscript
- [`lowry_multipleqtl.R`](code/lowry_multipleqtl.R): run multiple QTL
  analysis on individual phenotypes
