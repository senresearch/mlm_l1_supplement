# Fast algorithms for fitting L1-penalized linear models to matrix-valued structured high-throughput data

This repository contains code and data that can be used to reproduce
the results presented in the paper.

## Runtime and data size

---
Woodruff, T. J., Zota, A. R. and Schwartz, J. M. (2011),
'Environmental chemicals in pregnant women in the United States:
NHANES 2003-2004', Environmental health perspectives 119(6), 878.
---

Simulated data inspired by the structure of the above paper.

- [`scaling_times.jl`](scaling_times.jl) provides the code to examine
how runtime increases with dimension size in simulated data
- [`woodruff_sim.jl`](woodruff_sim.jl): runs L1-penalized matrix linear model
- [`woodruff_sim.R`](woodruff_sim.R): reproduces plot used in manuscript

## Agren data


---
Agren, J., Oakley, C. G., McKay, J. K., Lovell, J. T. and Schemske,
D. W.  (2013), ‘Genetic mapping of adaptation reveals fitness
tradeoffs in Arabidopsis thaliana’, Proceedings of the National
Academy of Sciences 110(52), 21077–21082.
---

Download `RIL_DataForSelectionAnalyses3yrs.xls` and `geno.csv` from 
[https://datadryad.org/resource/doi:10.5061/dryad.77971](https://datadryad.org/resource/doi:10.5061/dryad.77971) [^fn1] [^fn2]. 

[^fn1]: Agren, J., Oakley, C. G., Lundemo, S. and Schemske, D. W. (2017), 
    ‘Adaptive divergence in flowering time among natural populations of
    arabidopsis thaliana: Estimates of selection and qtl mapping’, Evolution 
    71(3), 550–564.

[^fn2]: Agren, J., Oakley, C. G., Lundemo, S. and Schemske, D. W. (2016), 
    ‘Adaptive divergence in flowering time among natural populations of 
    arabidopsis thaliana: estimates of selection and qtl mapping’. Data from: 
    Dryad Digital Repository. https://doi.org/10.5061/ dryad.77971.

- [`agren_preprocess.R`](agren_preprocess.R): preprocess data
- [`agren.jl`](agren.jl): run L1-penalized matrix linear model
- [`agren.R`](agren.R): reproduce plot used in manuscript
- [`agren_times.jl`](agren_times.jl): compare runtimes for different
  L1-penalized algorithms

## Lowry data

---
Lowry, D. B., Logan, T. L., Santuari, L., Hardtke, C. S., Richards, J. H., 
DeRose-Wilson, L. J., McKay, J. K., Sen, S. and Juenger, T. E. (2013), 
'Expression quantitative trait locus mapping across water availability 
environments reveals contrasting associations with genomic features in 
Arabidopsis', The Plant Cell 25(9), 3266–3279.
---

Download the series matrix file from the paper's GEO site
[https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42408] 
and Supplemental Dataset 1b from 
[http://www.plantcell.org/content/27/4/969/tab-figures-data] [^fn3]. 

[^fn3]: Lovell, J. T., Mullen, J. L., Lowry, D. B., Awole, K., Richards, 
    J. H., Sen, S., Verslues, P. E., Juenger, T. E. and McKay, J. K. (2015), 
    'Exploiting differential gene expression and epistasis to discover 
    candidate genes for drought-associated qtls in arabidopsis thaliana', The 
    Plant Cell 27(4), 969{983.

Download annotation file here:  
[https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff](https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff)
and the files `TKrils_Marker_PhysPos.csv` and `dd2014_cytocovar.csv`,
which are provided.

- [`lowry_preprocess.R`](lowry_preprocess.R): preprocess data
- [`lowry.jl`](lowry.jl): run L1-penalized matrix linear model
- [`lowry.R`](lowry.R): reproduce plot used in manuscript
- [`lowry_multipleqtl.R`](lowry_multipleqtl.R): run multiple QTL
  analysis on individual phenotypes
