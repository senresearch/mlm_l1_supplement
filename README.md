Repository for storing code and data to include in the supplement for the L1 
paper. 

---

Agren, J., Oakley, C. G., McKay, J. K., Lovell, J. T. and Schemske, D. W. 
    (2013), ‘Genetic mapping of adaptation reveals fitness tradeoffs in 
    Arabidopsis thaliana’, Proceedings of the National Academy of Sciences 
    110(52), 21077–21082

Data downloaded from: [https://datadryad.org/resource/doi:10.5061/dryad.77971]. 
Currently using `cross.ge.rda`  

- `agren_preprocess.R`: Preprocess data

- `agren.jl`: Run L1-penalized matrix linear model

- `agren.R`: Reproduce plot used in manuscript

---

Lowry, D. B., Logan, T. L., Santuari, L., Hardtke, C. S., Richards, J. H., 
    DeRose-Wilson, L. J., McKay, J. K., Sen, S. and Juenger, T. E. (2013), 
    ‘Expression quantitative trait locus mapping across water availability 
    environments reveals contrasting associations with genomic features in 
    arabidopsis’, The Plant Cell 25(9), 3266–3279.

Data downloaded from: [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42408] 
and supplemental data set 5 in the text. Requires annotation file downloaded 
here: 	 https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff 
Currently using `agren2013_fullGenoMatrix4qtl_withParents.csv`. 

- `lowry_preprocess.R`: Preprocess data

- `lowry.jl`: Run L1-penalized matrix linear model

- `lowry.R`: Reproduce plot used in manuscript

- `lowry_multipleqtl.R`: Run multiple QTL analysis on individual phenotypes

---

Woodruff, T. J., Zota, A. R. and Schwartz, J. M. (2011), ‘Environmental 
    chemicals in pregnant women in the United States: NHANES 2003-2004’, 
    Environmental health perspectives 119(6), 878.

Ran simulations based on this study. 

- `woodruff_sim.jl`: Run L1-penalized matrix linear model

- `woodruff_sim.R`: Reproduce plot used in manuscript

---

`runtime_scaling_dimensions.jl` provides to code to examine how run time 
increases with dimension size in simulated data. 


Helper functions for the analysis: 

- `dummy_fun.jl`: Functions to generate contrasts/dummy variables for 
categorical data

- `sim_funs.jl`: Functions to simulate data
