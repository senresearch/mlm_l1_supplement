# Sparse matrix linear models for structured high-throughput data

This repository contains code to reproduce the results presented in the paper ["Sparse matrix linear models for structured high-throughput data"](https://arxiv.org/abs/1712.05767). 

Analysis was primarily performed in [Julia](https://julialang.org)<sup>[1](#myfootnote1)</sup>, and visualizations were created using [R](https://www.r-project.org/)<sup>[2](#myfootnote2)</sup>. The Julia package associated with this paper is [`MatrixLMnet`](https://github.com/senresearch/MatrixLMnet.jl), which extends the [`MatrixLM`](https://github.com/senresearch/MatrixLM.jl) package. 


## Simulations examining dependence of runtimes on data size

- [`scaling_times.jl`](code/scaling_times.jl): Examine how runtime increases with dimension size in simulated data for two algorithms, FISTA with backtracking and ADMM. 
- [`scaling_times.R`](code/scaling_times.R): Visually compare runtimes for FISTA with backtracking and ADMM. 


## Simulations inspired by environmental screening data (Woodruff et al., 2011)

Woodruff, T. J., Zota, A. R., & Schwartz, J. M. (2011). Environmental chemicals in pregnant women in the United States: NHANES 2003–2004. *Environmental health perspectives*, 119(6), 878-885.

- [`woodruff_sim.jl`](code/woodruff_sim.jl): Simulate data and run L<sub>1</sub>-penalized matrix linear model. 
- [`woodruff_sim.R`](code/woodruff_sim.R): Run univariate linear models on simulated data, and reproduce ROC plot comparing approaches. 


## E. coli genetic screening data (Nichols et al., 2011)

Nichols, R. J., Sen, S., Choo, Y. J., Beltrao, P., Zietek, M., Chaba, R., Lee, S., Kazmierczak, K. M., Lee, K. J., Wong, A., et al. (2011). Phenotypic landscape of a bacterial cell. *Cell*, 144(1):143–156. 

Download the data [here](https://figshare.com/s/f7da693dee83595eafd7)<sup>[3](#myfootnote3)</sup>. Once downloaded, the files should be saved in the `data/raw_KEIO_data/` directory. 

- [`nichols_preprocess.R`](code/nichols_preprocess.R): Preprocess data. 
- [`nichols.jl`](code/nichols.jl): Run L<sub>1</sub>-penalized matrix linear model, with cross-validation. 
- [`nichols.R`](code/nichols.R): Reproduce dot plot and ROC plot for analyzing auxotrophs. 
- [`nichols_sim.jl`](code/nichols_sim.jl): Simulate data and run matrix linear models (least squares and L<sub>1</sub>-penalized). 
- [`nichols_sim.R`](code/nichols_sim.R): Reproduce ROC plots for comparing matrix linear models with and without the L<sub>1</sub> penalty. 


## Arabidopsis fitness adaptation QTL data (Ågren et al., 2013)

Ågren, J., Oakley, C. G., McKay, J. K., Lovell, J. T., & Schemske, D. W. (2013). Genetic mapping of adaptation reveals fitness tradeoffs in Arabidopsis thaliana. *Proceedings of the National Academy of Sciences*, 110(52), 21077-21082.

Download `RIL_DataForSelectionAnalyses3yrs.xls` and `geno.csv` [here](https://datadryad.org/resource/doi:10.5061/dryad.77971)<sup>[4](#myfootnote4), [5](#myfootnote5)</sup>. Once downloaded, the files should be saved in the `data/` directory. 

- [`agren_preprocess.R`](code/agren_preprocess.R): Preprocess data. 
- [`agren.jl`](code/agren.jl): Run L<sub>1</sub>-penalized matrix linear model, with cross-validation. 
- [`agren.R`](code/agren.R): Reproduce plot of main effect and interaction QTL peaks. 
- [`agren_times.jl`](code/agren_times.jl): Compare the runtimes for different L<sub>1</sub>-penalized algorithms. 


## Arabidopsis eQTL experiment data (Lowry et al., 2013)

Lowry, D. B., Logan, T. L., Santuari, L., Hardtke, C. S., Richards, J. H., DeRose-Wilson, L. J., ... & Juenger, T. E. (2013). Expression quantitative trait locus mapping across water availability environments reveals contrasting associations with genomic features in Arabidopsis. *The Plant Cell*, 25(9), 3266-3279.

Download the series matrix file from the paper's GEO site [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42408) and Supplemental Dataset 1b [here](http://www.plantcell.org/content/27/4/969/tab-figures-data)<sup>[6](#myfootnote6)</sup>. Once downloaded, the files should be saved in the `data/` directory. 

Download the annotation file [here](https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff). Once downloaded, the file should be saved in the `data/` directory. 

The files for the marker positions, [`TKrils_Marker_PhysPos.csv`](data/TKrils_Marker_PhysPos.csv), and the key for the TxK RIL IDs, [`dd2014_cytocovar.csv`](data/dd2014_cytocovar.csv), are provided in the `data` subdirectory.

- [`lowry_preprocess.R`](code/lowry_preprocess.R): Preprocess data. 
- [`lowry.jl`](code/lowry.jl): Run L<sub>1</sub>-penalized matrix linear model. 
- [`lowry.R`](code/lowry.R): Reproduce scatterplot of main effect and interaction QTL. 
- [`lowry_multipleqtl.R`](code/lowry_multipleqtl.R): Run multiple QTL analysis on individual phenotypes using the [R/qtl](https://rqtl.org/) package<sup>[7](#myfootnote7)</sup>.
- [`lowry_times.jl`](code/lowry_times.jl): Compare the runtimes for different L<sub>1</sub>-penalized algorithms on random subsets of the data. 
- [`lowry_times.R`](code/lowry_times.R): Plot runtimes of each L<sub>1</sub>-penalized algorithm against number of genes in subset. 


---

<a name="myfootnote1">1</a>. Bezanson, J., Edelman, A., Karpinski, S., and Shah, V. B. (2017). Julia: A fresh approach to numerical computing. *SIAM review*, 59(1):65–98. 

<a name="myfootnote2">2</a>. R Core Team (2018). *R: A Language and Environment for Statistical Computing*. R Foundation for Statistical Computing, Vienna, Austria.

<a name="myfootnote3">3</a>. Nichols, R. J., Sen, S., Choo, Y. J., Beltrao, P., Zietek, M., Chaba, R., Lee, S., Kazmierczak, K. M., Lee, K. J., Wong, A., et al. (2011). Phenotypic landscape of a bacterial cell. *Cell*, 144(1):143–156. 

<a name="myfootnote4">4</a>. Ågren, J., Oakley, C. G., Lundemo, S., & Schemske, D. W. (2017). Adaptive divergence in flowering time among natural populations of Arabidopsis thaliana: estimates of selection and QTL mapping. *Evolution*, 71(3), 550-564.

<a name="myfootnote5">5</a>. Ågren, J., Oakley, C. G., Lundemo, S., & Schemske, D. W. (2016), Adaptive divergence in flowering time among natural populations of Arabidopsis thaliana: estimates of selection and QTL mapping. Data from: Dryad Digital Repository. https://doi.org/10.5061/dryad.77971.

<a name="myfootnote6">6</a>. Lovell, J. T., Mullen, J. L., Lowry, D. B., Awole, K., Richards, J. H., Sen, S., ... & McKay, J. K. (2015). Exploiting differential gene expression and epistasis to discover candidate genes for drought-associated QTLs in Arabidopsis thaliana. *The Plant Cell*, 27(4), 969-983.

<a name="myfootnote7">7</a>. Broman, K. W., Wu, H., Sen, Ś., & Churchill, G. A. (2003). R/qtl: QTL mapping in experimental crosses. *Bioinformatics*, 19(7), 889-890.
