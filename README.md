# Imaging-Based Machine Learning Analysis of Patient-Derived Tumor Organoid Drug Response

R code and data to replicate the figures and analyses in Spiller ER, Ung N, Kim S, Patsch K, Lau R, Strelez C, Doshi C, Choung SJ, Choi B, Juarez Rosales EF, Lenz H-J, Matasci N and Mumenthaler SM (2021) Imaging-Based Machine Learning Analysis of Patient-Derived Tumor Organoid Drug Response. Front. Oncol. 11:771173. [doi: 10.3389/fonc.2021.771173](https://doi.org/10.3389/fonc.2021.771173)

This repository contains the main R Markdown document used to generate figures for the paper. Note that only Figure 2A, 3A, 3C, 3D, 3E, 4A, 4B, 4C and Supplementary Figures S2, S3 and 5B are directly generated (whereas the other figures are composites that integrated images or other elements). 

The file `figures_organoids-vital-status.Rmd` can be `knit`ted to generate the report. Accessory functions (analysis and plotting) are in `functions.R`, that is sourced from the main report.

The consolidated organoids measurements are found in `data/consolidated.csv`. The ground truth and vital dye comparisons are found under `data/ground_truth`.

The project uses [`renv`](https://rstudio.github.io/renv/index.html) for dependency management.
