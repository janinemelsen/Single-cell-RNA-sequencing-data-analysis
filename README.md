# Single-cell-RNA-sequencing-data-analysis
R code for the analysis of single-cell RNA sequencing analysis of human bone marrow natural killer cells as published in [Frontiers in Immunology](https://www.frontiersin.org/articles/10.3389/fimmu.2022.1044398/full). In addition a [script](https://github.com/janinemelsen/Single-cell-RNA-sequencing-data-analysis/tree/main/scripts/NKatlas) can be found to analyze publically available data on human circulating NK cells (~50.000 cells).

## Frontiers publication
### Datasets
Melsen: [GSE199411](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199411)

Crinier: [GSE159624](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159624)

Yang: [GSE130430](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130430)


### Reproduction of figures
[Quality control Melsen dataset](https://github.com/janinemelsen/Single-cell-RNA-sequencing-data-analysis/blob/main/scripts/Frontiers_2022/qualitycontrol.Rmd)

[Clustering of Melsen dataset](https://github.com/janinemelsen/Single-cell-RNA-sequencing-data-analysis/blob/main/scripts/Frontiers_2022/clustering.Rmd)

[Subclustering of CD56dimGZMK- population Melsen dataset](https://github.com/janinemelsen/Single-cell-RNA-sequencing-data-analysis/blob/main/scripts/Frontiers_2022/CD56dimsubclustering.Rmd) 

[Integration Yang dataset](https://github.com/janinemelsen/Single-cell-RNA-sequencing-data-analysis/blob/main/scripts/Frontiers_2022/Integration_Yang.Rmd)

[Integration Crinier dataset](https://github.com/janinemelsen/Single-cell-RNA-sequencing-data-analysis/blob/main/scripts/Frontiers_2022/Integration_Crinier.Rmd)

[Integration Melsen, Crinier and Yang dataset](https://github.com/janinemelsen/Single-cell-RNA-sequencing-data-analysis/blob/main/scripts/Frontiers_2022/Integration_Crinier_Yang_Melsen.Rmd)

[Integration Melsen and HCA](https://github.com/janinemelsen/Single-cell-RNA-sequencing-data-analysis/blob/main/scripts/Frontiers_2022/Integration_HCA_Melsen.Rmd)

[Pseudotime analysis Melsen dataset](https://github.com/janinemelsen/Single-cell-RNA-sequencing-data-analysis/blob/main/scripts/Frontiers_2022/pseudotime.Rmd)


### Sce objects
Click [here](https://github.com/janinemelsen/Single-cell-RNA-sequencing-data-analysis/tree/main/sce%20objects) to open the sce objects, as generated in the paper (big integrated objects were excluded). 
The objects with 'sce' in title, refer to the Melsen dataset


## NK atlas peripheral blood
### Datasets
Crinier et al, 2018 [GSE119562](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119562)

Yang et al, 2021 [GSE130430](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130430)

Witkowski et al, 2021 [GSE184329](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184329)

Ruckert et al, 2022 [GSE197037](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197037)
