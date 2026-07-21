------------------------------------------------------------------------

# **Cultural affiliation serves to explain spatiotemporal patterns in burial rite practices**

Describing and interpreting spatiotemporal patterns in human culture has been a central focus of anthropology and archaeology for over a century. Recent ethnographic studies have highlighted the complexity of processes generating these patterns, including isolation-by-distance, homophily, and common descent. However, investigating these processes in prehistoric archaeology remains challenging. Using data from the Big Interdisciplinary Archaeological Database (BIAD), linking mortuary information from \~4,200 individuals (\~10,000-2000 BP) to genetic ancestry and mobility from \>1,300 human genomes from Western Eurasia, we analyse the relationship between spatiotemporal patterns in cultural and genomic variation. Using Gaussian process models, we test whether broadly defined genomic affinity clusters correspond to spatiotemporal changes in burial rites, while controlling for other factors.. , . For burial orientation, cultural affiliation was the main explanatory factor, with little to no role for ancestry, whereas body position showed a more mixed pattern in which cultural affiliation remained important. By integrating and modelling these diverse datasets, we provide a detailed understanding of how genomic history intersects with cultural evolution, offering new insights into the dynamics behind these complex processes, and the extent to which genes and culture are transmitted in parallel.

## Overview

This repository contains code and data associated to the manuscript Canteri et al. (2026) *Cultural affiliation accounts for most of the spatiotemporal variation in burial rite practices* **bioRxiv** <https://doi.org/10.64898/2026.05.25.725982>. The R scripts and associated files used for data processing, analysis, modelling, and visualization are organized into separate directories to promote reproducibility, readability, and ease of use.

## Repository Structure

The repository contains a Data folder, a Results folder and a Scripts folder. The .Rmd files in the Scripts folder follow all the different steps, from data processing to data analysis. The INLA models of burial side were run in a HPC computer, following the `FinalModels.R` file. Model performance and cross-validation is run in the files Permutations.R and `Culture_permutations.R`. The generation of figures was done using the `Figures.R` and the `SupplementaryFigures.R` files.

```         
.
├── Burials.Rproj
├── Data/
│   ├── Burial.csv
│   ├── Raw/
│   │   ├── Dataset_S1.csv
│   │   ├── combined.tsv
│   │   ├── mobility_estimates_250y_retrospecive_distance.csv
│   │   └── neo.impute.1000g.sampleInfo_clusterInfo.txt
│   ├── Spatial/
│   │   ├── IntersectionVector.cpg
│   │   ├── IntersectionVector.dbf
│   │   ├── IntersectionVector.prj
│   │   ├── IntersectionVector.shp
│   │   └── IntersectionVector.shx
│   ├── burial rites_culture_period.csv
│   └── culture_keep.RDS
├── README.md
├── Results/
│   ├── Back_Permutations.csv
│   ├── BestModelPerm_Back.csv
│   ├── BestModelPerm_Left.csv
│   ├── BestModelPerm_Right.csv
│   ├── BestModelsAfter5k_results.csv
│   ├── BestModelsBefore5k_results.csv
│   ├── BestModels_Permutations.RDS
│   ├── BestModels_predictions.RData
│   ├── BestModels_results.csv
│   ├── BurialOrientationDifferences_culture.csv
│   ├── BurialOrientationDifferences_culture_ance.csv
│   ├── BurialOrientationMeans_AllIndividuals.csv
│   ├── BurialOrientationModelsFit.csv
│   ├── BurialOrientation_AllIndividuals.RDS
│   ├── Culture/
│   │   ├── BestModel_results.csv
│   │   ├── BestModel_results_old.csv
│   │   ├── CulturePerformance.csv
│   │   ├── CulturePerm_Back.csv
│   │   ├── CulturePerm_Left.csv
│   │   ├── CulturePerm_Right.csv
│   │   ├── Culture_ExplainedVariance.csv
│   │   ├── Culture_Permutations.RDS
│   │   ├── FullModel_results.csv
│   │   ├── FullModel_results_old.csv
│   │   ├── Models_fit.csv
│   │   ├── Models_fit_old.csv
│   │   ├── Models_fit_test2.csv
│   │   └── Results_culture_test2.csv
│   ├── ExplainedVariance.csv
│   ├── FullModelBefore5k_results.csv
│   ├── FullModelNoWHG_results.csv
│   ├── FullModelSteppe_results.csv
│   ├── FullModel_results.csv
│   ├── FullModels_predictions.RData
│   ├── Left_Permutations.csv
│   ├── ModelPerformance.csv
│   ├── ModelPerformanceSummary.csv
│   ├── ModelsFit_BeforeSteppe.csv
│   ├── ModelsFit_Steppe.csv
│   ├── Models_fit.csv
│   ├── Right_Permutations.csv
│   └── noWHGModels_predictions.RData
└── Scripts/
    ├── 01-Graves.Rmd
    ├── 02-DataExploration.rmd
    ├── 03-BurialOrientation.Rmd
    ├── 04-CircularAnalysis.Rmd
    ├── Culture_Permutations.r
    ├── Figures.r
    ├── FinalModels.r
    ├── Permutation.R
    └── SupplementaryFigures.R
```

## Requirements

This project was developed in R.

The required packages are loaded within the individual scripts. If a package is missing, install it using:

``` r
install.packages("package_name")
```

Alternatively, install all required packages before running the analyses.

## Running the Code

The scripts are intended to be run sequentially, executing them in numerical order. Some scripts depend on outputs generated by previous scripts.

## Input Data

All original genetic datasets are placed in the Data/Raw/ directory.

## Output

Results generated by the analyses are saved in the Results/ directory.

## Reproducibility

To ensure reproducibility:

- Use the same version of R whenever possible.
- Install the required R packages before running the scripts.
- Keep the folder structure unchanged so that relative file paths remain valid.

## Notes

- File paths are assumed to be relative to the project directory.
- Individual scripts contain comments describing their purpose and any required inputs.

**Author**

Elisabetta Canteri

University of Copenhagen

elisabetta.canteri\@sund.ku.dk

**License**

[CC-BY-NC-ND 4.0 International license](http://creativecommons.org/licenses/by-nc-nd/4.0/)
