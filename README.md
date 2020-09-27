# seuratpipeline
R Shiny website for running the R package [Seurat.](https://satijalab.org/seurat/) package to analyse single cell RNA-seq data
Seurat is also hosted on GitHub. You can view the repository at

- https://github.com/satijalab/seurat

For more documentation, visit their page https://satijalab.org/seurat/

## Introduction
Seuratpipeline is a website that uses the R shiny framework. It is user-friendly tool to run a basic analysis on single cell RNA-Seq data using the Seurat package. It reads in the output from 10x and does the pre-processing of the data such as filtering, normalizing and scaling and then proceeds to do the dimensionality reduction, clustering and determining differential markers between cell groups. The expression data, sample data, feature annotation, dimensionality reduction/ clustering, and marker gene information is stored as an RDS object which can then be viewed using [SeuratViewer](https://github.com/Morriseylab/SeuratViewer)

## Requirements
- R (version > 3.4)
- RStudio Server
- Shiny Server (if you need to host it online)
- UMAP python package

If you need help installing Shiny Server or getting started, refer to [this](https://deanattali.com/2015/05/09/setup-rstudio-shiny-server-digital-ocean/#install-r)

For help installing UMAP, refer to this [link](https://github.com/lmcinnes/umap)

## Installation

Run the following command in R to install all required packages
```
install.packages(c("devtools","shiny","shinydashboard","shinyjs","shinyBS","shinyFiles","RColorBrewer",,"ggplot2","dashboardthemes",
                   "dplyr","tidyr","DT","Seurat","visNetwork"))

```

## Input Data 
 Specify the path to the cellranger files as the input. Follow instructions on each page to proceed. For information on the pipeline, refer to this [Seurat documentation](https://satijalab.org/seurat/v3.2/pbmc3k_tutorial.html)
 
## Output data
After processing, the last tab will allow you to save the dataset as an RDS file which can then be uploaded to [SeuratViewer](https://github.com/Morriseylab/SeuratViewer)
