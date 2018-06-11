# disease-similarity-fusion

This repository contains all data and code necessary to faithfully recreate the analysis presented in 'Understanding and predicting disease relationships through similarity fusion' (Oerton et al., 2018).
 
A worked example is given [here](https://rawgit.com/e-oerton/disease-similarity-fusion/master/Disease_Similarity_Fusion_Main.html).

## Instructions

All code and evaluation functions are supplied in markdown format ([Disease Similarity Fusion Main.Rmd](https://github.com/e-oerton/disease-similarity-fusion/blob/master/Disease%20Similarity%20Fusion%20Main.Rmd)).

The output of this file is also supplied as a [premade .html file](https://rawgit.com/e-oerton/disease-similarity-fusion/master/Disease_Similarity_Fusion_Main.html), so code and results can be viewed without downloading and running the code locally.  

To download and run the full code, navigate to the directory you wish to install disease-similarity-fusion in and in Mac/Linux terminal run ```git clone https://github.com/e-oerton/disease-similarity-fusion/``` (recommended) or download/extract the zip from GitHub webpage.  Note the repository data files are large (488 Mb).

The disease similarity fusion code is written in R (https://www.r-project.org/) and tested using RStudio. 
The following R packages are required to run 'Disease Similarity Fusion Main.rmd':
1. *limma* (Bioconductor package) for quantile normalization
2. *plotly* for displaying interactive graphs
3. *rmarkdown* to run the .Rmd file
4. *gplots* to display the similarity heatmap
5. *RColorBrewer* to display the similarity heatmap 

The following packages are optional (the .Rmd file can be run without installing these packages, but they are required for certain optional analyses where indicated within the .Rmd file):
1. *igraph* to easily output the graph to e.g. Cytoscape
2. *DOSE* (Bioconductor package) to compute a new ontological similarity matrix
3. *randomForest* to rerun evaluation of classification performance of disease similarities
4. *ROCR* to rerun evaluation of classification performance of disease similarities

To install all packages, run [this script](https://github.com/e-oerton/disease-similarity-fusion/blob/master/Scripts/Package%20Installation%20Commands.R).

To run the full code, open the R Markdown file (making sure the R working directory is set to disease-similarity-fusion: this should happen automatically on opening the Rmd file) and then knit the file (e.g. use the 'Knit' button in RStudio), which will illustrate the steps necessary to create the disease map, and display the output (initially as a hierarchical cluster plot; instructions to output the map in a form that can be used with graph analysis software such as Cytoscape are also given) and evaluation functions which generate the figures used in the manuscript.  

Alternatively, the R script ['Perform Similarity Fusion'](https://github.com/e-oerton/disease-similarity-fusion/blob/master/Scripts/Perform%20Similarity%20Fusion.R) is the key function that performs the actual similarity fusion on input similarity matrices (example matrices are given in the ['Similarity Matrices' folder](https://github.com/e-oerton/disease-similarity-fusion/blob/master/Data/Similarity%20Matrices)).  To simply run 'Perform Similarity Fusion' based on the instructions given in 'Disease Similarity Fusion Main', *limma* is the only package that is required.  
