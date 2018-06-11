# disease-similarity-fusion

This repository contains all data and code necessary to faithfully recreate the analysis presented in 'Understanding and predicting disease relationships through similarity fusion' (Oerton et al., 2018).
 
A worked example is given [here](https://github.com/e-oerton/disease-similarity-fusion/blob/master/Disease%20Similarity%20Fusion%20Main.html).

## Instructions

Note: the disease similarity fusion code is written in R (https://www.r-project.org/) and tested using RStudio. 

Navigate to the directory you wish to install disease-similarity-fusion in and in Mac/Linux terminal run ```git clone https://github.com/e-oerton/disease-similarity-fusion/``` (recommended) or download/extract the zip from GitHub webpage.
Note the repository data files are large (488 Mb).

Follow the instructions in [Disease Similarity Fusion Main](https://github.com/e-oerton/disease-similarity-fusion/blob/master/Disease%20Similarity%20Fusion%20Main.Rmd), which illustrates the construction of a disease map using disease similarity fusion.  This file also includes evaluation functions and instructions to output the map in a form that can be used with graph analysis software such as Cytoscape.  

The output of this file is also supplied as a .html file, so code and results can be viewed without necessarily needing to download and run the whole file.  If you wish to run 'Disease Similarity Fusion Main.rmd', the following R packages are required:
1. *limma* (Bioconductor package)
2. *plotly*
3. *rmarkdown*
4. *gplots*
5. *RColorBrewer*
6. *igraph*
7. *randomForest*
8. *ROCR*
9. *DOSE* (Bioconductor package)

Alternatively, the R script ['Perform Similarity Fusion'](https://github.com/e-oerton/disease-similarity-fusion/blob/master/Scripts/Perform%20Similarity%20Fusion.R) is the key function that performs the actual similarity fusion on input similarity matrices (example matrices are given in the ['Similarity Matrices' folder](https://github.com/e-oerton/disease-similarity-fusion/blob/master/Data/Similarity%20Matrices)).  To simply run 'Perform Similarity Fusion' based on the instructions given in 'Disease Similarity Fusion Main', *limma* is the only package that is required.  
