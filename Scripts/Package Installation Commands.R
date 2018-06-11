# These commands provide a guideline for package installation
# Depending on your R environment, additional steps may be necessary to install package dependencies

install.packages("plotly")
install.packages("rmarkdown")
#Dependencies of plotly/rmarkdown which may be necessary:
#install.packages("httr")
#install.packages("stringi")
#install.packages("ggplot2")
install.packages("gplots")
install.packages("RColorBrewer")
install.packages("igraph")
install.packages("ROCR")
install.packages("igraph")
install.packages("randomForest")
source("https://bioconductor.org/biocLite.R")
biocLite("DOSE")
biocLite("limma")
