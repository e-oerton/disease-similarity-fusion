### Ontological similarity from Disease Ontology 
#### Feature vector: top 100 most similar DO terms by Resnik similarity

# For each disease, which are the most similar diseases according to the DO ontology structure?  
# Use Lin's measure of semantic similarity, which is based on the information content of the two terms and their most informative common ancestor (see Pesquita 2009 Semantic Similarity in Biomedical Ontologies] (http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000443))

# We use the DOSE package (Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an R/Bioconductor package for Disease Ontology Semantic and Enrichment analysis. Bioinformatics 2015, 31(4):608-609) to measure semantic similarity:
# instructions: https://www.bioconductor.org/packages/release/bioc/vignettes/DOSE/inst/doc/semanticAnalysis.html
createDiseaseOntologySimilarityMatrix = function(permuted = FALSE){
  
  # Load matrix describing similarity between 6553 DO terms, according to Lin's similarity measure (calculated using DOSE package):
  load("Creating the Similarity Matrices/Input Data/linSimilarityDiseaseOntology.Rda")

  # Create feature vector for the 84 diseases in our dataset:
  topLinSimilarities = vector("list",length(unique(diseaseDatasetInfo$disont.name)))
  names(topLinSimilarities) = unique(diseaseDatasetInfo$disont.name)
  for(disname in diseaseDatasetInfo$disont.name){
    linSimilaritySorted = sort(linSimilarity[disname,],decreasing = TRUE)
    thresh = linSimilaritySorted[100]
    topLinSimilarities[disname] = list(names(linSimilaritySorted[which(linSimilaritySorted>=thresh)]))
  }
  
  # The two diseases of class 'syndrome' have no similar diseases: 
  topLinSimilarities[["polycystic ovary syndrome"]] = character(0)
  topLinSimilarities[["irritable bowel syndrome"]] = character(0)
  
  # Code to create random matrices
  if(permuted==TRUE){
    allDOIDs = sort(unique(unlist(topLinSimilarities)))
    DOIDprob = table(unlist(topLinSimilarities))/sum(table(unlist(topLinSimilarities)))
    shuffledLengths = sample(sapply(topLinSimilarities,length))
    for(i in 1:length(topLinSimilarities)){
      topLinSimilarities[[i]] <- sample(allDOIDs,shuffledLengths[[i]],prob = DOIDprob)
    }
  }

  # Create Jaccard similarity matrix
  combs = combn(unique(diseaseDatasetInfo$disont.name),2)
  ontologicalSimilarity = matrix(nrow=nrow(diseaseDatasetInfo),ncol=nrow(diseaseDatasetInfo))
  rownames(ontologicalSimilarity) = unique(diseaseDatasetInfo$disont.name)
  colnames(ontologicalSimilarity) = unique(diseaseDatasetInfo$disont.name)
  
  for(i in 1:ncol(combs)){
    
    ontologicalSimilarity[combs[,i][1],combs[,i][2]] = jaccard(topLinSimilarities[[combs[,i][1]]],topLinSimilarities[[combs[,i][2]]])
    ontologicalSimilarity[combs[,i][2],combs[,i][1]] = ontologicalSimilarity[combs[,i][1],combs[,i][2]]
    
  }
  rm(i)
  # Replace NA values for the two diseases which don't have any ontological partners with 0s
  ontologicalSimilarity[is.na(ontologicalSimilarity)]<-0
  
  diag(ontologicalSimilarity) <- NA 
  
  rownames(ontologicalSimilarity) = diseaseDatasetInfo$condition
  colnames(ontologicalSimilarity) = diseaseDatasetInfo$condition
  

  
  return(ontologicalSimilarity)
}

