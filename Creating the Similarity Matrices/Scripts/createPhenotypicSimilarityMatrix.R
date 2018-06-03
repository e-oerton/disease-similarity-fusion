### Phenotypic data from Human Phenotype Ontology
#### Feature vector: top 100 most similar phenotypes by NPMI score 

# Disease-phenotype association was calculated from HPO as described in Hoehndorf (2015) Analysis of the human diseasome using phenotype similarity between common, genetic, and infectious diseases

createPhenotypicSimilarityMatrix = function(permuted = FALSE){
  
  # Load the list of associated phenotypes for each disease
  # This is based on filtered-doid-pheno-21.txt: filtered phenotype associations data supplying the highest-ranking 21 phenotypes (based on NPMI)
  # from http://aber-owl.net/aber-owl/diseasephenotypes/data/
  # There are occasional duplicated phenotypes in the supplied data (e.g. where the same term has been used from Human and Mammalian Phenotype Ontologies)
  # This means that there are not exactly 21 phenotypes for each disease

  load("Creating the Similarity Matrices/Input Data/synsetlistforphenotypes.Rda")

  # Code to create random matrix
  if(permuted==TRUE){
    allphenotypes = sort(unique(unlist(synsetlistforphenotypes)))
    phenotypeProbs = table(unlist(synsetlistforphenotypes))/sum(table(unlist(synsetlistforphenotypes)))
    shuffledLengths = sample(sapply(synsetlistforphenotypes,length))
    for(i in 1:length(synsetlistforphenotypes)){
      synsetlistforphenotypes[[i]] <- sample(allphenotypes,shuffledLengths[[i]],prob = phenotypeProbs)
    }
  }
  

  # Create Jaccard similarity matrix
  combs = combn(unique(diseaseDatasetInfo$disont.name),2)
  phenotypicSimilarity = matrix(nrow=nrow(diseaseDatasetInfo),ncol=nrow(diseaseDatasetInfo))
  rownames(phenotypicSimilarity) = unique(diseaseDatasetInfo$disont.name)
  colnames(phenotypicSimilarity) = unique(diseaseDatasetInfo$disont.name)
  
  for(i in 1:ncol(combs)){
    temp1 = synsetlistforphenotypes[[combs[,i][1]]]
    temp2 = synsetlistforphenotypes[[combs[,i][2]]]
    phenotypicSimilarity[combs[,i][1],combs[,i][2]]= jaccard(temp1,temp2) 
    phenotypicSimilarity[combs[,i][2],combs[,i][1]] = phenotypicSimilarity[combs[,i][1],combs[,i][2]]
    
  }
  rm(i)
  diag(phenotypicSimilarity) <- NA 
  
  rownames(phenotypicSimilarity) = diseaseDatasetInfo$condition
  colnames(phenotypicSimilarity) = diseaseDatasetInfo$condition
  
  return(phenotypicSimilarity)
}