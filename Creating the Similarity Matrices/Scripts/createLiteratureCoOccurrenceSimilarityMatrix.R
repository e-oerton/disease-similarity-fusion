### Literature co-occurence matrices from PubMed abstracts
#### Feature vector: top 100 most frequently co-occurring diseases by NPMI score

# This space is based on a matrix of 8690 document-level co-occurrence scores (NPMI) from literature mining of PubMed abstracts by Patrick S.H. Lewis (2016)

createCoOccurrenceSimilarityMatrix = function(permuted=FALSE){
  
  # Load the 8690x8690 co-occurrence scores
  # The colnames are the MeSH IDs, the rownames have the MeSH IDs translated into the equivalent disease names 
  load("Creating the Similarity Matrices/Input Data/fulllitsim.Rda")
  
  # Create feature vector from literature co-occurrence matrix
  co_occurrence = vector("list",length(unique(diseaseDatasetInfo$mesh.name)))
  names(co_occurrence) = unique(diseaseDatasetInfo$mesh.name)
  
  for(disname in unique(diseaseDatasetInfo$mesh.name)){
    # These three diseases have less than 100 co-occurring diseases:
    if(disname%in%c("allergic contact dermatitis","rheumatoid arthritis, systemic juvenile","male infertility")){
      co_occurrence[disname] = list(names(fulllitsim[disname,][which(fulllitsim[disname,]!=0)]))
    } else {
      co_occurrence[disname] = list(names(sort(fulllitsim[disname,],decreasing = TRUE)[1:100])) # No ties in this dataset, unlike with ontological and genetic spaces
    }
  }
  
  # Code to create random matrix
  if(permuted == TRUE){ 
    shuffledLengths = sample(sapply(co_occurrence,length))
    allLitDiseases = sort(unique(unlist(co_occurrence)))
    litDiseaseProbs = table(unlist(co_occurrence))/sum(table(unlist(co_occurrence)))
    for(i in 1:length(co_occurrence)){
      co_occurrence[[i]] <- sample(allLitDiseases,shuffledLengths[[i]],prob = litDiseaseProbs)
    }
  }
  
  
  # Create Jaccard similarity matrix
  combs = combn(unique(diseaseDatasetInfo$mesh.name),2)
  litCoOccurrenceSimilarity = matrix(nrow=nrow(diseaseDatasetInfo),ncol=nrow(diseaseDatasetInfo))
  rownames(litCoOccurrenceSimilarity) = unique(diseaseDatasetInfo$mesh.name)
  colnames(litCoOccurrenceSimilarity) = unique(diseaseDatasetInfo$mesh.name)
  
  for(i in 1:ncol(combs)){
    
    litCoOccurrenceSimilarity[combs[,i][1],combs[,i][2]] = jaccard(co_occurrence[[combs[,i][1]]],co_occurrence[[combs[,i][2]]])
    litCoOccurrenceSimilarity[combs[,i][2],combs[,i][1]] = litCoOccurrenceSimilarity[combs[,i][1],combs[,i][2]]
    
  }
  rm(i)
  diag(litCoOccurrenceSimilarity)<- NA
  
  rownames(litCoOccurrenceSimilarity) = diseaseDatasetInfo$condition
  colnames(litCoOccurrenceSimilarity) = diseaseDatasetInfo$condition
  
  return(litCoOccurrenceSimilarity)
  
}

