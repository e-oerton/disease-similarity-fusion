### Genetic data from DisGeNet
#### Feature vector: top 100 genes with strongest evidence of variation associated with the disease

# Genetic variant data from DisGeNet: http://www.disgenet.org/web/DisGeNET/menu/downloads

createGeneticSimilarityMatrix = function(permuted = FALSE){
  
  # Load DisGeNet data
  # Associations of type AlteredExpression were removed to avoid overlap with the transcriptomic feature space
  # Associations to the non-gene NEWENTRY were also removed
  load("Creating the Similarity Matrices/Input Data/disgene.Rda")
  
  # Create feature vector: get genetic variants (and evidence score of gene association) for our 85 diseases
  disgenes = vector("list",length(unique(diseaseDatasetInfo$disgene.name)))
  names(disgenes) = unique(diseaseDatasetInfo$disgene.name)
  for(disease in diseaseDatasetInfo$disgene.name){
   temp = disgene[which(tolower(disgene$diseaseName)==disease),"score"]
   names(temp) = disgene[which(tolower(disgene$diseaseName)==disease),"geneSymbol"]
   disgenes[[disease]] <- temp
  }

  # If there are more than 100 genes associated, keep the top 100 by evidence score
  disgenes = sapply(disgenes,function(x) if (length(x)>100) {thresh = sort(x,decreasing = TRUE)[100];x[which(x>=thresh)]} else x)
  names(disgenes) = tolower(names(disgenes))
  
  # Code to create random matrix
  if(permuted==TRUE){
    alldisgenes = sort(unique(unlist(sapply(disgenes,names))))
    disgeneProbs = table(unlist(sapply(disgenes,names)))/sum(table(unlist(sapply(disgenes,names))))
    shuffledGeneEvidenceScores = sample(disgenes)
    for(i in 1:length(disgenes)){
      temp = shuffledGeneEvidenceScores[[i]]
      names(temp) = sample(alldisgenes,length(temp),prob = disgeneProbs)
      disgenes[[i]] <- temp
    }
  }
  
  # Create Jaccard similarity matrix
  combs = combn(unique(diseaseDatasetInfo$disgene.name),2)
  geneticSimilarity = matrix(nrow=nrow(diseaseDatasetInfo),ncol=nrow(diseaseDatasetInfo))
  rownames(geneticSimilarity) = unique(diseaseDatasetInfo$disgene.name)
  colnames(geneticSimilarity) = unique(diseaseDatasetInfo$disgene.name)
  
  for(i in 1:ncol(combs)){
    
    geneticSimilarity[combs[,i][1],combs[,i][2]]= jaccard(names(disgenes[[combs[,i][1]]]),names(disgenes[[combs[,i][2]]]))
    geneticSimilarity[combs[,i][2],combs[,i][1]] = geneticSimilarity[combs[,i][1],combs[,i][2]]
    
  }
  rm(i)
  diag(geneticSimilarity)<- NA 

  rownames(geneticSimilarity) = diseaseDatasetInfo$condition
  colnames(geneticSimilarity) = diseaseDatasetInfo$condition
  
  return(geneticSimilarity)
  
}