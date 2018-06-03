### Drug data from ChEMBL
#### Feature vector: drugs prescribed or currently or previously in clinical trials

# This information is downloaded from ChEMBL, which lists the maximum clinical trial phase that was reached for a drug for a particular indication.  

createDrugSimilarityMatrix = function(permuted = FALSE){
  
  # Load ChEMBL drug data
  # ChEMBL drug trial phase information downloaded from https://www.ebi.ac.uk/chembl/downloads and mapped to diseases using MeSH term
  # The R object 'drugindicationslist' was created by mapping between the files chembl_drug_indication and chembl_mol_dict obtained from ChEMBL
  load("Creating the Similarity Matrices/Input Data/drugindicationslist.Rda")
  
  # Create feature vector: names of drugs which are approved (max_phase_for_ind >3)
  drugs = sapply(drugindicationslist, function(x) unique(tolower(x[which(x$max_phase_for_ind>3),"pref_name"])))
  
  # Code to create random matrix:
  if(permuted==TRUE){
    alldrugs = sort(unique(unlist(drugs)))
    drugProb = table(unlist(drugs))/sum(table(unlist(drugs)))
    shuffledLengths = sample(sapply(drugs,length))
    for(i in 1:length(drugs)){
      drugs[[i]] = sample(alldrugs,shuffledLengths[[i]],prob = drugProb)
    }
  }


  # Create Jaccard similarity matrix: Jaccard overlap of prescribed drugs
  combs = combn(unique(diseaseDatasetInfo$condition),2)
  drugSimilarity = matrix(nrow=nrow(diseaseDatasetInfo),ncol=nrow(diseaseDatasetInfo))
  rownames(drugSimilarity) = unique(diseaseDatasetInfo$condition)
  colnames(drugSimilarity) = unique(diseaseDatasetInfo$condition)
  
  for(i in 1:ncol(combs)){
    
    drugSimilarity[combs[,i][1],combs[,i][2]] = jaccard(drugs[[combs[,i][1]]],drugs[[combs[,i][2]]])
    drugSimilarity[combs[,i][2],combs[,i][1]] = drugSimilarity[combs[,i][1],combs[,i][2]]
    
  }
  rm(i)
  
  #replace NA values for diseases which don't have any drugs prescribed with 0s
  drugSimilarity[is.na(drugSimilarity)]<-0
  diag(drugSimilarity)<-NA 
  
  return(drugSimilarity)

}