### Evaluation functions ###

jaccard = function (x,y) length(intersect(x,y))/length(union(x,y))

# Apply a threshold below which all similarity will be set to 0
applySignificanceThreshold = function(input,sig=0.95){
  if(is.null(dim(input))){print("STOP: METHOD DOES NOT ACCEPT A VECTOR")}
  else {
    thresh = sort(input)[sig*length(sort(input))]
    binarizedinput = input
    binarizedinput[binarizedinput<thresh] <- 0
    return(binarizedinput)
  }
}

# Define a matrix which indicates whether two diseases have the same top-level Disease Ontology class:
sameDOClass = matrix(nrow = nrow(diseaseDatasetInfo),ncol = nrow(diseaseDatasetInfo))
rownames(sameDOClass) = diseaseDatasetInfo$condition
colnames(sameDOClass) = diseaseDatasetInfo$condition
for(disease1 in rownames(sameDOClass)){
  for(disease2 in colnames(sameDOClass)){
    if(diseaseDatasetInfo[which(diseaseDatasetInfo$condition==disease1),"disont.toplevel"]==diseaseDatasetInfo[which(diseaseDatasetInfo$condition==disease2),"disont.toplevel"]){
      sameDOClass[disease1,disease2] <- 1
    } else {
      sameDOClass[disease1,disease2] <- 0
    }
  }
}
diag(sameDOClass) = NA 
sameDOClassValues = sameDOClass[lower.tri(sameDOClass)]

# This function sets similarity values to 0 for those links which are in the same Disease Ontology class, therefore returning only novel links
getNovelLinks = function(similarityValues){
  if(!(is.null(dim(similarityValues)))){print("STOP: METHOD DOES NOT ACCEPT A MATRIX")}
  else {
    similarityValues[which(sameDOClassValues==1)] <- 0
    return(similarityValues)
  }
}


# Define function to quantify number of drug sharing relationships at specified significance threshold
getLinksThatAreNovelAndOrShareDrugs = function(proportionCutOff = 1-proportionOfValuesWhichAreSignificantInFullMatrix, drugSimilarityInfo = shareApprovedAndPhaseThreeDrugs,verbose = TRUE){
  
  drugSimilarityValues = drugSimilarityInfo[lower.tri(drugSimilarityInfo)]
  
  ontologicalSimilarityCutOff = applySignificanceThreshold(ontologicalSimilarity,proportionCutOff)
  phenotypicSimilarityCutOff = applySignificanceThreshold(phenotypicSimilarity,proportionCutOff)
  litCoOccurrenceSimilarityCutOff = applySignificanceThreshold(litCoOccurrenceSimilarity,proportionCutOff)
  geneticSimilarityCutOff = applySignificanceThreshold(geneticSimilarity,proportionCutOff)
  transcriptomicSimilarityCutOff = applySignificanceThreshold(transcriptomicSimilarity,proportionCutOff)
  drugSimilarityCutOff = applySignificanceThreshold(drugSimilarity,proportionCutOff)
  fusedMatrixMinusDrugCutOff = applySignificanceThreshold(fusedMatrixMinusDrug,proportionCutOff)
  fusedMatrixMinusDOCutOff = applySignificanceThreshold(fusedMatrixMinusDO,proportionCutOff)
  fusedMatrixCutOff = applySignificanceThreshold(fusedMatrix,proportionCutOff)
  
  allMatricesCutOff = list(ontologicalSimilarityCutOff[lower.tri(ontologicalSimilarityCutOff)], phenotypicSimilarityCutOff[lower.tri(phenotypicSimilarityCutOff)], litCoOccurrenceSimilarityCutOff[lower.tri(litCoOccurrenceSimilarityCutOff)], geneticSimilarityCutOff[lower.tri(geneticSimilarityCutOff)], transcriptomicSimilarityCutOff[lower.tri(transcriptomicSimilarityCutOff)],drugSimilarityCutOff[lower.tri(drugSimilarityCutOff)], fusedMatrixMinusDrugCutOff[lower.tri(fusedMatrixMinusDrugCutOff)], fusedMatrixMinusDOCutOff[lower.tri(fusedMatrixMinusDOCutOff)], fusedMatrixCutOff[lower.tri(fusedMatrixCutOff)])
  
  novelMatricesCutOff = lapply(allMatricesCutOff, function(x) getNovelLinks(x))
  
  meanJaccardScores = sapply(allMatricesCutOff,function(x) mean(drugSimilarityValues[which(x!=0)]))
  proportionThatAreNovel = sapply(allMatricesCutOff,function(x) length(which(sameDOClassValues[which(x!=0)]==0))/length(which(x!=0)))
  meanJaccardScoresNovel = sapply(novelMatricesCutOff,function(x) mean(drugSimilarityValues[which(x!=0)]))
  
  # This version in the manuscript was evaluated on the full 1000 random matrices, but has been changed to 100 for speed
  meanJaccardScoresRandom = vector("numeric",length(fusedRandomMatrices[1:100]))
  proportionThatAreNovelRandom = vector("numeric",length(fusedRandomMatrices[1:100]))
  meanJaccardScoresNovelRandom = vector("numeric",length(fusedRandomMatrices[1:100]))
  
  for(i in 1:length(fusedRandomMatrices[1:100])){
    aRandomMatrix = applySignificanceThreshold(fusedRandomMatrices[[i]],sig = proportionCutOff)
    aRandomMatrix = aRandomMatrix[lower.tri(aRandomMatrix)]
    tempnovel = getNovelLinks(aRandomMatrix)
    meanJaccardScoresRandom[[i]] = mean(drugSimilarityValues[which(aRandomMatrix!=0)]) 
    proportionThatAreNovelRandom[[i]] = length(which(sameDOClassValues[which(aRandomMatrix!=0)]==0))/length(which(aRandomMatrix!=0))
    if(length(tempnovel[which(tempnovel!=0)])==0){
      meanJaccardScoresNovelRandom[[i]] = NA
    } else { 
      meanJaccardScoresNovelRandom[[i]] = mean(drugSimilarityValues[which(tempnovel!=0)])
    }
  }
  
  if(verbose){
    
    writeLines(c(paste("Mean jaccard scores, DO-Transcriptomic matrices (5 individual matrices):",round(mean(meanJaccardScores[1:5]),3)),
                 paste("Individual jaccard scores, DO-Transcriptomic matrices (5 individual matrices):"),
                 paste(meanJaccardScores[1:5]),
                 paste("Mean jaccard score, kernel created from 5 spaces:",round(meanJaccardScores[7],3))))
    
    writeLines(c(paste("Mean jaccard score, full disease map:",round(meanJaccardScores[9],3)),
                 paste("Mean jaccard score, novel links, full disease map:",round(meanJaccardScoresNovel[9],3))))

  }
  
  return(list(c(mean(meanJaccardScoresRandom),meanJaccardScores),c(mean(meanJaccardScoresNovelRandom),meanJaccardScoresNovel),c(mean(proportionThatAreNovelRandom),proportionThatAreNovel),length(which(fusedMatrixMinusDrugCutOff[lower.tri(fusedMatrixMinusDrugCutOff)]!=0))))
  
}


# Mode function from https://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode
myMode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Create random forest classifier to predict top level Disease Ontology classes
RFDim = function(classToPredict,kernel,nruns = 1000,returnROC=FALSE){
  
  diag(kernel) = 1
  features = as.data.frame(kernel)
  classes = diseaseDatasetInfo$disont.toplevel
  
  binaryclasses <- classes
  binaryclasses[which(binaryclasses!=as.character(classToPredict))]<-0
  binaryclasses[which(binaryclasses==as.character(classToPredict))]<-1
  
  aucs = vector("numeric",nruns)
  FPRs = vector("numeric",nruns)
  TPRs = vector("numeric",nruns)
  for(i in 1:nruns){
    
    sampDisease <- sample(which(binaryclasses==1), 0.8 * nrow(features[which(binaryclasses==1),]))
    sampNonDisease <- sample(which(binaryclasses==0), 0.8 * nrow(features[which(binaryclasses==0),]))
    samp = c(sampDisease,sampNonDisease)
    train <- features[samp, ]
    test <- features[-samp, ]
    
    output.forest <- randomForest(train,as.factor(binaryclasses[samp]),na.action=na.omit,importance=TRUE)
    
    pred <- predict(output.forest, newdata = test,type = "prob")
    pr <- prediction(pred[,2], as.numeric(binaryclasses[-samp]))
    
    auc <- performance(pr, measure = "auc")
    auc <- auc@y.values[[1]]
    aucs[i] <- auc
    
    if(returnROC){
      prf <- performance(pr, measure = "tpr", x.measure = "fpr")
      FPRs[i] <- prf@x.values
      TPRs[i] <- prf@y.values
    }
    
  }
  if(returnROC){
    FPRs = FPRs[which(sapply(FPRs,length)==myMode(sapply(FPRs,length)))]
    FPRs = sapply(1:length(FPRs[[1]]),function(n) mean(sapply(FPRs,function(x)x[n])))
    TPRs = TPRs[which(sapply(TPRs,length)==myMode(sapply(TPRs,length)))]
    TPRs = sapply(1:length(TPRs[[1]]),function(n) mean(sapply(TPRs,function(x)x[n])))
    return(list(aucs,FPRs,TPRs))
  } else {
  return(aucs)
  }
}