### Perform Similarity Fusion ###
#### Normalize and fuse symmetric similarity matrices

library(limma) # for quantileNormalize function
library(randomForest) # for creating Random Forest classifiers

# Main function:
createFusedMatrix = function(filePathToSimilarityMatrices = "Data/Similarity Matrices", 
                             inputMatrices = c("ontologicalSimilarity","phenotypicSimilarity","litCoOccurrenceSimilarity","geneticSimilarity","transcriptomicSimilarity","drugSimilarity"),
                             randomFeatureVectors = FALSE,normalize = TRUE,weights = c(1,1,1,1,1,1)){
  
  # filePathToSimilarityMatrices: directory where the similarity matrices are stored
  # inputMatrices: object/file names of the matrices to be fused (the R object names should be the same as the filenames minus path and extension)
  # randomFeatureVectors: set to TRUE to use create random similarity matrices based on randomly sampled feature vectors
  # normalize: set to FALSE to skip quantile normalization (not recommended)
  # weights: adjust the weighting of the input matrices to contribute greater or lesser to the fused similarities. Weights need not sum to 1. 
  #          By default the weights are equal (equivalent to taking the mean of similarities in each space), adjusting the weights allows different
  #          spaces to have greater or lesser influence.
  
  if(randomFeatureVectors==FALSE){
    # Load the specified similarity matrices and gather them into a list:
    allFiles = list.files(path = filePathToSimilarityMatrices, pattern="*.Rda")
    for(file in allFiles){
      if(sub(".Rda","",file)%in%inputMatrices){
        load(paste(filePathToSimilarityMatrices,file,sep = "/"))
      }
    }
  } else { # Create randomly sampled feature vectors for 'random' similarity matrices
    
    if("ontologicalSimilarity"%in%inputMatrices){
      source("Creating the Similarity Matrices/Scripts/createDiseaseOntologySimilarityMatrix.R")
      temp <- createDiseaseOntologySimilarityMatrix(permuted = TRUE)
      ontologicalSimilarity <- temp
    }
    
    if("phenotypicSimilarity"%in%inputMatrices){
      source("Creating the Similarity Matrices/Scripts/createPhenotypicSimilarityMatrix.R")
      temp <- createPhenotypicSimilarityMatrix(permuted = TRUE)
      phenotypicSimilarity <- temp
    }
    
    if("litCoOccurrenceSimilarity"%in%inputMatrices){
      source("Creating the Similarity Matrices/Scripts/createLiteratureCoOccurrenceSimilarityMatrix.R")
      temp <- createCoOccurrenceSimilarityMatrix(permuted = TRUE)
      litCoOccurrenceSimilarity <- temp
    }
    
    if("geneticSimilarity"%in%inputMatrices){
      source("Creating the Similarity Matrices/Scripts/createGeneticSimilarityMatrix.R")
      temp <- createGeneticSimilarityMatrix(permuted = TRUE)
      geneticSimilarity <- temp
    }
    
    if("transcriptomicSimilarity"%in%inputMatrices){
      source("Creating the Similarity Matrices/Scripts/createTranscriptomicSimilarityMatrix.R")
      temp <- createTranscriptomicSimilarityMatrix(permuted = TRUE)
      transcriptomicSimilarity <- temp
    }
    
    if("drugSimilarity"%in%inputMatrices){
      source("Creating the Similarity Matrices/Scripts/createDrugSharingSimilarityMatrix.R")
      temp <- createDrugSimilarityMatrix(permuted = TRUE)
      drugSimilarity <- temp
    }
    
  }
  matrixList = lapply(inputMatrices,function(x) get(x))
  names(matrixList) = inputMatrices

  # Quantile normalize the similarity matrices:
  if(normalize){
    matrixListNormalized = normalizeMatrices(matrixList)
  } else {
    matrixListNormalized = matrixList
  }
  names(matrixListNormalized) = inputMatrices
  
  
  # Now combine the normalized matrices:    
  fusedMatrix = fuseMatrices(matrixListNormalized,weights)
  
  return(fusedMatrix)
}



normalizeMatrices = function(matrixListNonNormalized){
  # Create n*k matrix of n similarity scores from each of k spaces:
  similarityScores = as.matrix(sapply(matrixListNonNormalized,function(x) x[lower.tri(x)]))
  # Normalize scores across each space using limma normalizeQuantiles function
  normalizedValues = normalizeQuantiles(similarityScores) 
  for(k in 1:length(matrixListNonNormalized)){
    # Convert each of k vectors of normalized similarity scores back into symmetric matrix
    # Code from: http://r.789695.n4.nabble.com/how-to-convert-the-lower-triangle-of-a-matrix-to-a-symmetric-matrix-td823271.html
    temp = diag(nrow(matrixListNonNormalized[[k]]))
    temp[lower.tri(temp)] = normalizedValues[,k]
    temp <- temp + t(temp)
    diag(temp) <- NA
    rownames(temp) = rownames(matrixListNonNormalized[[k]])
    colnames(temp) = rownames(matrixListNonNormalized[[k]])
    assign(paste0(names(matrixListNonNormalized)[[k]],"Normalized"),temp)
  }
  # Return list of normalized symmetric matrices
  return(lapply(names(matrixListNonNormalized),function(x) get(paste0(x,"Normalized"))))
}



fuseMatrices = function(matrixListToFuse,weights){
  # Intialize empty matrix
  fusedMatrix = matrix(nrow=nrow(matrixListToFuse[[1]]),ncol=nrow(matrixListToFuse[[1]]))
  rownames(fusedMatrix) = rownames(matrixListToFuse[[1]])
  colnames(fusedMatrix) = rownames(matrixListToFuse[[1]])
  
  # For each matrix item, find its mean across the k normalized matrices
  combs = combn(rownames(matrixListToFuse[[1]]),2)
  for(i in 1:ncol(combs)){
    fusedMatrix[combs[,i][1],combs[,i][2]] <- sum(weights*sapply(matrixListToFuse, function(x) x[combs[,i][1],combs[,i][2]]))/sum(weights)
    fusedMatrix[combs[,i][2],combs[,i][1]] <- fusedMatrix[combs[,i][1],combs[,i][2]]
  }
  
  return(fusedMatrix)
}
