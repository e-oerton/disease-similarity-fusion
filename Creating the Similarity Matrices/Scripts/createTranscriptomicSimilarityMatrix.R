### Transcriptomic data from Gene Expression Omnibus
#### Feature vector: top 100 most significantly changed genes by log fold change

# Gene expression data from GEO (see diseaseDatasetInfo for experiment details) in the form of 'topTables' output from Limma

createTranscriptomicSimilarityMatrix = function(permuted = FALSE){
  
  # Load gene expression data
  load("Creating the Similarity Matrices/Input Data/geneExpressionData.Rda")
  sharedGenes = Reduce(intersect,lapply(geneExpressionData,function(x) x$gene))
  geneExpressionData = lapply(geneExpressionData, function(x) x[which(x$gene%in%sharedGenes),])
  
  # Create feature vector:
  # Remove genes with non-significant p values of log fold change:
  geneExpressionDataSignificant = lapply(geneExpressionData, function(x) x[which(x$P.Value<0.05&(x$gene!="")),])
  # Get the top 100 most differentially expressed genes by absolute log fold change
  topregs = lapply(geneExpressionDataSignificant, function (x) x[order(-abs(x$logFC)),][1:100,])
  #Note direction of regulation:
  upregs = lapply(topregs, function (x) unlist(x[which(x$logFC>0),"gene"]))
  downregs = lapply(topregs, function (x) unlist(x[which(x$logFC<0),"gene"]))
  names(upregs) = sapply(names(geneExpressionData),function(x) diseaseDatasetInfo[which(rownames(diseaseDatasetInfo)==x),"condition"]) 
  names(downregs) = names(upregs)
  
  # Code to create random matrix:
  if(permuted==TRUE){
    alltopregs = sort(unique(unlist(sapply(topregs,function(x) x$gene))))
    topregsProb = table(unlist(sapply(topregs,function(x) x$gene)))/sum(table(unlist(sapply(topregs,function(x) x$gene))))
    upregdist = sample(sapply(upregs, length))
    
    for(i in 1:length(geneExpressionData)){
      temp = sample(alltopregs,100,prob = topregsProb)
      upregs[[i]] = sample(temp,upregdist[[i]])
      downregs[[i]] = sample(temp,100-upregdist[[i]])
    }
  }
  
  # Create Jaccard similarity matrices
  combs = combn(unique(diseaseDatasetInfo$condition),2)
  transcriptomicSimilarity = matrix(nrow=nrow(diseaseDatasetInfo),ncol=nrow(diseaseDatasetInfo))
  rownames(transcriptomicSimilarity) = diseaseDatasetInfo$condition
  colnames(transcriptomicSimilarity) = diseaseDatasetInfo$condition
  
  for(i in 1:ncol(combs)){

    jup = jaccard(upregs[[combs[,i][1]]],upregs[[combs[,i][2]]])
    jdown = jaccard(downregs[[combs[,i][1]]],downregs[[combs[,i][2]]])
    jaccardavg = mean(c(length(c(upregs[[combs[,i][1]]],upregs[[combs[,i][2]]]))*jup,length(c(downregs[[combs[,i][1]]],downregs[[combs[,i][2]]]))*jdown),na.rm = TRUE)/100
    
    transcriptomicSimilarity[combs[,i][1],combs[,i][2]] = jaccardavg
    transcriptomicSimilarity[combs[,i][2],combs[,i][1]] = transcriptomicSimilarity[combs[,i][1],combs[,i][2]]
    
  }
  rm(i)
  diag(transcriptomicSimilarity) <- NA 
  
  return(transcriptomicSimilarity)
   
}
