phenotypeVectorFilename = paste(outputFolderName, filePrefix, "CategoricalPhenotypeVector.rds",sep="") #make a filename based on the prefix
phenotypeVector = readRDS(phenotypeVectorFilename)

lowCategoryGeneDropper = function(mainTrees, phenotypeVector){
  genesToDrop = vector()
  for(i in 1:length(mainTrees$trees)){
    currentTree = mainTrees$trees[[i]]
    currentTreeName = names(mainTrees$trees[i])
    currentTips = currentTree$tip.label
    
    phenotypedTips = currentTips[which(currentTips %in% names(phenotypeVector))]
    
    phenotypeValues = phenotypeVector[match(phenotypedTips, names(phenotypeVector))]
    phenotypeNumbers = table(phenotypeValues)
    if(any(phenotypeNumbers <3)){
      message(currentTreeName)
      print(phenotypeNumbers)
      genesToDrop = append(genesToDrop, currentTreeName)
    }
  }
  return(genesToDrop)
}  
