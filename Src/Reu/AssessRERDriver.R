library(RERconverge)
# ----------- Defaults for testing ----------
#mainPrefix = "CategoricalInsVertivoreTree"
#mainPairwise = "Herbivore-Carnivore"
#binaryPrefixOne = "CategoricalBinaryHerbivoreTree"
#binaryPrefixTwo = "CategoricalBinaryCarnivoreTree"
#binaryPhenOne = "Herbivore"
#binaryPhenTwo = "Carnivore"
#strictBoth = F

#mainPrefix = "CategoricalPrunedCarnivoreTree"
#mainPairwise = "Carvniore-Herbivore"
#binaryPrefixOne = "CategoricalBinaryHerbivoreTree"
#binaryPhenOne = "Herbivore"
#binaryPrefixTwo = "CategoricalBinaryCarnivoreTree"
#binaryPhenTwo = "Carnivore"

#mainPrefix = "CategoricalInsVertivoreTree"
#mainPairwise = "Insectivore-Vertivore"
#binaryPrefixOne = "CategoricalBinaryInsectivoreTree"
#binaryPhenOne = "Insectivore"
#binaryPrefixTwo = "CategoricalBinaryVertivoreTree"
#binaryPhenTwo = "Vertivore"

AssessRERDriver = function(mainPrefix, mainPairwise, binaryPrefixOne, binaryPhenOne, binaryPrefixTwo, binaryPhenTwo, strictBoth = F){
  # -- handle pairwises without spaces pre-included --
  if(!grepl(" ", mainPairwise)){
    mainPairwise <- gsub("-", " - ", mainPairwise)
  }
    
  # ---- load the relevant files -----
  mainCorrelationLocation = paste0("Output/",mainPrefix, "/", mainPrefix, "PairwiseCorrelationFile.rds")
  mainCorrelationCombined = readRDS(mainCorrelationLocation)
  pairwisePostion = which(names(mainCorrelationCombined) == mainPairwise)
  mainCorrelation = mainCorrelationCombined[[pairwisePostion]]
  rm(mainCorrelationCombined)
  rm(pairwisePostion)
  
  binaryOneCorrelationLocation = paste0("Output/",binaryPrefixOne, "/", binaryPrefixOne, "PairwiseCorrelationFile.rds")
  binaryOneCorrelationCombined = readRDS(binaryOneCorrelationLocation)
  binaryOneCorrelation = binaryOneCorrelationCombined[[1]]
  rm(binaryOneCorrelationCombined)
  
  binaryTwoCorrelationLocation = paste0("Output/",binaryPrefixTwo, "/", binaryPrefixTwo, "PairwiseCorrelationFile.rds")
  binaryTwoCorrelationCombined = readRDS(binaryTwoCorrelationLocation)
  binaryTwoCorrelation = binaryTwoCorrelationCombined[[1]]
  rm(binaryTwoCorrelationCombined)
  
  
  # ----- check that all genes are identical -- 
  if(!all.equal(rownames(mainCorrelation), rownames(binaryOneCorrelation)) && all.equal(rownames(mainCorrelation), rownames(binaryTwoCorrelation))){
    stop("The genes in the supplied correlations do not match. Aborting.")
  }
  
  colnames(mainCorrelation) = paste0("main_",colnames(mainCorrelation))
  colnames(binaryOneCorrelation) = paste0(binaryPhenOne,"_",colnames(binaryOneCorrelation))
  colnames(binaryTwoCorrelation) = paste0(binaryPhenTwo,"_",colnames(binaryTwoCorrelation))
  
  combinedCorrelations = cbind(mainCorrelation, binaryOneCorrelation, binaryTwoCorrelation)
  combinedCorrelations = combinedCorrelations[,c(1,4,7,2,5,8,3,6,9)]
  
  combinedCorrelations$Driver = rep(NA, nrow(combinedCorrelations))
  combinedCorrelations$DriverNumeric = rep(0, nrow(combinedCorrelations))
  
  for(i in 1:nrow(combinedCorrelations)){
    workingRow = combinedCorrelations[i,]
    if(any(is.na(workingRow[c(1,2,3)]))){
      combinedCorrelations$Driver[i] = "MissingData"
      combinedCorrelations$DriverNumeric[i] = 5
      next
    }
    mainRho = abs(workingRow[1])
    binOneRho = abs(workingRow[2])
    binTwoRho = abs(workingRow[3])
    binSumRho = binOneRho + binTwoRho

    if(strictBoth){
      if(binOneRho > 1.5*binTwoRho){
        combinedCorrelations$Driver[i] = binaryPhenOne
        combinedCorrelations$DriverNumeric[i] = 1
      }else if(binTwoRho > 1.5* binOneRho){
        combinedCorrelations$Driver[i] = binaryPhenTwo
        combinedCorrelations$DriverNumeric[i] = 2
      }else if(mainRho > binSumRho){
        combinedCorrelations$Driver[i] = "DoubleBoth"
        combinedCorrelations$DriverNumeric[i] = 4
      }else if(mainRho > binOneRho && mainRho > binTwoRho){
        combinedCorrelations$Driver[i] = "Both"
        combinedCorrelations$DriverNumeric[i] = 3
      }else{
        combinedCorrelations$Driver[i] = "Unclear"
      }
    }else{    
      if(mainRho > binSumRho){
        combinedCorrelations$Driver[i] = "DoubleBoth"
        combinedCorrelations$DriverNumeric[i] = 4
      }else if(mainRho > binOneRho && mainRho > binTwoRho){
        combinedCorrelations$Driver[i] = "Both"
        combinedCorrelations$DriverNumeric[i] = 3
      }else if(binOneRho > 1.5*binTwoRho){
        combinedCorrelations$Driver[i] = binaryPhenOne
        combinedCorrelations$DriverNumeric[i] = 1
      }else if(binTwoRho > 1.5* binOneRho){
        combinedCorrelations$Driver[i] = binaryPhenTwo
        combinedCorrelations$DriverNumeric[i] = 2
      }else{
        combinedCorrelations$Driver[i] = "Unclear"
      }
    }
  }
  combinedCorrelations
}
