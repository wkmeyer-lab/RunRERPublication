makePhenMasterTree = function(geneName = NA, filePrefix, convertToCommon = F, tipCol = "tipName", manualPhenotypeTreeLocation = NULL, foregroundCategory = NULL, binarizeTree = NULL){
  binarizeTree = as.logical(binarizeTree)
  useManualTree = F
  if(!is.null(manualPhenotypeTreeLocation)){useManualTree = T}
  outputFolderName = paste("Output/",filePrefix,"/", sep = "")
  
    
  index = geneName
  if(!exists("mainTrees")){mainTrees = readRDS(mainTreesLocation)}
  masterTree = mainTrees$masterTree
  masterTree$node.label = NULL
  
  if(useManualTree){
    phenotypeTree = readRDS(manualPhenotypeTreeLocation)
  }else{
    #pathsFilename = paste(outputFolderName, filePrefix, "CategoricalPathsFile.rds", sep= "") #make a filename based on the prefix
    #pathsObject = readRDS(pathsFilename)
    #pathsTree = paths2Tree(mainTrees, pathsObject, index)
    #paths tree actually currently being unused because of how the master tree phenotype matching works. Because it's relying on the trees beingthe same shape and therefore having matching node numbers, I can't use the paths -- or, at least, it's very messy to try, so I'm not.
    
    phenotypeTreeCategoricalLocation = paste(outputFolderName, filePrefix, "CategoricalTree.rds", sep="") #make a filename based on the prefix
    phenotypeTreeBinaryLocation = paste(outputFolderName, filePrefix, "BinaryTree.rds", sep="") #make a filename based on the prefix
    if(file.exists(phenotypeTreeCategoricalLocation)){
      fullPhenotypeTree = readRDS(phenotypeTreeCategoricalLocation)
    }else if(file.exists(phenotypeTreeBinaryLocation)){
      fullPhenotypeTree = readRDS(phenotypeTreeBinaryLocation)
    }else{
      stop("The prefix has neither a categorical or binary phenotype tree")
    }
    
  }
  phenotypeTree = fullPhenotypeTree
  
  
  phenMasterTree = masterTree
  phenMasterTree = drop.tip(phenMasterTree, phenMasterTree$tip.label[!phenMasterTree$tip.label %in% phenotypeTree$tip.label])
  
  #add category as label to nodes
  allLabels = rep("", (length(phenMasterTree$tip.label)+phenMasterTree$Nnode))
  for(i in 1:length(allLabels)){
    parentEdge = which(phenMasterTree$edge[,2]==i)
    if(!length(parentEdge)==0){
      allLabels[i] = parentEdge
      allLabels[i] = phenotypeTree$edge.length[parentEdge]
    }
  }
  
  if(!binarizeTree){
    allLabels[which(allLabels == foregroundCategory)] = "Foreground"
  }
  if(binarizeTree){
    allLabels[which(allLabels != foregroundCategory)] = "Background"
  }
  
  tipLabels = allLabels[c(1:length(phenotypeTree$tip.label))]
  originalTipValues = phenMasterTree$tip.label
  originalTipValues <<- phenMasterTree$tip.label
  if(convertToCommon){
    source("Src/Reu/ZoonomTreeNameToCommon.R")
    phenMasterTree$tip.label = ZonomNameConvertVectorCommon(phenMasterTree$tip.label, tipColumn = tipCol)
  }
  
  phenMasterTree$tip.label = paste0(phenMasterTree$tip.label, "{", tipLabels, "}")
  internalLabels = allLabels[-c(1:length(phenotypeTree$tip.label))]
  phenMasterTree$node.label = paste0("{", internalLabels, "}")
  phenMasterTree

}
