ZonomNameConvertVectorCommonMid = function(namesVector, common = T, manualAnnotLocation = "Data/manualAnnotationsSheet.csv"){
  names = namesVector                                                    #make a vector of the names
  manualAnnot = read.csv(manualAnnotLocation)                     #improt manual annots file
  for(i in 1:length(names)){                                                    #for each name: 
    currentName = names[i]                                                      #use the 'i'th name in the list
    currentRow = manualAnnot[manualAnnot$FaName %in% currentName, ]             #find a row with the zonom name that matches the current name 
    currentSize = dim(currentRow)                                               #Part of "does row exist check": get the dimensions of the currentRow dataframe
    obsNumber = currentSize[1]                                                  #set "size" equal to the number of observations in 'currentRow'; which is the number of matches to the current name. If none exist it will be 0, if more than one it will be greater than 1. 
    if(obsNumber == 1){                                                         #if only one match exists:
      if(common){
        currentName = currentRow$Common.Name.or.Group                             #get the name from that row
      }else{
        currentName = currentRow$Species.Name   
      }
    }else{
      currentName = names[i]                                                    #Otherwise, keep the name the same
    }
    names[i] = currentName                                                      #update the main name list with the name chose
  }
  #colnames(nMatrix) = names                                                     #update the matrix with the new names. 
  return(names)
}

ZoonomTreeNameToCommonOld = function(tree, plot = T, isForegroundTree = T, manualAnnotLocation = "Data/manualAnnotationsSheet.csv", hlcol = "blue", bgcol = "black", fontSize = 0.8, scientific = F){
  
  manualAnnot = read.csv(manualAnnotLocation) 
  inputTree = tree
  tipNames = inputTree$tip.label
  for(i in 1:length(tipNames)){                                                    #for each name: 
    currentName = tipNames[i]                                                      #use the 'i'th name in the list
    currentRow = manualAnnot[manualAnnot$FaName %in% currentName, ]             #find a row with the zonom name that matches the current name 
    currentSize = dim(currentRow)                                               #Part of "does row exist check": get the dimensions of the currentRow dataframe
    obsNumber = currentSize[1]                                                  #set "size" equal to the number of observations in 'currentRow'; which is the number of matches to the current name. If none exist it will be 0, if more than one it will be greater than 1. 
    if(obsNumber == 1){                                                         #if only one match exists:
      if(scientific){
        currentName = currentRow$Tip_Label..Red.is.included.in.CMU.enhancer.dataset..but.missing.alignment.                             #get the name from that row
      }else{
        currentName = currentRow$Common.Name.or.Group                             #get the name from that row 
      }
    }else{
      currentName = tipNames[i]                                                    #Otherwise, keep the name the same
    }
    tipNames[i] = currentName                                                      #update the main name list with the name chose
  }
  inputTree$tip.label = tipNames
  
  if(plot){
    if(isForegroundTree){
      readableTree = inputTree
      readableTree$edge.length[readableTree$edge.length == 0] = 1
      plotTreeHighlightBranches2(readableTree, hlspecies = which(inputTree$edge.length == 1), hlcols = hlcol, bgcol = bgcol, fontSize = fontSize)
    }else{
      plotTree(inputTree)
    }
  }
  return(inputTree)
}