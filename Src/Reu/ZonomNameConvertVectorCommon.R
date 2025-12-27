#USAGE: 
# Convert names of one type to another type based on two columns in a specified spreadsheet. 


ZonomNameConvertVectorCommon = function(namesVector, toScientific = F, annotationLocation = "Data/mergedData.csv", scientificColumn = "ScientificName", commonColumn = "CommonName", tipColumn = "tipName"){
  names = namesVector                                                    #make a vector of the names
  manualAnnot = read.csv(annotationLocation)                     #import manual annots file
  tipColumnNumber = which(names(manualAnnot) == tipColumn)
  for(i in 1:length(names)){                                                    #for each name: 
    currentName = names[i]                                                      #use the 'i'th name in the list
    currentRow = manualAnnot[manualAnnot[,tipColumnNumber] %in% currentName, ]             #find a row with the zonom name that matches the current name 
    currentSize = dim(currentRow)                                               #Part of "does row exist check": get the dimensions of the currentRow dataframe
    obsNumber = currentSize[1]                                                  #set "size" equal to the number of observations in 'currentRow'; which is the number of matches to the current name. If none exist it will be 0, if more than one it will be greater than 1. 
    if(obsNumber == 1){                                                         #if only one match exists:
      if(toScientific){
        scientificColumnNumber = which(names(manualAnnot) == scientificColumn)
        currentName = currentRow[[scientificColumnNumber]]                            #get the name from that row
      }else{
        commonColumnNumber = which(names(manualAnnot) == commonColumn)
        currentName = currentRow[[commonColumnNumber]]   
      }
    }else{
      currentName = names[i]                                                    #Otherwise, keep the name the same
    }
    names[i] = currentName                                                      #update the main name list with the name chose
  }
  #colnames(nMatrix) = names                                                     #update the matrix with the new names. 
  return(names)
}