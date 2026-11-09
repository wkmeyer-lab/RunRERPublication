source("Src/Reu/convertColNameToIndex.R")

detectScientificCol = function(inputScientificCol, data){
  if(inputScientificCol == "auto"){
    message("attempting to automatically detect scientific name column.")
    processedColNames = tolower(colnames(data))
    processedColNames = gsub(" ", "", processedColNames)
    processedColNames = gsub("_", "", processedColNames)
    processedColNames = gsub("-", "", processedColNames)
    
    if("scientific" %in% processedColNames){
      indexScientificCol = which(processedColNames == "scientific")
      message("New scientific column set to 'scientific'")
    }else if("scientificname" %in% processedColNames){
      indexScientificCol = which(processedColNames == "scientificname")
      message("New scientific column set to 'scientificname'")
    }else if("scientificnamefull" %in% processedColNames){
      indexScientificCol = which(processedColNames == "scientificnamefull")
      message("New scientific column set to 'scientificnamefull'")
    }else if("species" %in% processedColNames){
      indexScientificCol = which(processedColNames == "species")
      message("New scientific column set to 'species'")
    }else if("scientificbinomial" %in% processedColNames){
      indexScientificCol = which(processedColNames == "scientificbinomial")
      message("New scientific column set to 'scientificbinomial'")
    }else{
      stop("Could not find scientific name column, please manually specify.")
    }
  }else{
    if(!is.integer(inputScientificCol)){
      indexScientificCol = convertColNameToIndex(inputScientificCol, data)
    }else{
      indexScientificCol = inputScientificCol
    }
  }
  indexScientificCol
}


createDebugMatchCol = function(vector){
  vector = gsub(" ", "_", vector) #replace spaces with underscores 
  vector = sub('^([^_]+_[^_]+).*', '\\1', vector) #remove anything  after a second underscore
  vector = tolower(vector)
  
  
}


addColumnToData = function(mainData, newData, colsToAdd = NULL, colRenames = NULL, newDataScientificCol = "auto", mainDataScientificCol = "auto"){
  
  #get the column index of the scientific data 
  newDataScientificCol = detectScientificCol(newDataScientificCol, newData)
  mainDataScientificCol = detectScientificCol(mainDataScientificCol, mainData)
  
  
  #Create columns with identical formatting of scientific names 
  newData$DebugNewMatchingName = createDebugMatchCol(newData[,newDataScientificCol])
  mainData$DebugNewMatchingName = createDebugMatchCol(mainData[,mainDataScientificCol])
  
  #create a column in the new data with the rownumber of the species in the main data
  newData$DebugMatchColumn = match(newData$DebugNewMatchingName, mainData$DebugNewMatchingName)
  newSpecies = newData$DebugNewMatchingName[is.na(newData$DebugMatchColumn)]
  message("Species in new Data not in merged data:")
  message(paste(length(newSpecies), "New Species"))
  if(length(newSpecies) < 50){
    for(i in 1:length(newSpecies)){
      message(newSpecies[i])
    }
  }else{
    message("Too many new species, not printing individual species.")
  }
  
  #Make a data frame with the same species as the main data 
  toAddData = as.data.frame(matrix(NA, ncol = ncol(newData), nrow = nrow(mainData)))
  names(toAddData) = names(newData)
  toAddData$DebugNewMatchingName = mainData$DebugNewMatchingName
  
  #fill the rows of the new datafram with the rows from match species in the new data
  for(i in 1:nrow(toAddData)){
    if(i %in% newData$DebugMatchColumn){
      matchedRow = which(newData$DebugMatchColumn == i)
      toAddData[i,] = newData[matchedRow,]
    }
  }
  
  
  #if specified adding only certain columns, trim the toAddData to those columns
  if(!is.null(colsToAdd)){
    if(length(colsToAdd)==1 & !is.null(colsToAdd)){
      addingSingleColumn = T
    }
    
    #Add a catch to ensure that single additional columns re still treated as dataframes 
    if(addingSingleColumn){
      colsToAdd = c(colsToAdd, colsToAdd)
      if(!is.null(colRenames)){colRenames = c(colRenames, colRenames)}
    }
    
    
    
    #get the indexes of columsn to add
    addColIndexes = integer()
    for(i in 1:length(colsToAdd)){
      addColIndexes[i] = convertColNameToIndex(colsToAdd[i], toAddData)
    }
    toAddData = toAddData[,addColIndexes]
    
    
    if(!is.null(colRenames)){
      if(length(colRenames) != length(colsToAdd)){
        stop("Length of coumn renames must match length of added columns. \n This function replaces the names of ALL added columns. \n To preserve a name from the new data, include that name in the column renames.")
      }else{
        colnames(toAddData) = colRenames
      }
    }
  }else{
    if(!is.null(colRenames)){
      stop("Column renaming only supported when adding specified columns. \n It is not supported when adding all columns.")
    }
  }
  
  #Remove the debug matching column from the main data 
  mainDataDebugMatchColIndex = convertColNameToIndex("DebugNewMatchingName", mainData)
  mainData = mainData[,-mainDataDebugMatchColIndex]
  
  #combine the data, and return the combined DF
  
  if(addingSingleColumn){
    combinedData = cbind(mainData, toAddData[,1])
    names(combinedData)[length(names(combinedData))] = names(toAddData)[1]
  }else{
    combinedData = cbind(mainData, toAddData)
  }
  combinedData
}
