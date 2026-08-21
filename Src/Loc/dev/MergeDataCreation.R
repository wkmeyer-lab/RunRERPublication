library(tidyr)
library(tidyverse)
createSortedTable = function(dataSet, scientificNameColumn){
  #Make a column with the scientific name for matching
  inScientificNameColumn = which(names(dataSet) == scientificNameColumn)
  dataSet$Scientific_Binomial = dataSet[,inScientificNameColumn]
  
  #Clean up that name into the same formatting 
  for(i in 1:length(dataSet$Scientific_Binomial)){
    dataSet$Scientific_Binomial[i] = gsub(" ", "_", dataSet$Scientific_Binomial[i]) #replace spaces with underscores 
    dataSet$Scientific_Binomial[i] = sub('^([^_]+_[^_]+).*', '\\1', dataSet$Scientific_Binomial[i]) #remove anything  after a second underscore
    dataSet$Scientific_Binomial[i] = tolower(dataSet$Scientific_Binomial[i])
  }
  dataSet$Scientific_Binomial
  
  dataSet = dataSet[,c(ncol(dataSet), 1:(ncol(dataSet)-1))]
  Divider = rep("-", nrow(dataSet))
  dataSet$Divider = Divider
  dataSet
}
  

CombineDatasets = function(combinedDataInput, newDatasetInput, newDataScientificNameColumn = "ScientificName", newDataCommonNameColumn = NA,  
                            readData = F, addNewSpeciesValue = T, attachAllColumns = T, 
                            nameColumns = NA, manualAddColumns = NA, manualColumnRenames = NA, manualNameColumnRenames = NA,  manualColumnsToIgnore = NA){
  
  combinedData = combinedDataInput
  # -- Set up the two dataframes -- 
  if(readData){
  combinedData = read.csv(combinedDataInput)
  newDataset = read.csv(newDatasetInput)
  }else{
    combinedData = combinedDataInput
    newDataset = newDatasetInput
  }
  
  orderedNewData = orderNewData(newDataset, newDataScientificNameColumn, mergedData = combinedData, addNewSpecies = addNewSpeciesValue, commonNameColumn = newDataCommonNameColumn)
  
  if(addNewSpeciesValue){
    combinedData = addMainSpecies(speciesToAdd = newSpeciesDataset)
  }
  
  # --- attach the columns --- 
  
  if(!attachAllColumns){
    #Data columns
    if(!all(is.na(manualColumnRenames))){
      lengthMatchCheck = length(manualAddColumns) == length(manualColumnRenames)
      if(!lengthMatchCheck){
        stop("The list of main columns to add and renames are not the same length. Makes sure they are, and re-run the script. If you do not want to rename a column, use the columns original name in that postion.")
      }
      usingMainRenames = T
    }else{usingMainRenames = F}
    columnsToAdd = manualAddColumns
    
    if(usingMainRenames){
      for(i in 1:length(columnsToAdd)){
        combinedData = addColumn(orderedNewData, columnsToAdd[i], columnRename = manualColumnRenames[i])
      }
    }else{
      if(!all(is.na(columnsToAdd))){
        for(i in 1:length(columnsToAdd)){
          combinedData = addColumn(orderedNewData, columnsToAdd[i])
        }
      }
    }
    #Name columns 
    if(!all(is.na(nameColumns))){
      if(!all(is.na(manualNameColumnRenames))){
        lengthMatchCheck = length(nameColumns) == length(manualNameColumnRenames)
        if(!lengthMatchCheck){
          stop("The list of name columns to add and renames are not the same length. Makes sure they are, and re-run the script. If you do not want to rename a column, use the columns original name in that postion.")
        }
        usingNameRenames = T
        
      }else{
        usingNameRenames = F
      }
      
      
      if(usingNameRenames){
        for(i in 1:length(nameColumns)){
          combinedData = addColumn(orderedNewData, nameColumns[i], columnRename = manualNameColumnRenames[i], nameColumn = T)
        }
      }else{
        for(i in 1:length(nameColumns)){
          combinedData = addColumn(orderedNewData, nameColumns[i], nameColumn = T)
        }
      }
    }
    
  }else{
    columnsToAdd = names(orderedNewData)
    
    #Making a list of collumn which have already been proccessed by or created by the script 
    collumsAlwaysIgnored = c("DebugNewMatchingName", "X", "DebugMatchColumn")
    collumsAlwaysIgnored = append(collumsAlwaysIgnored, newDataScientificNameColumn)
    collumsAlwaysIgnored = append(collumsAlwaysIgnored, newDataCommonNameColumn)
    collumsAlwaysIgnored = append(collumsAlwaysIgnored, manualColumnsToIgnore)
    
    #Remove those columns from the set to be added 
    columnsToAdd = columnsToAdd[-which(columnsToAdd %in% collumsAlwaysIgnored)]
    
    #Remove any columns specified as name columns, as those will be handled later 
    if(any(columnsToAdd %in% nameColumns)){
    columnsToAdd = columnsToAdd[-which(columnsToAdd %in% nameColumns)]
    }
    #Add non-name columns
    for(i in columnsToAdd){
      combinedData = addColumn(orderedNewData, i)
    }
    
    #Add name columns
    if(!all(is.na(manualNameColumnRenames))){
      lengthMatchCheck = length(nameColumns) == length(manualNameColumnRenames)
      if(!lengthMatchCheck){
        stop("The list of name columns to add and renames are not the same length. Makes sure they are, and re-run the script. If you do not want to rename a column, use the columns original name in that postion.")
      }
      usingNameRenames = T
    }else{usingNameRenames = F}
    
    if(!is.na(nameColumns)){
      if(usingNameRenames){
        for(i in 1:length(nameColumns)){
          combinedData = addColumn(orderedNewData, nameColumns[i], columnRename = manualNameColumnRenames[i], nameColumn = T)
        }
      }else{
        for(i in 1:length(nameColumns)){
          combinedData = addColumn(orderedNewData, nameColumns[i], nameColumn = T)
        }
      }
  }
  }
  combinedData
}  








orderNewData = function(dataSet, scientificNameColumn, mergedData = combinedData, addNewSpecies = F, commonNameColumn = NA){
  #Make a column with the scientific name for matching
  inScientificNameColumn = which(names(dataSet) == scientificNameColumn)
  dataSet$DebugNewMatchingName = dataSet[,inScientificNameColumn]
  
  #Clean up that name into the same formatting 
  for(i in 1:length(dataSet$DebugNewMatchingName)){
    dataSet$DebugNewMatchingName[i] = gsub(" ", "_", dataSet$DebugNewMatchingName[i]) #replace spaces with underscores 
    dataSet$DebugNewMatchingName[i] = sub('^([^_]+_[^_]+).*', '\\1', dataSet$DebugNewMatchingName[i]) #remove anything  after a second underscore
    dataSet$DebugNewMatchingName[i] = tolower(dataSet$DebugNewMatchingName[i])
  }
  dataSet$DebugMatchColumn = match(dataSet$DebugNewMatchingName, mergedData$Scientific_Binomial)
  dataSet = dataSet %>% select(DebugNewMatchingName, everything())
  
  newSpecies = dataSet$DebugNewMatchingName[is.na(dataSet$DebugMatchColumn)]
  message("Species in new Data not in merged data:")
  for(i in 1:length(newSpecies)){
    message(newSpecies[i])
  }
  newData = as.data.frame(matrix(NA, ncol = ncol(dataSet), nrow = nrow(mergedData)))
  names(newData) = names(dataSet)
  newData$DebugNewMatchingName = mergedData$Scientific_Binomial
  
  for(i in 1:nrow(newData)){
    if(i %in% dataSet$DebugMatchColumn){
      matchedRow = which(dataSet$DebugMatchColumn == i)
      newData[i,] = dataSet[matchedRow,]
    }
  }
  
  if(addNewSpecies){
    newSpeciesRows = which(is.na(dataSet$DebugMatchColumn))
    for(i in newSpeciesRows){
      newData = rbind(newData, dataSet[i,])
    }
    
    #Make a global envinroment dataset of the new species with scientific and optionally common names 
    
    newSpeciesData = dataSet[newSpeciesRows,]
    newScientificNameColumn = which(names(newSpeciesData) == scientificNameColumn)
    if(!is.na(commonNameColumn)){
      newCommonNameColumn = which(names(newSpeciesData) == commonNameColumn)
      if(length(newCommonNameColumn) == 0){stop("Specified Common Name Column not found. Make sure it is spelled correctly, and run the code again.")}
    }else{
      newCommonNameColumn = which(names(newSpeciesData) == "DebugMatchColumn") #This is a column that garuntees that all values will be NA, because if they had a match they wouldn't be new.
    }
    names(newSpeciesData)[1] = "Scientific_Binomial"
    names(newSpeciesData)[newCommonNameColumn] = "CommonName"
    names(newSpeciesData)[newScientificNameColumn] = "ScientificNameFull"
    newSpeciesDataset <<- newSpeciesData[,c(1,newCommonNameColumn,newScientificNameColumn)]
  }
  
  newData
}


# -- Function to add new species to main data --

addMainSpecies = function(speciesToAdd = newSpeciesDataset, mergedData = combinedData){
  newSpeciesDataLong = as.data.frame(matrix(NA, ncol = ncol(mergedData), nrow = nrow(speciesToAdd)))
  newSpeciesDataLong[,1] = speciesToAdd[,1]
  newSpeciesDataLong[,2] = speciesToAdd[,2]
  newSpeciesDataLong[,3] = speciesToAdd[,3]
  names(newSpeciesDataLong) = names(mergedData)
  
  mergedData = rbind(mergedData, newSpeciesDataLong)
  
  mergedData
}

# -- Adding column function --
#This function adds a column with the correct order into the merged data.
#It supports renaming the column, and supports adding "name" columns before the divider. 
addColumn = function(dataSet, column, columnRename = NA, mergedData = combinedData, nameColumn = F){
  inColumnIndex = which(names(dataSet) == column)
  columnToAdd = dataSet[[inColumnIndex]]
  inColumnName = names(dataSet)[inColumnIndex]
  if(!is.na(columnRename)){
    outColumnName = columnRename
  }else{
    outColumnName = inColumnName
  }
  
  if(nameColumn){
    mergedData = add_column(mergedData, debugInsertedColumn = columnToAdd, .before = "Divider")
  }else{
    mergedData = add_column(mergedData, debugInsertedColumn = columnToAdd)
  }
  
  outColumnIndex = which(names(mergedData) == "debugInsertedColumn")
  names(mergedData)[outColumnIndex] = outColumnName
  
  mergedData
}
