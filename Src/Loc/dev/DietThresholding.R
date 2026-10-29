#---------------------------------------------------------------------
# --- Futzing with classification --- 
# --------------------------------------------------------------------

mergedData = read.csv("Data/mergedData.csv")

mergedData = mergedData[,-c(3,4,6,7,9,10,11,12,13,14,15,16,17,19,20,21,23)]
mergedData = mergedData[,c(1:26)]
#Set up columsn which are combinations of other diet info 
mergedData$Diet.VertAll = mergedData$Diet.Vend + mergedData$Diet.Vect + mergedData$Diet.Vfish+mergedData$Diet.Vunk
mergedData$VertAllPlusScav = mergedData$Diet.Scav + mergedData$Diet.VertAll
mergedData = mergedData[,c(1:25, 27, 26)]
mergedData = mergedData[which(mergedData$ZoonomiaTip != "vs_NA"),]
centralData = mergedData

#These are a few variables used in functions but created in the global environment as to be modified centrally. 
intermediatePredtoryDietLabels = c("EqualSubsetsPredator", "EqualPredatorPiscivore", "OtherPredator")
allPredatorDiets = append(intermediatePredtoryDietLabels, c("Invertivore", "Vertivore", "zMixedPredator"))
pantheriaColName = "panTheriaTrophicLevelCharacter"
walkerColName = "Meyer.Lab.Classification.Compressed"


#---------------------------------------------------------------------
# --- Categorization functions--- 


#This function categorizes diet in two steps, the three-diet step and the four-diet step. 
#The three diet step ensures that all species with prey diets > the main threshold are in one of: 
  #Vertivore, Invertivore, PredatorHandlingDiet 
#These categories (and any derivatives) can then all be combined back into a single category for a three diet analysis
#The four diet analysis then subsets PredatorHandlingDiet into more diet labels, which can be handled in different ways by later steps


categorizeDiet = function(mainThreshold = 90, secondaryThreshold = 50, centralData = mergedData){
  
  newCategorization = rep(NA)
  
  newCategorization[which(is.na(mergedData$Diet.Inv))] = "NoData" # Set all of the ones with an actual NA value to this placeholder
  
  #Handle the direct conversions based on the primary threshold
  newCategorization[which(mergedData$Diet.PlantAll >= mainThreshold)] = "Herbivore"
  newCategorization[which((mergedData$VertAllPlusScav + mergedData$Diet.Inv) >= mainThreshold)] = "PredatorHandlingDiet"
  newCategorization[which(mergedData$Diet.Inv >= mainThreshold)] = "Invertivore"
  newCategorization[which(mergedData$VertAllPlusScav >= mainThreshold)] = "Vertivore"
  newCategorization[is.na(newCategorization)] = "Omnivore" #Anything that still has an NA (because the NoDatas have been moved) is one that does not fit the above filters
  
  #Subset mixedPredators 
  if(secondaryThreshold == "LPD" | secondaryThreshold == "LargerPredatoryDiet"){
    newCategorization[which( newCategorization == "PredatorHandlingDiet" & mergedData$Diet.Inv > mergedData$VertAllPlusScav)] = "Invertivore"
    newCategorization[which( newCategorization == "PredatorHandlingDiet" & mergedData$VertAllPlusScav > mergedData$Diet.Inv)] = "Vertivore"
  }else{
    secondaryThreshold = as.integer(secondaryThreshold)
    newCategorization[which( newCategorization == "PredatorHandlingDiet" & mergedData$Diet.Inv > secondaryThreshold)] = "Invertivore"
    newCategorization[which( newCategorization == "PredatorHandlingDiet" & mergedData$VertAllPlusScav > secondaryThreshold)] = "Vertivore"
    #Get a category for when vertibrates and invertebrates are exactly equal 
    newCategorization[which( newCategorization == "PredatorHandlingDiet" & mergedData$VertAllPlusScav == mergedData$Diet.Inv)] = "EqualSubsetsPredator"
    
    #Make a category for when the balanced predators are half fish, so that they can optionally be handled differently 
    newCategorization[which( newCategorization == "EqualSubsetsPredator" & mergedData$Diet.Vfish == mergedData$Diet.Inv)] = "EqualPredatorPiscivore"
  
    #Make a catch-all category for all remaining combined predators
    newCategorization[which( newCategorization == "PredatorHandlingDiet")] = "OtherPredator"
  }
  
  newCategorization[which(newCategorization == "NoData")] = NA  #Swap the no datas back to actual NAs. 
  
  newCategorization
}




makeCategorySetForThreshold = function(mainThreshold, secondaryThresholdSet = c(50, 70, mainThreshold-20), centralData = mergedData){
  
  newData = centralData[,c(1,2)] #make a demo dataframe of the same length of the original main data 
  
  #secondaryThresholdSet = append(mainThreshold, secondaryThresholdSet)
  secondaryThresholdSet = unique(secondaryThresholdSet) #Remove any duplicates (if the mainThreshold-20 is equal to one of the presets)
  
  newData$StrictRaw = categorizeDiet(mainThreshold, secondaryThreshold = mainThreshold, centralData)
  newData$ThreeCategory = newData$StrictRaw
  newData$ThreeCategory[which(newData$ThreeCategory %in% allPredatorDiets)] = "CombinedPredator"
  newData$Strict = newData$StrictRaw
  newData$Strict[which(newData$Strict %in% intermediatePredtoryDietLabels)] = "zMixedPredator"
 
  
  #First, make a category with the raw categorizations
  for(i in secondaryThresholdSet){
    if(i < mainThreshold | i=="LPD"){
    currentData = newData[,c(1,2)] #make a dataframe of the same size
    currentData$RawCategorizations = categorizeDiet(mainThreshold, secondaryThreshold = i, centralData)
    
    currentData$Temp = currentData$RawCategorizations # make a column that has the shared conversions for EqualPiscToVert and EqualPiscToMixed
    currentData$Temp[which(currentData$Temp == "EqualSubsetsPredator")] = "zMixedPredator"
    currentData$Temp[which(currentData$Temp == "OtherPredator")] = "zMixedPredator"
  
    if(i <= 50){
      currentData$SubthresholdEqualPiscToVert = currentData$Temp
      currentData$SubthresholdEqualPiscToVert[which(currentData$SubthresholdEqualPiscToVert == "EqualPredatorPiscivore")] = "Vertivore"
      
      currentData$SubthresholdEqualPiscToMixed = currentData$Temp
      currentData$SubthresholdEqualPiscToMixed[which(currentData$SubthresholdEqualPiscToMixed == "EqualPredatorPiscivore")] = "zMixedPredator"  
    }else{
      currentData$Subthreshold = currentData$Temp
      currentData$Subthreshold[which(currentData$Subthreshold == "EqualPredatorPiscivore")] = "Vertivore"
    }
    
    

  
    currentData = currentData[,-c(1,2,4)]
    
    names(currentData) = paste0(names(currentData), i)
    newData = cbind(newData, currentData)
    }
    
    
  }
  
  newData = newData[,-c(1,2)]
  names(newData) = paste0("a",mainThreshold,names(newData))
  
  newData
  
}




#---------------------------------------------------------------------
# --- Run Categorizations --- 

thresholdSet = c(90, 80, 70)

thresholdSet = c(90)
thresholdSet = c(50, 90)

CompareDiets = mergedData
for(i in thresholdSet){
  thresholdData = makeCategorySetForThreshold(i)
  CompareDiets = cbind(CompareDiets, thresholdData)
}

CompareDiets = mergedData
for(i in thresholdSet){
  thresholdData = makeCategorySetForThreshold(i, i)
  CompareDiets = cbind(CompareDiets, thresholdData)
}
CompareDiets = CompareDiets[,-grep("Strict", names(CompareDiets))]

#---------------------------------------------------------------------
# --- Check against external Data  --- 
# --------------------------------------------------------------------


#These functions convert the various phrasing of diets into diet a three-diet 
# or a four diet classification for comparision against external data 

convertToSymbol3category = function(vector){
  vector = trimws(vector)
  vector = gsub("Carnivore", "C", vector)
  vector = gsub("Piscivore", "C", vector)
  vector = gsub("Hematophagy","C", vector)
  vector = gsub("Vertivore", "C", vector)
  vector = gsub("CombinedPredator", "C", vector)
  
  vector = gsub("Planktivore", "C", vector)
  vector = gsub("Insectivore", "C", vector)
  vector = gsub("Invertivore", "C", vector)
  
  vector = gsub("zMixedPredator", "C", vector)
  
  vector = gsub("_Omnivore", "O", vector)
  vector = gsub("Omnivore", "O", vector)
  vector = gsub("Ambiguous", "O", vector)
  
  
  
  vector = gsub("Herbivore", "H", vector)
  
}

convertToSymbol4category = function(vector){
  vector = trimws(vector)
  vector = gsub("Carnivore", "V", vector)
  vector = gsub("Piscivore", "V", vector)
  vector = gsub("Hematophagy", "V", vector)
  vector = gsub("Vertivore", "V", vector)
  
  vector = gsub("Planktivore", "I", vector)
  vector = gsub("Insectivore", "I", vector)
  vector = gsub("Invertivore", "I", vector)
  
  
  vector = gsub("_Omnivore", "O", vector)
  vector = gsub("Omnivore", "O", vector)
  vector = gsub("Ambiguous", "O", vector)
  
  vector = gsub("zMixedPredator", NA, vector)
  
  vector = gsub("Herbivore", "H", vector)
  
}

#these functions run comparisons against the external data, using the global-environment name target column
pantheriaCompare3 = function(vector, centralData = mergedData){
  pantheriaCol = which(names(centralData) == pantheriaColName)
  message(pantheriaCol)
  matchValue = (convertToSymbol3category(vector) == convertToSymbol3category(centralData[,pantheriaCol]))
  matchValue
}
pantheriaCompare4 = function(vector, centralData = mergedData){
  pantheriaCol = which(names(centralData) == pantheriaColName)
  matchValue = (convertToSymbol4category(vector) == convertToSymbol4category(centralData[,pantheriaCol]))
  matchValue
}

walkerCompare3 = function(vector, centralData = mergedData){
  walkerCol = which(names(centralData) == walkerColName)
  matchValue = (convertToSymbol3category(vector) == convertToSymbol3category(centralData[,walkerCol]))
  matchValue
}
walkerCompare4 = function(vector, centralData = mergedData){
  walkerCol = which(names(centralData) == walkerColName)
  matchValue = (convertToSymbol4category(vector) == convertToSymbol4category(centralData[,walkerCol]))
  matchValue
}


colsToCompare = which(!colnames(CompareDiets) %in% colnames(mergedData))
for(i in colsToCompare){
  CompareDiets$newCol = walkerCompare4(CompareDiets[,i])
  names(CompareDiets)[length(names(CompareDiets))] = paste0(names(CompareDiets)[i], "VsWalker4")
  
}
for(i in colsToCompare){
  CompareDiets$newCol = pantheriaCompare3(CompareDiets[,i])
  names(CompareDiets)[length(names(CompareDiets))] = paste0(names(CompareDiets)[i], "VsPantheria3")
}
for(i in colsToCompare){
  CompareDiets$newCol = walkerCompare3(CompareDiets[,i])
  names(CompareDiets)[length(names(CompareDiets))] = paste0(names(CompareDiets)[i], "VsWalker3")
}
for(i in colsToCompare){
  CompareDiets$newCol = pantheriaCompare4(CompareDiets[,i])
  names(CompareDiets)[length(names(CompareDiets))] = paste0(names(CompareDiets)[i], "VsPantheria4")
}



library(tidyverse)


makeBarPlot = function(colNameSuffix){
  plotData <- CompareDiets %>%
    select(contains(colNameSuffix)) %>%
    select(!contains("Raw")) %>%
    pivot_longer(
      cols = everything(),
      names_to = "Column",
      values_to = "Value"
    ) %>%
    filter(!is.na(Value)) %>%
    mutate(Column = str_remove(Column, colNameSuffix))
  
  percentData <- plotData %>%
    group_by(Column) %>%
    summarise(
      TruePercent = mean(Value) * 100,
      Total = n(),
      .groups = "drop"
    )
  
  sizeData <- CompareDiets %>%
    select(contains(percentData[[1]])) %>%
    select(!contains("Raw")) %>%
    select(!contains("Vs")) %>%
    pivot_longer(
      cols = everything(),
      names_to = "Column",
      values_to = "Value"
    )%>%
    
    group_by(Column) %>%
    summarise(
      
      TruePercent = paste0(substring(names(table(Value)),1, 1), "=", table(Value)),
      Total = n(),
      .groups = "drop"
    )  %>%
    group_by(Column) %>%
    summarise(
      TruePercent = paste(TruePercent, collapse = "\n"),
      Total = first(Total),
      .groups = "drop"
    )
  
  plot = ggplot(plotData, aes(x = Column, fill = Value)) +
    geom_bar() +
    geom_text(
      data = percentData,
      aes(
        x = Column,
        y = Total,
        label = sprintf("%.1f%%", TruePercent)
      ),
      inherit.aes = FALSE,
      vjust = -0.5
    ) +
    geom_text(
      data = sizeData,
      aes(
        x = Column,
        y = 50,
        label = TruePercent
      ),
      inherit.aes = FALSE,
      vjust = -0.5
    ) +
    labs(
      x = "Column",
      y = "Count",
      fill = "Value",
      title = colNameSuffix
    ) +
    scale_fill_manual(
      values = c("TRUE" = "steelblue", "FALSE" = "tomato")
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.12))
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  plot
}


Walker3PLot = makeBarPlot("VsWalker3")
Walker4PLot = makeBarPlot("VsWalker4")
Pantheria3PLot = makeBarPlot("vsPantheria3")
Pantheria4PLot = makeBarPlot("vsPantheria4")

Walker4PLot

table(CompareDiets$a90Subthreshold70) / table(CompareDiets$a90Strict)
table(CompareDiets$a90SubthresholdEqualPiscToVert50) / table(CompareDiets$a90Strict)[1:4]
table(CompareDiets$a90SubthresholdEqualPiscToMixed50) / table(CompareDiets$a90Strict)
