clusterRun = F
clusterRun = T
library(RERconverge)
library(ggvenn)
library(stats)
library(combinat)
library(ggpointdensity)
library(viridis)
library(ggpmisc)
library(eulerr)
library(gridExtra)
library(grid)
library(gridGraphics)
source("Src/Reu/cmdArgImport.R")

# -- argument setup  -- 
significanceCutoff = 0.05
prefix = "CategoricalInsvertivoreTreeLiamInference"
pairwiseSets = c("Herbivore-Insectivore", "Herbivore-Vertivore", "Carnivore-Herbivore", "Herbivore-Omnivore", "Insectivore-Vertivore", "Omnivore-Vertivore", "Invertivore-Omnivore")
geneSet = "KeggReactome"
vennDiagramSet = c("Herbivore-Invertivore", "Herbivore-Vertivore", "Carnivore-Herbivore")  
vennColorset = c("darkblue", "red", "orange")
usingGo = !is.null(geneSet)
saveCombinedData = T
saveCombinedData = F
bothAxis = T
saveData = T
saveData = F



args = c("r=CategoricalInsvertivoreTree")
args = c("r=CategoricalInsvertivoreTreeLiamInference")

# -- Standard Startup code -- 
if(clusterRun)args = commandArgs(trailingOnly = TRUE)
{  # Bracket used for collapsing purposes
  #File Prefix
  if(!is.na(cmdArgImport('r'))){
    filePrefix = cmdArgImport('r')
  }else{
    stop("THIS IS AN ISSUE MESSAGE; SPECIFY FILE PREFIX")
  }
  
  #  Output Directory 
  if(!dir.exists("Output")){                                      #Make output directory if it does not exist
    dir.create("Output")
  }
  outputFolderNameNoSlash = paste("Output/",filePrefix, sep = "") #Set the prefix sub directory
  if(!dir.exists(outputFolderNameNoSlash)){                       #create that directory if it does not exist
    dir.create(outputFolderNameNoSlash)
  }
  outputFolderName = paste("Output/",filePrefix,"/", sep = "")
  
  #  Force update argument
  forceUpdate = FALSE
  if(!is.na(cmdArgImport('v'))){                                 #Import if update being forced with argument 
    forceUpdate = cmdArgImport('v')
    forceUpdate = as.logical(forceUpdate)
  }else{
    message("Force update not specified, not forcing update")
  }
}



#------------------------------------------------
# -- OVerlap Figure exclusive code -- 
#------------------------------------------------


# -- Read Data 
combinedGeneDataFilename = paste0(outputFolderName, filePrefix, "combinedGeneResults.rds")
combinedResults = readRDS(combinedGeneDataFilename)
significanceColumns = names(combinedResults)[grep("significant", names(combinedResults))]
geneSignificanceResults = combinedResults[, names(combinedResults) %in% significanceColumns]


combinedGODataFilename = paste0(outputFolderName, filePrefix, "combinedGOResults-", geneSet, ".rds")
GoCombinedResults = readRDS(combinedGODataFilename)
GoSignificanceColumns = names(GoCombinedResults)[grep("significant", names(GoCombinedResults))]
GoSignificanceResults = GoCombinedResults[, names(GoCombinedResults) %in% GoSignificanceColumns]


# -- make resources to prefix-phenotype conversion 
prefixSet = NULL
prefixList = NULL
for(i in 1:length(pairwiseSets)){
  currentSet = pairwiseSets[i]
  correlationSubsetName = gsub("-", " - ", currentSet)
  correlationPrefix = paste(substr(strsplit(currentSet, split = "-")[[1]],1,1), collapse = '')
  prefixEntry = correlationPrefix; names(prefixEntry) = currentSet; prefixSet = append(prefixSet, prefixEntry)
  for(j in 1:2){
    prefixSingle = strsplit(correlationPrefix, split = "")[[1]][j]; names(prefixSingle) = strsplit(currentSet, split="-")[[1]][j]; prefixList = append(prefixList, prefixSingle)
  }
}


prefixList = unlist(prefixList); prefixList = prefixList[!duplicated(prefixList)]

# Debug code line for personal use 
names(prefixList)[which(names(prefixList) == "Insectivore")] = "Invertivore"

replacePrefixWithName = function(x) {
  inversePrefixList = setNames(names(prefixList), prefixList); parts = unlist(strsplit(x, "-")); fullNames = inversePrefixList[parts]; paste(fullNames, collapse = "-")
}

addDashes = function(vector) {
  sapply(vector, function(s) paste(strsplit(s, "")[[1]], collapse = "-"))
}

# ----- Make Plots ------ 

# -- Make rho value corrleation plots --- 
{
  
corrleationColumnType = "-Rho"
  
  
grep(corrleationColumnType, names(combinedResults))
rhoValues = combinedResults[,grep(corrleationColumnType, names(combinedResults))]

rhoComparisions = names(rhoValues)
rhoPhenotypes = strsplit(gsub(corrleationColumnType, "", rhoComparisions), split = "")
commonBackground = Reduce(intersect, rhoPhenotypes)

if(length(commonBackground) == 1){
  for(i in 1:length(rhoPhenotypes)){ #invert tho if background in second postion so rho has consistent meaning relative to background
    if(rhoPhenotypes[[i]][1] != commonBackground){
      cat("Inverting rho of ", rhoPhenotypes[[i]] , "becuase background is in first position.")
      rhoValues[i] = -1*rhoValues[i]
    }
  }
}

densityScaleSet = NULL

#generate the plots slim-ly to get the desired density scale 
for(i in 1:length(rhoValues)){
  xName = names(rhoValues)[i]
  if(i+1 <= length(rhoValues)){
    for(j in (i+1):length(rhoValues)){
      yName = names(rhoValues)[j]
      
      rhoCorrellPlot = ggplot(rhoValues, aes(x = .data[[xName]], y = .data[[yName]])) + 
        geom_point() + geom_pointdensity() + scale_color_viridis()
      
      
      denstiyScaleValue = ggplot_build(rhoCorrellPlot)$plot$scales$scales[[1]]$get_limits()[2]
      densityScaleSet = append(densityScaleSet, denstiyScaleValue)
      rm(rhoCorrellPlot)
    }
  }
}
densityScale = c(1, max(densityScaleSet))


rhoPlotSet = list()
netIndex= 0
for(i in 1:length(rhoValues)){
  xName = names(rhoValues)[i]
  if(i <= length(rhoValues)){
    if(bothAxis){jStart = 1}else{jStart = i+1}
    for(j in (jStart):length(rhoValues)){
      yName = names(rhoValues)[j]
      yLabel =  paste0(replacePrefixWithName(addDashes(gsub(corrleationColumnType, "", yName))), " Dunn Z Statistic")
      xLabel =  paste0(replacePrefixWithName(addDashes(gsub(corrleationColumnType, "", xName))), " Dunn Z Statistic")
      
      rhoCorrellPlot = ggplot(rhoValues, aes(x = .data[[xName]], y = .data[[yName]])) + 
        geom_point() + geom_pointdensity() + scale_color_viridis(name = "Density of genes", limits = densityScale) + 
        stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),formula = y ~ x,parse = TRUE, size = 6) +
        theme_classic()+
        xlab(xLabel) + ylab(yLabel)+
        theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))
      
      netIndex = netIndex +1
      rhoPlotSet[[netIndex]] = rhoCorrellPlot
      names(rhoPlotSet)[netIndex] = paste(xName, yName, sep="-")
      rm(rhoCorrellPlot)
    }
  }
}
# add a null plot by using two using permulations as a comparision 

#pdf()
print(rhoPlotSet)
#dev.off()
}



# -- Make proportional venn diagram via eulerr ---
{ 
  vennPhenotypes = (strsplit(vennDiagramSet, "-"))
  commonBackground = Reduce(intersect, vennPhenotypes)
  commonBackground = substr(commonBackground, 1,1)
  # Make required functions 
  
  trimSignificanceToVenn = function(significanceResults){
    trimableComparisions = gsub("-significant", "", names(significanceResults))
    trimableComparisions = addDashes(trimableComparisions)
    trimableComparisions = sapply(trimableComparisions, replacePrefixWithName)
    vennSignificanceResults = significanceResults[match(vennDiagramSet, trimableComparisions)]
  }
  
  
  makeVennPlot = function(vennInputDataframe, mainTitle, plot = T){
    # Build logical vectors for each set
    set1 <- vennInputDataframe[1] == TRUE
    set2 <- vennInputDataframe[2] == TRUE
    set3 <- vennInputDataframe[3] == TRUE
    
    # Create the Venn counts for each region
    vennCounts = c(
      "Column One" = sum(set1 & !set2 & !set3, na.rm = T),
      "Column Two" = sum(!set1 & set2 & !set3, na.rm = T),
      "Column Three" = sum(!set1 & !set2 & set3, na.rm = T),
      "Column One&Column Two" = sum(set1 & set2 & !set3, na.rm = T),
      "Column One&Column Three" = sum(set1 & !set2 & set3, na.rm = T),
      "Column Two&Column Three" = sum(!set1 & set2 & set3, na.rm = T),
      "Column One&Column Two&Column Three" = sum(set1 & set2 & set3, na.rm = T)
    )
    
    comparisonPrefixes = gsub("-significant", "", names(vennInputDataframe))
    comparisonPrefixes = gsub(commonBackground, "", comparisonPrefixes)
    comparisonNames = sapply(comparisonPrefixes, replacePrefixWithName)
    names(vennCounts)=c(comparisonNames[1], comparisonNames[2], comparisonNames[3], paste(comparisonNames[1], comparisonNames[2], sep="&"),paste(comparisonNames[1], comparisonNames[3], sep="&"),paste(comparisonNames[2], comparisonNames[3], sep="&"), paste(comparisonNames[1], comparisonNames[2], comparisonNames[3], sep="&"))
    
    # - make venn diagram labels, with the combination section on two lines and the solo sections on one 
    totalVennValues = sum(vennCounts)
    vennLabels = paste0(
      vennCounts, "\n (", round(vennCounts / totalVennValues * 100, 1), "%)"
    )
    for(i in 1:3){
      vennLabels[i] = paste0(
        vennCounts[i], " (", round(vennCounts[i] / totalVennValues * 100, 1), "%)"
      )
    }
    
    # Create Euler diagram
    fit = euler(vennCounts)
    outPlot = plot(fit,
         fills = list(fill = vennColorset, alpha = 0.5),
         labels = list(font = 4),
         quantities = list(labels = vennLabels, font = 3),
         main = mainTitle, theme)
    if(plot){print(outPlot)}
    return(outPlot)
  }
}


if(length(pairwiseSets)==3){ #can simply run directly if only running on three categories. 
  geneVenn = makeVennPlot(geneSignificanceResults, paste0("Genes (p.adj < ", significanceCutoff, ")"))
  if(usingGo){goVenn = makeVennPlot(GoSignificanceResults, paste0("GO Categories (p.adj < ", significanceCutoff, ")"))}
}else if(length(vennDiagramSet)==3){ #if a correct set of three categories has been given for the venn diagram set, narrow down a larger selection to that set, then run the vennDiagram. 
  vennGeneSignificanceResults = trimSignificanceToVenn(geneSignificanceResults)
  geneVenn = makeVennPlot(vennGeneSignificanceResults, paste0("Genes (p.adj < ", significanceCutoff, ")"))
  if(usingGo){
    vennGoSignificanceResults = trimSignificanceToVenn(GoSignificanceResults)
    goVenn = makeVennPlot(vennGoSignificanceResults, paste0("GO Categories (p.adj < ", significanceCutoff, ")"))
  }
  
}


#---
#---
#---
#--


# --- Hypergeometric and chi squared tests -----
overlapResults = list()
totalSamples = nrow(significanceResults)
for(i in 1:length(overlapValues)){
  comboCol = overlapValues[[i]]
  comboName = names(overlapValues)[i]
  
  # Get columns that were used in this combination
  involvedSets = unlist(strsplit(comboName, "-"))
  involvedSets = involvedSets[involvedSets != "Overlap"]
  if(length(involvedSets)==2){
    if(length(which(grepl(paste0("^", involvedSets[1], "-"), names(significanceResults))))> 1){
      sigColumnOne = significanceResults[, grepl(paste0("^", involvedSets[1], "-"), names(significanceResults))][[1]]
    }else{
      sigColumnOne = significanceResults[, grepl(paste0("^", involvedSets[1], "-"), names(significanceResults))]
    }
    
    if(length(which(grepl(paste0("^", involvedSets[2], "-"), names(significanceResults))))> 1){
      sigColumnTwo = significanceResults[, grepl(paste0("^", involvedSets[2], "-"), names(significanceResults))][[1]]
    }else{
      sigColumnTwo = significanceResults[, grepl(paste0("^", involvedSets[2], "-"), names(significanceResults))]
    }
    
    numberSignificantOne = sum(sigColumnOne, na.rm = T)
    numberSignificantTwo = sum(sigColumnTwo, na.rm = T) 
    
    observedOverlap = sum(comboCol, na.rm = T)
    
    sigOneNotOverlap = numberSignificantOne - observedOverlap
    sigTwoNotOverlap = numberSignificantTwo - observedOverlap
    sigNone = totalSamples - numberSignificantOne - numberSignificantTwo + observedOverlap #adding the overlap back accounts for those being in both subtractions
    
    
    
    contingency_matrix <- matrix(c(observedOverlap, sigOneNotOverlap, sigTwoNotOverlap, sigNone), nrow = 2)
    print(contingency_matrix)
    # Fisher's exact test
    fisherData = fisher.test(contingency_matrix, alternative = "greater")
    overlapResults[[i]] = fisherData
    names(overlapResults)[i] = comboName
  }else{
    if(length(which(grepl(paste0("^", involvedSets[1], "-"), names(significanceResults))))> 1){
      sigColumnOne = significanceResults[, grepl(paste0("^", involvedSets[1], "-"), names(significanceResults))][[1]]
    }else{
      sigColumnOne = significanceResults[, grepl(paste0("^", involvedSets[1], "-"), names(significanceResults))]
    }
    
    if(length(which(grepl(paste0("^", involvedSets[2], "-"), names(significanceResults))))> 1){
      sigColumnTwo = significanceResults[, grepl(paste0("^", involvedSets[2], "-"), names(significanceResults))][[1]]
    }else{
      sigColumnTwo = significanceResults[, grepl(paste0("^", involvedSets[2], "-"), names(significanceResults))]
    }
    
    if(length(which(grepl(paste0("^", involvedSets[3], "-"), names(significanceResults))))> 1){
      sigColumnThree = significanceResults[, grepl(paste0("^", involvedSets[3], "-"), names(significanceResults))][[1]]
    }else{
      sigColumnThree = significanceResults[, grepl(paste0("^", involvedSets[3], "-"), names(significanceResults))]
    }
    
    df <- data.frame(sigColumnOne, sigColumnTwo, sigColumnThree)
    
    
    # Tabulate all combinations (TRUE/FALSE for each of the 3 sets)
    three_way_table <- table(df$sigColumnOne, df$sigColumnTwo, df$sigColumnThree)
    
    
    contingencyDf <- as.data.frame(three_way_table)
    names(contingencyDf) <- c(involvedSets[1], involvedSets[2], involvedSets[3], "Count")
    
    loglin_result <- loglin(three_way_table, margin = list(c(1), c(2), c(3)), fit = TRUE)
    
    # Print test statistics
    out1 = paste("Likelihood Ratio Statistic (G²):", loglin_result$lrt, "\n")
    out2 = paste("Degrees of Freedom:", loglin_result$df, "\n")
    out3 = paste("p-value:", pchisq(loglin_result$lrt, df = loglin_result$df, lower.tail = FALSE), "\n")
    
    logLInSummary = paste(out1, out2, out3)
    cat(logLInSummary)
    
    overlapResults[[i]] = logLInSummary
    names(overlapResults)[i] = comboName
  }
  
}

overlapResults









#------------------------------------------------
# -- Make combined data - now in other script -- 
#------------------------------------------------
{
  # -- Make central RERData object 
  
  getComparisionDifference = function(dataframe, colOne, colTwo){
    colOneIndex = names(dataframe)[which(names(dataframe) == colOne)]
    colTwoIndex = names(dataframe)[which(names(dataframe) == colTwo)]
    distanceFromEqual = abs(dataframe[colOneIndex] - dataframe[colTwoIndex]) / sqrt(2)
    distanceFromEqual
  }
  
  pairwiseCorrelationFileName = paste(outputFolderName, filePrefix, "PairwiseCorrelationFile.rds", sep= "") #make a name for the pairwise comparisons based on prefix
  correlationResults = readRDS(pairwiseCorrelationFileName)
  
  {
    combinedResults = NA
    combinedDrivers = NA
    combinedBinaries = NA
    prefixSet = NULL
    prefixList = NULL
    
    { # load in the data 
      for(i in 1:length(pairwiseSets)){
        currentSet = pairwiseSets[i]
        correlationSubsetName = gsub("-", " - ", currentSet)
        correlationPrefix = paste(substr(strsplit(currentSet, split = "-")[[1]],1,1), collapse = '')
        prefixEntry = correlationPrefix; names(prefixEntry) = currentSet; prefixSet = append(prefixSet, prefixEntry)
        for(j in 1:2){
          prefixSingle = strsplit(correlationPrefix, split = "")[[1]][j]; names(prefixSingle) = strsplit(currentSet, split="-")[[1]][j]; prefixList = append(prefixList, prefixSingle)
        }
        currentDataframe = which(names(correlationResults) == correlationSubsetName) 
        
        currentResults = correlationResults[[currentDataframe]]
        currentResults$significant = currentResults$p.adj < significanceCutoff
        
        
        names(currentResults) = paste0(correlationPrefix, "-", names(currentResults))
        combinedResults = cbind(combinedResults, currentResults)
        
        driverFilename = paste0(outputFolderName, currentSet, "/", filePrefix, currentSet, "DriverTable.rds")
        if(file.exists(driverFilename)){
          driverTable = readRDS(driverFilename)
          driverTable = driverTable[,-grep("main", names(driverTable))]
          names(driverTable)[which(names(driverTable) == "Driver")] = paste0(correlationPrefix, "-", "Driver")
          names(driverTable)[which(names(driverTable) == "DriverNumeric")] = paste0(correlationPrefix, "-", "DriverNumeric")    
          
          combinedDrivers = cbind(combinedDrivers, driverTable[,grep(paste0(correlationPrefix,"-"), names(driverTable))])
          combinedBinaries = cbind(combinedBinaries, driverTable[,-grep(paste0(correlationPrefix,"-"), names(driverTable))])
          rm(driverTable)
        }
        rm(currentResults)
        
      }
    }
    
    combinedResults = combinedResults[,-1]
    combinedDrivers = combinedDrivers[,-1]
    combinedBinaries = combinedBinaries[,-1]
    combinedBinaries = combinedBinaries[,-grep(".1", names(combinedBinaries))]
    
    combinedResults = cbind(combinedResults, combinedDrivers)
    combinedResults = cbind(combinedResults, combinedBinaries)
    rm(combinedDrivers); rm(combinedBinaries)
    
    
    
    # -- Add overlap information -- 
    significanceColumns = names(combinedResults)[grep("significant", names(combinedResults))]
    geneSignificanceResults = combinedResults[, names(combinedResults) %in% significanceColumns]
    
    for(i in 2:length(significanceColumns)){
      combinations = combn(significanceColumns, i, simplify = FALSE)
      for(j in 1:length(combinations)){
        currentCombination = combinations[[j]]
        headers = gsub("-.*","",  currentCombination)
        comboName = paste0(paste0(headers, collapse = "-"), "-Overlap")
        
        colsToCompare = geneSignificanceResults[,names(geneSignificanceResults) %in% currentCombination]
        
        comboValue = apply(colsToCompare, 1, function(row) all(row == TRUE) == 1)
        which(comboValue)
        combinedResults$newOverlapColumn = comboValue
        names(combinedResults)[length(names(combinedResults))] = comboName
        rm(colsToCompare)
      }
      rm(combinations)
    }
    
    # -- Add Delta information -- 
    rhoColumns = names(combinedResults)[grep(corrleationColumnType, names(combinedResults))]
    geneRhoColumns = combinedResults[, names(combinedResults) %in% rhoColumns]
    for(i in 2:length(rhoColumns)){
      combinations = combn(rhoColumns, i, simplify = FALSE)
      for(j in 1:length(combinations)){
        currentCombination = combinations[[j]]
        if(length(currentCombination) > 2){next}
        headers = gsub("-.*","",  currentCombination)
        comboName = paste0(paste0(headers, collapse = "-"), "-Delta")
        
        colsToCompare = geneRhoColumns[,names(geneRhoColumns) %in% currentCombination]
        
        deltaValue = abs(colsToCompare[1] - colsToCompare[2]) / sqrt(2)
        
        
        
        combinedResults$newDeltaColumn = deltaValue[,1]
        names(combinedResults)[length(names(combinedResults))] = comboName
        rm(colsToCompare)
      }
      rm(combinations)
    }
    
    # -- save combination -- 
    if(saveCombinedData){
      combinedDataFilename = paste0(outputFolderName, filePrefix, "combinedGeneResults")
      write.csv(combinedResults, paste0(combinedDataFilename, ".csv"))
      if(saveData){saveRDS(combinedResults, paste0(combinedDataFilename, ".rds"))}
    }
  }
  
  
  # --- Import GO Data --- 
  if(usingGo){
    GOResults = list()
    for(i in 1:length(pairwiseSets)){
      currentSet = pairwiseSets[i]
      correlationSubsetName = gsub("-", " - ", currentSet)
      correlationPrefix = paste(substr(strsplit(currentSet, split = "-")[[1]],1,1), collapse = '')
      
      goFilename = paste0(outputFolderName, currentSet, "/", filePrefix, currentSet,"Enrichment-", geneSet, ".rds")
      currentGoData = readRDS(goFilename)[[1]]
      currentGoData$significant = currentGoData$p.adj < significanceCutoff
      
      
      names(currentGoData) = paste0(correlationPrefix, "-", names(currentGoData))
      
      
      driverFilename = paste0(outputFolderName, currentSet, "/", filePrefix, currentSet, "GoDriverTable-", geneSet, ".rds")
      if(file.exists(driverFilename)){
        driverTable = readRDS(driverFilename)
        if(all(rownames(driverTable) == rownames(currentGoData))){
          currentGoData$Driver = driverTable[which(names(driverTable) == "Driver")][[1]]
          currentGoData$DriverNumeric = driverTable[which(names(driverTable) == "DriverNumeric")][[1]]
          
          names(currentGoData)[which(names(currentGoData) == "Driver")] = paste0(correlationPrefix, "-", "Driver")
          names(currentGoData)[which(names(currentGoData) == "DriverNumeric")] = paste0(correlationPrefix, "-", "DriverNumeric")    
          
        }
        rm(driverTable)
      }
      GOResults[[i]] = currentGoData
      names(GOResults)[i] = correlationPrefix
      
    }
    rm(currentGoData)
    
    GoCombinedResults = NA
    for(i in 1:length(GOResults)){
      if(all(rownames(GOResults[[1]]) == rownames(GOResults[[i]]))){
        cat("Combining GO Data", i, "\n")
        GoCombinedResults = cbind(GoCombinedResults, GOResults[[i]])
      }else{
        stop("Rownames of GO results are not the same, there is an issue with the GO data.")
      }
    }
    GoCombinedResults = GoCombinedResults[,-1]
    rm(GOResults)
    
    # -- Add overlap information -- 
    GoSignificanceColumns = names(GoCombinedResults)[grep("significant", names(GoCombinedResults))]
    GoSignificanceResults = GoCombinedResults[, names(GoCombinedResults) %in% GoSignificanceColumns]
    
    for(i in 2:length(GoSignificanceColumns)){
      combinations = combn(GoSignificanceColumns, i, simplify = FALSE)
      for(j in 1:length(combinations)){
        currentCombination = combinations[[j]]
        headers = gsub("-.*","",  currentCombination)
        comboName = paste0(paste0(headers, collapse = "-"), "-Overlap")
        
        colsToCompare = GoSignificanceResults[,names(GoSignificanceResults) %in% currentCombination]
        
        comboValue = apply(colsToCompare, 1, function(row) all(row == TRUE) == 1)
        which(comboValue)
        GoCombinedResults$newOverlapColumn = comboValue
        names(GoCombinedResults)[length(names(GoCombinedResults))] = comboName
        rm(colsToCompare)
      }
      rm(combinations)
    }
    
    if(saveCombinedData){
      combinedGODataFilename = paste0(outputFolderName, filePrefix, "combinedGOResults-", geneSet)
      write.csv(GoCombinedResults, paste0(combinedGODataFilename, ".csv"))
      if(saveData){saveRDS(GoCombinedResults, paste0(combinedGODataFilename, ".rds"))}
    }
  }
}



  




















# ---- Old code ------ 




library("ggplotify")
ggVennGene = as.grob(geneVenn)
ggVennGo = as.grob(goVenn)



names(rhoPlotSet)
c(2,4,1,5) # single axis
rhoPlotsDisplay = grid.arrange(rhoPlotSet[[2]], rhoPlotSet[[4]],rhoPlotSet[[1]],rhoPlotSet[[5]], ncol = 2, padding = 4)
c(3,7,5,8) #both axis
rhoPlotsDisplay = grid.arrange(rhoPlotSet[[3]], rhoPlotSet[[7]],rhoPlotSet[[5]],rhoPlotSet[[8]], ncol = 2, padding = 4)

vennDiagramDisplay = grid.arrange(ggVennGene, ggVennGo, ncol = 2)


rhoPlotSet[2]

rhoPlotName = paste0(outputFolderName, filePrefix, "StatCorrelationPlots.pdf")
vennPlotName = paste0(outputFolderName, filePrefix, "VennPlots.pdf")
overlapPlotName = paste0(outputFolderName, filePrefix, "OverlapPlots.pdf")


pdf(rhoPlotName, 20,10)
grid.arrange(rhoPlotsDisplay)
dev.off()

pdf(vennPlotName, 30,15)
grid.arrange(vennDiagramDisplay)
dev.off()


pdf(overlapPlotName, 24, 24)
grid.arrange(vennDiagramDisplay, rhoPlotsDisplay, ncol = 1, padding = 4)
dev.off()





# --- make edited version for confrence talk ----


rhoPlotSet = list()
netIndex= 0
for(i in 1:length(rhoValues)){
  xName = names(rhoValues)[i]
  if(i <= length(rhoValues)){
    if(bothAxis){jStart = 1}else{jStart = i+1}
    for(j in (jStart):length(rhoValues)){
      yName = names(rhoValues)[j]
      yLabel =  paste0(replacePrefixWithName(addDashes(gsub(corrleationColumnType, "", yName))), "")
      xLabel =  paste0(replacePrefixWithName(addDashes(gsub(corrleationColumnType, "", xName))), "")
      
      rhoCorrellPlot = ggplot(rhoValues, aes(x = .data[[xName]], y = .data[[yName]])) + 
        geom_point() + geom_pointdensity() + scale_color_viridis(name = "Density of genes", limits = densityScale) + 
        stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")),formula = y ~ x,parse = TRUE, size = 10) +
        theme_classic()+
        xlab(xLabel) + ylab(yLabel)+
        theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))+ 
        theme(legend.position="none")
      
      netIndex = netIndex +1
      rhoPlotSet[[netIndex]] = rhoCorrellPlot
      names(rhoPlotSet)[netIndex] = paste(xName, yName, sep="-")
      rm(rhoCorrellPlot)
    }
  }
}


names(rhoPlotSet)
figurePlotName = "Figureplot.png"
png(figurePlotName)
rhoPlotSet[2]
dev.off()


png(figurePlotName)
rhoPlotSet[8]
dev.off()

png(figurePlotName)
rhoPlotSet[3]
dev.off()

png(figurePlotName)
rhoPlotSet[7]
dev.off()


#----
?pdf


which(GoCombinedResults$`HI-significant`)

goCategory = 26

targetGenes = gsub(":.*", "", strsplit(GoCombinedResults$`HI-gene.vals`[goCategory], split = ", ")[[1]])
GoCombinedResults$`HI-stat`[goCategory]
combinedResults$`HI-Rho`[rownames(combinedResults) %in% targetGenes]


# --- block 1 --- 
RERResults = readRDS("Output/CategoricalInsvertivoreTree/CategoricalInsVertivoreTreePairwiseCorrelationFile.rds")

HICorrelations = RERResults$`Herbivore - Insectivore`
HVCorrelations = RERResults$`Herbivore - Vertivore`
HOCorrelations = RERResults$`Herbivore - Omnivore`


HISignificantGenes = rownames(HICorrelations)[which(HICorrelations$p.adj < 0.05)]
HVSignificantGenes = rownames(HVCorrelations)[which(HVCorrelations$p.adj < 0.05)]
HOSignificantGenes = rownames(HOCorrelations)[which(HOCorrelations$p.adj < 0.05)]

sharedGenes = HISignificantGenes[which(HISignificantGenes %in% HVSignificantGenes)]

triSharedGenes = HOSignificantGenes[which(HOSignificantGenes %in% sharedGenes)]


HIDrivingData = readRDS("Output/CategoricalInsvertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreDirectionalityTable.rds")
HIDrivingData[which(rownames(HIDrivingData) %in% sharedGenes),]


HVDrivingData = read.csv("Output/CategoricalInsvertivoreTree/Herbivore-Vertivore/CategoricalInsVertivoreTreeHerbivore-VertivoreDirectionalityTable.csv")
rownames(HVDrivingData) = HVDrivingData$X
HVDrivingData[which(rownames(HVDrivingData) %in% sharedGenes),]

CombinedDrivingData = HIDrivingData
colnames(CombinedDrivingData)[10] = "HIDirectionality"
colnames(CombinedDrivingData)[11] = "HIDirectionalityNumeric"
CombinedDrivingData = cbind(CombinedDrivingData, HVDrivingData[,c(11,12)])
colnames(CombinedDrivingData)[12] = "HVDirectionality"
colnames(CombinedDrivingData)[13] = "HVDirectionalityNumeric"

all.equal(rownames(HVDrivingData), rownames(HIDrivingData))


SharedCombinedDrivingData = CombinedDrivingData[which(rownames(CombinedDrivingData) %in% sharedGenes),]

nrow(SharedCombinedDrivingData[which(SharedCombinedDrivingData$HIDirectionality == "Herbivore" & SharedCombinedDrivingData$HVDirectionality =="Herbivore"),])
nrow(SharedCombinedDrivingData[which(SharedCombinedDrivingData$HIDirectionality == "Herbivore" | SharedCombinedDrivingData$HVDirectionality =="Herbivore"),])


HISharedDriving = HIDrivingData[which(rownames(HIDrivingData) %in% sharedGenes),]
HVSharedDriving = HVDrivingData[which(rownames(HVDrivingData) %in% sharedGenes),]

rownames(HISharedDriving[which(HISharedDriving$directionality == "Herbivore"),])
rownames(HVSharedDriving[which(HVSharedDriving$directionality == "Herbivore"),])



length(which(rownames(HISharedDriving[which(HISharedDriving$directionality == "Herbivore"),]) %in% rownames(HVSharedDriving[which(HVSharedDriving$directionality == "Herbivore"),])))




# ----------- Group of code 2 

HiGeneData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreCorrelationFile.rds")
HvGeneData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Vertivore/CategoricalInsVertivoreTreeHerbivore-VertivoreCorrelationFile.rds")
HoGeneData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Omnivore/CategoricalInsVertivoreTreeHerbivore-OmnivoreCorrelationFile.rds")
IoGeneData = readRDS("Output/CategoricalInsVertivoreTree/Insectivore-Omnivore/CategoricalInsVertivoreTreeInsectivore-OmnivoreCorrelationFile.rds")
IvGeneData = readRDS("Output/CategoricalInsVertivoreTree/Insectivore-Vertivore/CategoricalInsVertivoreTreeInsectivore-VertivoreCorrelationFile.rds")
OvGeneData = readRDS("Output/CategoricalInsVertivoreTree/Omnivore-Vertivore/CategoricalInsVertivoreTreeOmnivore-VertivoreCorrelationFile.rds")
DhiGeneData = readRDS("Output/CategoricalDownsampledInsvertTree/Herbivore-Insectivore/CategoricalDownsampledInsvertTreeHerbivore-InsectivoreCorrelationFile.rds")

arrangeData = function(data, prefix){
  data$index = 1:nrow(data)
  data= data[order(data$p.adj),]
  data$rank = 1:nrow(data)
  data= data[order(data$index),]
  data$index=NULL
  colnames(data) = paste0(prefix, colnames(data))
  return(data)
}

HiGeneData = arrangeData(HiGeneData, "Hi_")
HvGeneData = arrangeData(HvGeneData, "Hv_")
HoGeneData = arrangeData(HoGeneData, "Ho_")
IoGeneData = arrangeData(IoGeneData, "Io_")
IvGeneData = arrangeData(IvGeneData, "Iv_")
OvGeneData = arrangeData(OvGeneData, "Ov_")
DhiGeneData = arrangeData(DhiGeneData, "Dhi_")


combinedData = cbind(HiGeneData, HvGeneData, HoGeneData, IoGeneData, IvGeneData, OvGeneData, DhiGeneData)
combinedData$index = 1:nrow(combinedData)


combinedDataPval = combinedData[,c(2,6,10,14,18,22,26)]
combinedDataPadjVal = combinedData[,c(3,7,11,15,19,23,27)]



length(which(combinedDataPval$Dhi_p.adj < 0.05))

apply(combinedDataPval, MARGIN = 2, mean)

apply(combinedDataPval, 2, function(column) length(which(column < 0.02)))

colnames(combinedDataPval) = c("Herbivore-Insectivore", "Herbivore-Vertivore", "Herbivore-Omnivore", "Insectivore-Omnivore", "Insectivore-Vertivore", "Omnivore-Vertivore", "Downsampled Insectivore-Herbivore")

colnames(combinedDataPval) = c("H-I", "H-V", "H-O", "I-O", "I-V", "O-V", "Downsampled H-I")
combinedDataPval$`Downsampled H-I` = NULL

colnames(combinedDataPadjVal) = c("H-I", "H-V", "H-O", "I-O", "I-V", "O-V", "Downsampled H-I")
combinedDataPadjVal$`Downsampled H-I` = NULL
sigGenesData =data.frame(category = colnames(combinedDataPadjVal), value = apply(combinedDataPadjVal, 2, function(column) length(which(column < 0.02)))) 


library(ggplot2)
library(ggpattern)

ggplot(sigGenesData, aes(x = category, y = value, fill = category, pattern)) +
  geom_bar(stat = "identity", color = "black", show.legend = FALSE) +  # Base bar color
  geom_bar_pattern(
    stat = "identity",
    pattern = "stripe",  # Options: "stripe", "crosshatch", "dots", etc.
    pattern_density = 0.25,
    pattern_fill = c("darkblue", "red", "black","black", "red", "red"),  # Pattern color
    aes(pattern = Category),  # Apply pattern per category
    show.legend = FALSE
  ) +
  theme_minimal() +
  labs(title = "Significant Genes per Pairwise Analysis",
       x = "Category", y = "Number of significant genes") +
  scale_fill_manual(values = c("darkgreen", "darkgreen", "darkgreen", "darkblue", "darkblue", "black"))+
  scale_pattern_fill_manual(values = c("darkblue", "black", "red", "black", "red", "red"))




?barplot
geom_ba
CarnivoreGeneData = readRDS("Output/CategoricalPrunedCarnivoreTree/Carnivore-Herbivore/CategoricalPrunedCarnivoreTreeCarnivore-HerbivoreCorrelationFile.rds")


InsectivoreGeneData$index = 1:nrow(InsectivoreGeneData)
VertivoreGeneData$index = 1:nrow(VertivoreGeneData)
CarnivoreGeneData$index = 1:nrow(CarnivoreGeneData)

InsectivoreGeneData = InsectivoreGeneData[order(InsectivoreGeneData$p.adj),]
VertivoreGeneData = VertivoreGeneData[order(VertivoreGeneData$p.adj),]
CarnivoreGeneData = CarnivoreGeneData[order(CarnivoreGeneData$p.adj),]

InsectivoreGeneData$rank = 1:nrow(InsectivoreGeneData)
VertivoreGeneData$rank = 1:nrow(VertivoreGeneData)
CarnivoreGeneData$rank = 1:nrow(CarnivoreGeneData)

InsectivoreGeneData = InsectivoreGeneData[order(InsectivoreGeneData$index),]
VertivoreGeneData = VertivoreGeneData[order(VertivoreGeneData$index),]
CarnivoreGeneData = CarnivoreGeneData[order(CarnivoreGeneData$index),]

colnames(InsectivoreGeneData) = paste0("I_", colnames(InsectivoreGeneData))
colnames(VertivoreGeneData) = paste0("V_", colnames(VertivoreGeneData))
colnames(CarnivoreGeneData) = paste0("C_", colnames(CarnivoreGeneData))

combinedData = cbind(InsectivoreGeneData, VertivoreGeneData, CarnivoreGeneData)


# ---- part 3


InsectivoreGoData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreEnrichment-GO_Biological_Process_2023.rds")
write.csv(InsectivoreGoData, "Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreEnrichment-GO_Biological_Process_2023.csv")

?correlateWithBinaryPhenotype
install.packages()

#----------------------------------------

InsectivoreGoData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreEnrichment-GO_Biological_Process_2023.rds")
VertivoreGoData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Vertivore/CategoricalInsVertivoreTreeHerbivore-VertivoreEnrichment-GO_Biological_Process_2023.rds")
CarnivoreGoData = readRDS("Output/CategoricalPrunedCarnivoryTree/Carnivore-Herbivore/CategoricalPrunedCarnivoreTreeCarnivore-HerbivoreEnrichment-GO_Biological_Process_2023.rds")


#----------------------------------------
library(ggvenn)
InsectivoreGoData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreEnrichment-GO_Biological_Process_2023.rds")[[1]]
VertivoreGoData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Vertivore/CategoricalInsVertivoreTreeHerbivore-VertivoreEnrichment-GO_Biological_Process_2023.rds")[[1]]
CarnivoreGoData = readRDS("Output/CategoricalPrunedCarnivoreTree/Carnivore-Herbivore/CategoricalPrunedCarnivoreTreeCarnivore-HerbivoreEnrichment-GO_Biological_Process_2023.rds")[[1]]

InsectivoreGoData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreEnrichment-KeggReactome.rds")[[1]]
VertivoreGoData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Vertivore/CategoricalInsVertivoreTreeHerbivore-VertivoreEnrichment-KeggReactome.rds")[[1]]
CarnivoreGoData = readRDS("Output/CategoricalPrunedCarnivoreTree/Carnivore-Herbivore/CategoricalPrunedCarnivoreTreeCarnivore-HerbivoreEnrichment-KeggReactome.rds")[[1]]



InsectivoreGoData = InsectivoreGoData[order(InsectivoreGoData$p.adj),]
VertivoreGoData = VertivoreGoData[order(VertivoreGoData$p.adj),]
CarnivoreGoData = CarnivoreGoData[order(CarnivoreGoData$p.adj),]


signficiantInsectivore = InsectivoreGoData[which(InsectivoreGoData$p.adj <0.05),]
signficiantVertivore = VertivoreGoData[which(VertivoreGoData$p.adj <0.05),]
signficiantCarnivore = CarnivoreGoData[which(CarnivoreGoData$p.adj <0.05),]


valuedInsectivore = InsectivoreGoData[which(InsectivoreGoData$pval <1),]
valuedVertivore = VertivoreGoData[which(VertivoreGoData$pval <1),]
valuedCarnivore = CarnivoreGoData[which(CarnivoreGoData$pval <1),]

vennData = list(
  Insectivore = rownames(valuedInsectivore),
  Vertivore = rownames(valuedVertivore),
  Carnivore = rownames(valuedCarnivore)
)

signficianterInsectivore = InsectivoreGoData[which(InsectivoreGoData$p.adj <0.05),]
signficianterVertivore = VertivoreGoData[which(VertivoreGoData$p.adj <0.05),]
signficianterCarnivore = CarnivoreGoData[which(CarnivoreGoData$p.adj <0.05),]

topInsectivore = InsectivoreGoData[1:100,]
topVertivore = VertivoreGoData[1:100,]
topCarnivore = CarnivoreGoData[1:100,]

vennData = list(
  Insectivore = rownames(signficiantInsectivore),
  Vertivore = rownames(signficiantVertivore),
  Carnivore = rownames(signficiantCarnivore)
)
ggvenn(vennData, fill_color = c("blue", "red", "orange"))


# --------

InsectivoreGeneData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreCorrelationFile.rds")
VertivoreGeneData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Vertivore/CategoricalInsVertivoreTreeHerbivore-VertivoreCorrelationFile.rds")
CarnivoreGeneData = readRDS("Output/CategoricalPrunedCarnivoreTree/Carnivore-Herbivore/CategoricalPrunedCarnivoreTreeCarnivore-HerbivoreCorrelationFile.rds")

InsectivoreGeneData$index = 1:nrow(InsectivoreGeneData)
VertivoreGeneData$index = 1:nrow(VertivoreGeneData)
CarnivoreGeneData$index = 1:nrow(CarnivoreGeneData)

InsectivoreGeneData = InsectivoreGeneData[order(InsectivoreGeneData$p.adj),]
VertivoreGeneData = VertivoreGeneData[order(VertivoreGeneData$p.adj),]
CarnivoreGeneData = CarnivoreGeneData[order(CarnivoreGeneData$p.adj),]

InsectivoreGeneData$rank = 1:nrow(InsectivoreGeneData)
VertivoreGeneData$rank = 1:nrow(VertivoreGeneData)
CarnivoreGeneData$rank = 1:nrow(CarnivoreGeneData)

InsectivoreGeneData = InsectivoreGeneData[order(InsectivoreGeneData$index),]
VertivoreGeneData = VertivoreGeneData[order(VertivoreGeneData$index),]
CarnivoreGeneData = CarnivoreGeneData[order(CarnivoreGeneData$index),]

colnames(InsectivoreGeneData) = paste0("I_", colnames(InsectivoreGeneData))
colnames(VertivoreGeneData) = paste0("V_", colnames(VertivoreGeneData))
colnames(CarnivoreGeneData) = paste0("C_", colnames(CarnivoreGeneData))

combinedData = cbind(InsectivoreGeneData, VertivoreGeneData, CarnivoreGeneData)
combinedData$V_index = NULL
combinedData$C_index = NULL

combinedData = combinedData[order(combinedData$I_p.adj),]


equationLinePlot = function(data, xIn, yIn){
  linearModel = lm(yIn ~ xIn, data = data) 
  
  equation = paste0("y=", round(coef(linearModel)[2], 2), "*x", round(coef(linearModel)[1], 2))
  rSquared = paste("R² =", round(summary(linearModel)$r.squared, 2))
  
  ggplot(data, aes(x = xIn, y = yIn)) + 
    geom_point()+
    geom_smooth(method = "lm")  
}

equationLinePlot(combinedData, "I_Rho", "C_Rho")

linearModel = lm(I_Rho ~ C_Rho, data = combinedData) 

library(gridExtra)



linearModel = lm(-C_Rho ~ I_Rho, data = combinedData) 
equation = paste0("y=", round(coef(linearModel)[2], 2), "*x", round(coef(linearModel)[1], 2))
rSquared = paste("R² =", round(summary(linearModel)$r.squared, 2))
plot1 = ggplot(combinedData, aes(x = I_Rho, y = -C_Rho)) + 
  geom_point()+
  geom_smooth(method = "lm")+
  annotate("text", x = 3, y = 9, label = paste(equation, rSquared, sep = "\n"), color = "blue", size = 10)+
  theme_minimal()

linearModel = lm(-C_Rho ~ V_Rho, data = combinedData) 
equation = paste0("y=", round(coef(linearModel)[2], 2), "*x", round(coef(linearModel)[1], 2))
rSquared = paste("R² =", round(summary(linearModel)$r.squared, 2))
plot2 = ggplot(combinedData, aes(x = V_Rho, y = -C_Rho)) + 
  geom_point()+
  geom_smooth(method = "lm")+
  annotate("text", x = 3, y = 9, label = paste(equation, rSquared, sep = "\n"), color = "blue", size = 10)+
  theme_minimal()


grid.arrange(plot1, plot2, ncol =2)



?lm

library(ggplot2)
ggplot(combinedData, aes(x = I_Rho, y = -C_Rho)) + 
  geom_point()+
  geom_smooth(method = "lm")


plot(combinedData$I_Rho, combinedData$V_Rho)
plot(-combinedData$C_Rho, combinedData$I_Rho)
plot(-combinedData$C_Rho, combinedData$V_Rho)




plot(combinedData$I_rank, combinedData$V_rank, xlim = c(0,3500), ylim = c(0,3500))
plot(combinedData$`I_rank`, combinedData$`C_rank`, xlim = c(0,4000), ylim = c(0,4000))
plot(combinedData$`V_rank`, combinedData$`C_rank`, xlim = c(0,4000), ylim = c(0,4000))
plot(combinedData$V_index, combinedData$V_p.adj)
hist(combinedData$V_p.adj)
length(combinedData$V_p.adj[which(combinedData$V_p.adj < 1)])
length(combinedData$I_p.adj[which(combinedData$I_p.adj < 1)])

?cbind
sigInsectGenes = InsectivoreGeneData[which(InsectivoreGeneData$p.adj <0.05),]
sigVertGenes = VertivoreGeneData[which(VertivoreGeneData$p.adj <0.05),]
sigCarnGenes = CarnivoreGeneData[which(CarnivoreGeneData$p.adj <0.05),]

geneVennData = list(
  Insectivore = rownames(sigInsectGenes),
  Vertivore = rownames(sigVertGenes),
  Carnivore = rownames(sigCarnGenes)
)
ggvenn(geneVennData, fill_color = c("blue", "red", "orange"))
require(venneuler)
v <- venneuler(c(Insectivore=194, Vertivore=145, Carnviore=598, "Insectivore&Vertivore"=0, "Insectivore&Carnviore"=294, "Vertivore&Carnviore"=125, "Carnviore&Vertivore&Insectivore"=75))
plot(v)

valuedInsectGenes = combinedData[which(combinedData$I_p.adj <1),]
valuedVertGenes = combinedData[which(combinedData$V_p.adj <1),]
valuedCarnGenes = combinedData[which(combinedData$C_p.adj <1),]

geneVennData = list(
  Insectivore = rownames(valuedInsectGenes),
  Vertivore = rownames(valuedVertGenes),
  Carnivore = rownames(valuedCarnGenes)
)

ggvenn(geneVennData, fill_color = c("blue", "red", "pink"))