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


palette(c("black", "#514c8e", "#b71568", "#bf7c00"))

vennColorset = c("darkblue", "red", "orange")

geneVennColorset = c("#7570B3", "#E7298A", "#FFA500")

goVennColorset= c("#514c8e", "#b71568", "#bf7c00")







args = c("r=CategoricalInsvertivoreTree")
args = c("r=CategoricalInsvertivoreTreeLiamInference")
args = c("r=CategoricalInsvertivoreTreeFamilyAgnostictLiamInference")
args = c("r=CategoricalInsvertivoreTreeNoYeastLiamInference")



args = c("r=CategoricalInsvertivoreTreeLiamInference", 
         "g=KeggReactome", 
         'n=c("Herbivore-Invertivore", "Herbivore-Vertivore", "Carnivore-Herbivore")',
         "d=T",'l=c("H>P", "P>H")')
args = c("r=CategoricalInsvertivoreTreeLiamInference", 
         "g=GO_Biological_Process_2023", 
         'n=c("Herbivore-Invertivore", "Herbivore-Vertivore", "Carnivore-Herbivore")',
         "d=T",'l=c("H>P", "P>H")')
args = c("r=CategoricalInsvertivoreTreeLiamInference", 
         "g=EnrichmentHsSymbolsFile2", 
         'n=c("Herbivore-Invertivore", "Herbivore-Vertivore", "Carnivore-Herbivore")',
         "d=T",'l=c("H>P", "P>H")')
args = c("r=CategoricalInsvertivoreTreeLiamInference", 
         "g=MGI_Mammalian_Phenotype_Level_4", 
         'n=c("Herbivore-Invertivore", "Herbivore-Vertivore", "Carnivore-Herbivore")',
         "d=T",'l=c("H>P", "P>H")')
args = c("r=CategoricalInsvertivoreTreeLiamInference", 
         "g=DisGeNET", 
         'n=c("Herbivore-Invertivore", "Herbivore-Vertivore", "Carnivore-Herbivore")',
         "d=T",'l=c("H>P", "P>H")')
args = c("r=CategoricalInsvertivoreTreeLiamInference", 
         "g=tissue_specific", 
         'n=c("Herbivore-Invertivore", "Herbivore-Vertivore", "Carnivore-Herbivore")',
         "d=T",'l=c("H>P", "P>H")')

args = c("r=ComplexDietCentralAnalysis", 
         "g=KeggReactome", 
         'n=c("Herbivore-Invertivore", "Herbivore-Vertivore", "Carnivore-Herbivore")',
         "d=T",'l=c("H>P", "P>H")', 'a=T')

args = c("r=ComplexDietCentralAnalysis", 
         'n=c("Herbivore-Invertivore", "Herbivore-Vertivore", "Carnivore-Herbivore")',
         "d=T",'l=c("H>P", "P>H")', 'z=50', 'a=T')


{
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



# --- Argument Imports ---
# Defaults
significanceCutoff = 0.05
geneSet = NULL
usingGo = F
pairwiseSets = NULL
vennDiagramSet = NULL
usePermulations = F
makeDirectional = T
positiveLabel = NULL
negativeLabel = NULL
bothAxis = T
useMainAndRobust = F





{ # Bracket used for collapsing purposes
  
  #significance -- -NOTE THAT THIS IS ONLY FOR THE LABEL, actual numbers are based on the columns
  if(!is.na(cmdArgImport('z'))){
    significanceCutoff = cmdArgImport('z')
    significanceCutoff = as.numeric(significanceCutoff)
  }else{
    message("Significance value not specified, using 0.05")
  }
  
  #geneset
  if(!is.na(cmdArgImport('g'))){
    geneSet = cmdArgImport('g')
  }else{
    message("No geneset specified, running only on genes")
  }
  
  #pairwise Directories
  if(!any(is.na(cmdArgImport('s')))){
    pairwiseSets = cmdArgImport('s')
  }else{
    subdirectories = list.dirs(outputFolderName, recursive = F, full.names = F)
    pairwiseSets = subdirectories[grep("-", subdirectories)]
    message("No pairwise directories specified, using all.")
  }
  
  #Venn Directories
  if(!any(is.na(cmdArgImport('n')))){
    vennDiagramSet = cmdArgImport('n')
  }else{
    message("No Venn directories specified, not suing venn directories.")
  }
  
  #Use Permulations
  if(!any(is.na(cmdArgImport('p')))){
    usePermulations = cmdArgImport('p')
  }else{
    message("No permulations use specified, not using permulations.")
  }
  
  #Make Directional
  if(!any(is.na(cmdArgImport('d')))){
    makeDirectional = cmdArgImport('d')
  }else{
    message("No makeDirectional specified, not making directional.")
  }
  
  #Make Directional
  if(!any(is.na(cmdArgImport('l')))){
    labelSet = cmdArgImport('l')
    positiveLabel = labelSet[1]
    negativeLabel = labelSet[2]
  }else{
    message("No directiona labels specified, not using positive and negative.")
  }
  
  #Use Main robust to alternates
  if(!any(is.na(cmdArgImport('a')))){
    useMainAndRobust = cmdArgImport('a')
    useMainAndRobust = as.logical(useMainAndRobust)
  }else{
    message("Not using main + robustness; either just main or just alternates.")
  }
  

  usingGo = !is.null(geneSet)
  if(is.null(positiveLabel)){
    positiveLabel = "Positive"
  }
  if(is.null(positiveLabel)){
    positiveLabel = "Negative"
  }

}
#------------------------------------------------
# -- OVerlap Figure exclusive code -- 
#------------------------------------------------


# -- Read Data 
combinedGeneDataFilename = paste0(outputFolderName, filePrefix, "combinedGeneResultsWithAlternates.rds")
if(file.exists(combinedGeneDataFilename)){
  useAlternates = T
  combinedResults = readRDS(combinedGeneDataFilename)
  
}else{
  combinedGeneDataFilename = paste0(outputFolderName, filePrefix, "combinedGeneResults.rds")
  combinedResults = readRDS(combinedGeneDataFilename)
}

if(useAlternates){
  if(useMainAndRobust){
    significanceColumns = names(combinedResults)[grep("significantRobust", names(combinedResults))]
  }else{
    significanceColumns = names(combinedResults)[grep("PadjNumSignificant", names(combinedResults))]
  }

}else{
  if(usePermulations){
    significanceColumns = names(combinedResults)[grep("permSignificant", names(combinedResults))]
  }else{
    significanceColumns = names(combinedResults)[grep("unpermSignificant", names(combinedResults))]
  }
  #backward compatibility to load in old files if the new names don't exist 
  if(length(grep("permSignificant", names(combinedResults))) == 0 && length(grep("unpermSignificant", names(combinedResults))) == 0){
    significanceColumns = names(combinedResults)[grep("significant", names(combinedResults))]
  }
}






geneSignificanceResults = combinedResults[, names(combinedResults) %in% significanceColumns]

if(usingGo){
  combinedGODataFilename = paste0(outputFolderName, filePrefix, "combinedGOResultsWithAlternates-", geneSet, ".rds")
  if(file.exists(combinedGODataFilename) &useAlternates){
    GoCombinedResults = readRDS(combinedGODataFilename)
  }else{
    combinedGODataFilename = paste0(outputFolderName, filePrefix, "combinedGOResults-", geneSet, ".rds")
    GoCombinedResults = readRDS(combinedGODataFilename)
  }
  

  if(useAlternates){
    if(useMainAndRobust){
      GoSignificanceColumns = names(combinedResults)[grep("significantRobust", names(combinedResults))]
    }else{
      GoSignificanceColumns = names(combinedResults)[grep("PadjNumSignificant", names(combinedResults))]
    }
  }else{
    if(usePermulations){
      GoSignificanceColumns = names(GoCombinedResults)[grep("permSignificant", names(GoCombinedResults))]
    }else{
      GoSignificanceColumns = names(GoCombinedResults)[grep("unpermSignificant", names(GoCombinedResults))]
    }
    #backward compatibility to load in old files if the new names don't exist 
    if(length(grep("permSignificant", names(GoCombinedResults))) == 0 && length(grep("unpermSignificant", names(GoCombinedResults))) == 0){
      GoSignificanceColumns = names(GoCombinedResults)[grep("significant", names(GoCombinedResults))]
    }
  }
  GoSignificanceResults = GoCombinedResults[, names(GoCombinedResults) %in% GoSignificanceColumns]
}







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

# -- Make proportional venn diagram via eulerr ---
{ 
  vennPhenotypes = (strsplit(vennDiagramSet, "-"))
  commonBackground = Reduce(intersect, vennPhenotypes)
  commonBackground = substr(commonBackground, 1,1)
  # Make required functions 
  
  trimSignificanceToVenn = function(significanceResults){
    trimableComparisions = gsub("-PadjNumSignificant", "", names(significanceResults))
    trimableComparisions = gsub("significantRobust", "", trimableComparisions)
    trimableComparisions = gsub("Significant", "", trimableComparisions)
    trimableComparisions = gsub("significant", "", trimableComparisions)
    trimableComparisions = gsub("unperm", "", trimableComparisions)
    trimableComparisions = gsub("perm", "", trimableComparisions)
    trimableComparisions = gsub("-", "", trimableComparisions)
    trimableComparisions = addDashes(trimableComparisions)
    trimableComparisions = sapply(trimableComparisions, replacePrefixWithName)
    vennSignificanceResults = significanceResults[match(vennDiagramSet, trimableComparisions)]
  }
  
  
  makeVennPlot = function(vennInputDataframe, mainTitle, plot = T){
    
    if(useAlternates & !useMainAndRobust){
      set1 <- vennInputDataframe[1] > significanceCutoff
      set2 <- vennInputDataframe[2] > significanceCutoff
      set3 <- vennInputDataframe[3] > significanceCutoff
    }else{
      # Build logical vectors for each set
      set1 <- vennInputDataframe[1] == TRUE
      set2 <- vennInputDataframe[2] == TRUE
      set3 <- vennInputDataframe[3] == TRUE
    }

    
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
    
    comparisonPrefixes = gsub("Significant", "", names(vennInputDataframe))
    comparisonPrefixes = gsub("-significantRobust", "", comparisonPrefixes)
    comparisonPrefixes = gsub("-PadjNum", "", comparisonPrefixes)
    comparisonPrefixes = gsub("significant", "", comparisonPrefixes)
    comparisonPrefixes = gsub("-unperm", "", comparisonPrefixes)
    comparisonPrefixes = gsub("-perm", "", comparisonPrefixes)
    comparisonPrefixes = gsub("-", "", comparisonPrefixes)
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
  
  makeDirectionalResults = function(directionData, usedStatCols){
    sigDirectionDataPositive = directionData
    sigDirectionDataNegative = directionData
    for(i in usedStatCols){
      currentPrefix = substr(i, 1, 2)
      
      # get matching significant column
      relevantCols = grep(currentPrefix, names(directionData))
      statCol = names(directionData)[relevantCols[which(!names(directionData)[relevantCols] %in% significanceColumns)]]
      sigCol = names(directionData)[relevantCols[which(names(directionData)[relevantCols] %in% significanceColumns)]]
      
      #make sure that the direction is reported the same because the common background is in the same order
      positionOfBackground = as.integer(regexpr(commonBackground, currentPrefix))
      if(positionOfBackground == 2){
        sigDirectionDataPositive[[statCol]] = -1*sigDirectionDataPositive[[statCol]] 
        sigDirectionDataNegative[[statCol]] = -1*sigDirectionDataNegative[[statCol]] 
      }
      
      # set significance to F when the stat is wrong 
      
      sigDirectionDataPositive[[sigCol]][sigDirectionDataPositive[[statCol]] < 0] <- F
      
      sigDirectionDataNegative[[sigCol]][sigDirectionDataNegative[[statCol]] > 0] <- F
      
    } 
    sigDirectionDataNegative = sigDirectionDataNegative[, -which(names(sigDirectionDataNegative) %in% usedStatCols)]
    sigDirectionDataPositive = sigDirectionDataPositive[, -which(names(sigDirectionDataPositive) %in% usedStatCols)]
    directionalResults = list(sigDirectionDataPositive, sigDirectionDataNegative)
    return(directionalResults)
  }
  

}


#Make the directional data if using it 
if(makeDirectional){
  statColumns = names(combinedResults)[grep("Rho$", names(combinedResults))]
  geneSignificanceResultsDirection = combinedResults[, names(combinedResults) %in% c(significanceColumns, statColumns)]
  
  geneDirectionalSignificanceResults = makeDirectionalResults(geneSignificanceResultsDirection, statColumns)
  genePositiveSignificance = geneDirectionalSignificanceResults[[1]]
  geneNegativeSignificance = geneDirectionalSignificanceResults[[2]]
  
  if(usingGo){
    GoStatColumns = names(GoCombinedResults)[grep("stat", names(GoCombinedResults))]
    GoSignificanceResultsDirection = GoCombinedResults[, names(GoCombinedResults) %in% c(GoSignificanceColumns, GoStatColumns)]
    
    GoDirectionalSignificanceResults = makeDirectionalResults(GoSignificanceResultsDirection, GoStatColumns)
    GoPositiveSignificance = GoDirectionalSignificanceResults[[1]]
    GoNegativeSignificance = GoDirectionalSignificanceResults[[2]]
  }
}

if(useAlternates){
  mainTitle = paste0("All Genes")
  mainGoTitle = paste0("All Gene sets")
}else{
  mainTitle = paste0("All Genes (p.adj < ", significanceCutoff, ")")
  mainGoTitle = paste0("All Gene sets (p.adj < ", significanceCutoff, ")")
}


#Make the main plots 
if(length(pairwiseSets)==3){ #can simply run directly if only running on three categories. 
  vennColorset = geneVennColorset
  geneVenn = makeVennPlot(geneSignificanceResults, mainTitle)
  vennColorset = goVennColorset
  if(usingGo){goVenn = makeVennPlot(GoSignificanceResults, mainGoTitle)}
}else if(length(vennDiagramSet)==3){ #if a correct set of three categories has been given for the venn diagram set, narrow down a larger selection to that set, then run the vennDiagram. 
  vennGeneSignificanceResults = trimSignificanceToVenn(geneSignificanceResults)
  vennColorset = geneVennColorset
  geneVenn = makeVennPlot(vennGeneSignificanceResults, mainTitle)
  if(usingGo){
    vennGoSignificanceResults = trimSignificanceToVenn(GoSignificanceResults)
    vennColorset = goVennColorset
    goVenn = makeVennPlot(vennGoSignificanceResults, mainGoTitle)
  }
  
}
if(usingGo){
  combinedVenn = grid.arrange(geneVenn, goVenn, nrow = 1, padding = unit(1, "line"))
}else{
  combinedVenn = geneVenn
}

if(makeDirectional){
 
  if(useAlternates){
    directionalTitle = paste0(" Genes")
    directionalGoTitle = paste0(" Gene sets")
  }else{
    directionalTitle = paste0(" Genes (p.adj < ", significanceCutoff, ")")
    directionalGoTitle = paste0(" Gene sets (p.adj < ", significanceCutoff, ")")
  }
  
   
  if(length(pairwiseSets)==3){ #can simply run directly if only running on three categories. 
    
    #Positive
    vennColorset = geneVennColorset
    geneVennPositive = makeVennPlot(genePositiveSignificance, paste0(positiveLabel, directionalTitle))
    vennColorset = goVennColorset
    if(usingGo){goVennPositive = makeVennPlot(GoPositiveSignificance, paste0(positiveLabel, directionalGoTitle))}
    
    #Negtive
    vennColorset = geneVennColorset
    geneVennNegative = makeVennPlot(geneNegativeSignificance, paste0(negativeLabel, directionalTitle))
    vennColorset = goVennColorset
    if(usingGo){goVennNegative = makeVennPlot(GoNegativeSignificance, paste0(negativeLabel,directionalGoTitle))}
  }else if(length(vennDiagramSet)==3){ #if a correct set of three categories has been given for the venn diagram set, narrow down a larger selection to that set, then run the vennDiagram. 
    
    #Positive
    vennGeneSignificanceResultsPositive = trimSignificanceToVenn(genePositiveSignificance)
    vennColorset = geneVennColorset
    geneVennPositive = makeVennPlot(vennGeneSignificanceResultsPositive, paste0(positiveLabel,directionalTitle))
    if(usingGo){
      vennGoSignificanceResultsPositive = trimSignificanceToVenn(GoPositiveSignificance)
      vennColorset = goVennColorset
      goVennPositive = makeVennPlot(vennGoSignificanceResultsPositive, paste0(positiveLabel,directionalGoTitle))
    }
    
    #Negative
    vennGeneSignificanceResultsNegative = trimSignificanceToVenn(geneNegativeSignificance)
    vennColorset = geneVennColorset
    geneVennNegative = makeVennPlot(vennGeneSignificanceResultsNegative, paste0(negativeLabel, directionalTitle))
    if(usingGo){
      vennGoSignificanceResultsNegative = trimSignificanceToVenn(GoNegativeSignificance)
      vennColorset = goVennColorset
      goVennNegative = makeVennPlot(vennGoSignificanceResultsNegative, paste0(negativeLabel, directionalGoTitle))
    }
    
  }
  
  if(usingGo){
    combinedDirectionVenn = grid.arrange(geneVennPositive, goVennPositive, geneVennNegative,  goVennNegative, nrow = 2, padding = unit(1, "line"))
    positiveResultsNumber = sum(apply(vennGoSignificanceResultsPositive, 1, function(x) any(x == TRUE, na.rm = TRUE)))
    negativeResultsNumber = sum(apply(vennGoSignificanceResultsNegative, 1, function(x) any(x == TRUE, na.rm = TRUE)))
    ratio = positiveResultsNumber/negativeResultsNumber
    combinedDirectionVenn = grid.arrange(goVennPositive, goVennNegative, nrow = 2, padding = unit(1, "line"), heights = c(ratio, 1))
  }else{
    positiveResultsNumber = sum(apply(vennGeneSignificanceResultsPositive, 1, function(x) any(x == TRUE, na.rm = TRUE)))
    negativeResultsNumber = sum(apply(vennGeneSignificanceResultsNegative, 1, function(x) any(x == TRUE, na.rm = TRUE)))
    ratio = positiveResultsNumber/negativeResultsNumber
    combinedDirectionVenn = grid.arrange(geneVennPositive, geneVennNegative, nrow = 2, padding = unit(1, "line"), heights = c(ratio, 1))
  }
  
  if(usingGo){
    allVenn = grid.arrange(geneVenn,goVenn, combinedDirectionVenn, nrow = 1, padding = unit(1, "line"))
    
    vennDiagramFilename = paste0(outputFolderName, filePrefix, "VennDiagram", geneSet, ".pdf")
    pdf(vennDiagramFilename, height = 12, width = 16)
    plot(combinedVenn)
    dev.off()
    
    vennDirectionalDiagramFilename = paste0(outputFolderName, filePrefix, "VennDirectionalDiagram", geneSet, ".pdf")
    pdf(vennDirectionalDiagramFilename, height = 12, width = 12)
    plot(combinedDirectionVenn)
    dev.off()
    
    vennCombinedDiagramFilename = paste0(outputFolderName, filePrefix, "VennCombinedDiagram", geneSet, ".pdf")
    pdf(vennCombinedDiagramFilename, height = 16, width = 36)
    plot(allVenn)
    dev.off()
    
    
    largerOfDirections = max(positiveResultsNumber, negativeResultsNumber)
    totalResultsNumber = sum(apply(vennGoSignificanceResults, 1, function(x) any(x == TRUE, na.rm = TRUE)))
    totalRatio = totalResultsNumber/largerOfDirections
    allVenn = grid.arrange(goVenn, combinedDirectionVenn, nrow = 1, padding = unit(1, "line"), widths = c(totalRatio, 1))
    
    vennDiagramFilename = paste0(outputFolderName, filePrefix, "VennDiagram", geneSet, ".pdf")
    pdf(vennDiagramFilename, height = 6, width = 6)
    plot(combinedVenn)
    dev.off()
    
    vennDirectionalDiagramFilename = paste0(outputFolderName, filePrefix, "VennDirectionalDiagram", geneSet, ".pdf")
    pdf(vennDirectionalDiagramFilename, height = 6, width = 6)
    plot(combinedDirectionVenn)
    dev.off()
    
    vennCombinedDiagramFilename = paste0(outputFolderName, filePrefix, "VennCombinedDiagram", geneSet, ".pdf")
    pdf(vennCombinedDiagramFilename, height = 12, width = 19.2)
    plot(allVenn)
    dev.off()
    
  }else{
    largerOfDirections = max(positiveResultsNumber, negativeResultsNumber)
    totalResultsNumber = sum(apply(vennGeneSignificanceResults, 1, function(x) any(x == TRUE, na.rm = TRUE)))
    totalRatio = totalResultsNumber/largerOfDirections
    allVenn = grid.arrange(combinedVenn, combinedDirectionVenn, nrow = 1, padding = unit(1, "line"), widths = c(totalRatio, 1))
    grid.text("A", x = unit(0.02, "npc"), y = unit(0.98, "npc"), just = c("left", "top"), gp = gpar(fontsize = 16, fontface = "bold"))
    grid.text("B", x = unit(totalRatio/ (totalRatio +1), "npc"), y = unit(0.98, "npc"), just = c("left", "top"), gp = gpar(fontsize = 16, fontface = "bold"))
    grid.text("C", x = unit(totalRatio/ (totalRatio +1), "npc"), y = unit(ratio/(ratio+1), "npc"), just = c("left", "top"), gp = gpar(fontsize = 16, fontface = "bold"))
    
    
    vennDiagramFilename = paste0(outputFolderName, filePrefix, "VennDiagram", geneSet, ".pdf")
    pdf(vennDiagramFilename, height = 6, width = 6)
    plot(combinedVenn)
    dev.off()
    
    vennDirectionalDiagramFilename = paste0(outputFolderName, filePrefix, "VennDirectionalDiagram", geneSet, ".pdf")
    pdf(vennDirectionalDiagramFilename, height = 6, width = 6)
    plot(combinedDirectionVenn)
    dev.off()
    
    vennCombinedDiagramFilename = paste0(outputFolderName, filePrefix, "VennCombinedDiagram", geneSet, ".pdf")
    pdf(vennCombinedDiagramFilename, height = 8, width = 13)
    plot(allVenn)
    grid.text("A", x = unit(0.02, "npc"), y = unit(0.99, "npc"), just = c("left", "top"), gp = gpar(fontsize = 24, fontface = "bold"))
    grid.text("B", x = unit(totalRatio/ (totalRatio +1), "npc"), y = unit(0.99, "npc"), just = c("left", "top"), gp = gpar(fontsize = 24, fontface = "bold"))
    grid.text("C", x = unit(totalRatio/ (totalRatio +1), "npc"), y = unit(ratio/(ratio+1), "npc"), just = c("left", "top"), gp = gpar(fontsize = 24, fontface = "bold"))
    dev.off()
  }
  
  
  
  

  
  
}else{
  vennDiagramFilename = paste0(outputFolderName, filePrefix, "VennDiagram", geneSet, ".pdf")
  pdf(vennDiagramFilename, height = 12, width = 16)
  plot(combinedVenn)
  dev.off() 
}
}



























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
