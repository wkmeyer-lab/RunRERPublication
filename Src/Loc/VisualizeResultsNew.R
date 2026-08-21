# Fixes before publication: Automatic number of geneset detection 

# -- Libraries 
clusterRun = F
clusterRun = T
if(clusterRun){.libPaths("/share/ceph/wym219group/shared/libraries/R4")} #add path to custom libraries to searched locations
library(RERconverge)
library(ggplot2)
library(ggpmisc)
library(cowplot)
library(qvalue)
library(insight)
source("Src/Reu/cmdArgImport.R")


# -- Usage:
# This script generates several visualizations of the RER outputs, and then organizes them into a single PDF file. 

# -- Command arguments list
# r = filePrefix                                     This is a prefix used to organize and separate files by analysis run. Always required. 
# v = <T or F>                                       This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# g = <T OR F>                                       This sets if Gene Enrichment analysis is included. Default TRUE. 
# l = <"enrichmentListName" OR c("name1","name2")>   This defines the enrichment list(s) used in the GO analysis the script should display
# s = "subdirectoryName"                             This is used to specify a subdirectory for the analysis to be run in. Primarily used for the components of categorical results. 
# o = "EnrichmentSortCollum"                         This is used to specify what the enrichment plots should be sorted by 
# u = < T or F or A>                                 This is used to set if the enrichment plots should be ordered by decreasing or not, or A for absolute value 

# p = <T or F or B or C>                             This sets if the code should use permulated or unpermualted values, or both. "C" indicates categorical permulations, which are stored in the main file. 
# f = "permulationPvalueFileLocation.rds"            This is a manual override to specify the script use a specific Permulation p-value file. 
      #If using any permulation p-value file other than "MainCombined" with no run instance number, it must be specified manually.
              


#----------------

args = c('r=ComplexDietCentralAnalysis', 'p=F', 'g=T', 's=c("Carnivore-Herbivore", "Omnivore-Vertivore", "Insectivore-Vertivore", "Herbivore-Vertivore", "Insectivore-Omnivore", "Herbivore-Omnivore", "Herbivore-Insectivore", "Overall")' )


# --- Standard start-up code ---
if(clusterRun){args = commandArgs(trailingOnly = TRUE)}
{  # Bracket used for collapsing purposes
  #File Prefix
  if(!is.na(cmdArgImport('r'))){                                                #This cmdArgImport script is a way to import arguments from the command line. 
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
permulationDefaultFilename = "MainCombinedPermulationsPValue.rds" #This argument isn't actually updated from the command line, it is used to easily swap what the default permulation filename the script expects.
permulationDefaultFilename= "CombinedPrunedFastAllPermulationsPValue.rds"

usePermulations = FALSE
permulationPValOverride = NULL 
useGeneEnrichment = TRUE
useBoth = FALSE
useCategoricalPerms = FALSE
useSubdirectory = FALSE
subdirectoryValueList = NULL
j=1
enrichmentSortColumn = "p.adj"
enrichmentOrderDecrease = F

{ # Bracket used for collapsing purposes
  
  #Permulation use
  if(!is.na(cmdArgImport('p'))){
    usePermulations = cmdArgImport('p')
    if(usePermulations %in% c("CB", "cb")){
      usePermulations = TRUE
      useCategoricalPerms = TRUE 
      useBoth = TRUE
    }else if(usePermulations %in% c("C", "c")){
      usePermulations = TRUE
      useCategoricalPerms = TRUE 
    } else if(usePermulations %in% c("B", "b")){
      useBoth = TRUE
      usePermulations = TRUE
    }else{
      usePermulations = as.logical(usePermulations)
      if(is.na(usePermulations)){
        usePermulations = TRUE
        useBoth = FALSE
        message("Use Permualtions value not interpretable as logical. Did you remember to capitalize? Using TRUE.")
      }
    }
  }else{
    message("Use Permulations value not specified, using TRUE.")
  }
  
  #Permulation File Override 
  if(usePermulations){
    if(!is.na(cmdArgImport('f'))){
      permulationPValOverride = cmdArgImport('f')
      message("Using Manually specified permulation p-value file. NOTE this setting ignores the meta-combination argument.")
    }else{
      message("No Permulation pValue override, using standard CombinedPrunedFastAllPermulationsPValue.rds.")
    } 
  }
  
  #Gene Enrichment Use 
  if(!is.na(cmdArgImport('g'))){
    useGeneEnrichment = cmdArgImport('g')
    useGeneEnrichment = as.logical(useGeneEnrichment)
    if(is.na(useGeneEnrichment)){
      useGeneEnrichment = TRUE
      message("Include geneset enrichment value not interpretable as logical. Did you remember to capitalize? Using TRUE.")
    }
  }else{
    message("Include geneset enrichment value not specified, using TRUE.")
  }
  
  #Import subdirectory
  if(!any(is.na(cmdArgImport('s')))){
    useSubdirectory = TRUE
    subdirectoryValueList = cmdArgImport('s')
    message(paste("Using subdirectories", subdirectoryValueList, "."))
    
    #outputFolderName = paste(outputFolderName, subdirectoryValue, "/", sep="")
    
  }else{
    message("No subdirectory specified.")
  }
  
  #Import enrichment column
  if(!any(is.na(cmdArgImport('o')))){
    useSubdirectory = TRUE
    enrichmentSortColumn = cmdArgImport('o')
    message(paste("Sorting erichments by", enrichmentSortColumn, "."))
    
    #outputFolderName = paste(outputFolderName, subdirectoryValue, "/", sep="")
    
  }else{
    message("No Column specified, using p.adj.")
  }
  
  #Import erichment sort direction
  if(!any(is.na(cmdArgImport('u')))){
    useSubdirectory = TRUE
    enrichmentOrderDecrease = cmdArgImport('u')
    message(paste("Ordering erichments by ", enrichmentOrderDecrease, "."))
    
    #outputFolderName = paste(outputFolderName, subdirectoryValue, "/", sep="")
    
  }else{
    message("No enrichment order speciefied, sorting by decreasing = F.")
  }
}

#                   ------- Code Body -------- 


# ------ Import the Data ------ 
for(j in 1:length(subdirectoryValueList)){
  outputFolderName = paste("Output/",filePrefix,"/", sep = "")
  message(paste("Using subdirectory", subdirectoryValueList[j], "."))
  if(useSubdirectory){
    outputFolderName = paste(outputFolderName, subdirectoryValueList[j], "/", sep="")
  }
  subdirectoryValue = subdirectoryValueList[j]
  
  if(useCategoricalPerms){
    correlationFileLocation = paste(outputFolderName, filePrefix, subdirectoryValue, "PermulationsCorrelationFile.rds", sep= "")
  }else{
    correlationFileLocation = paste(outputFolderName, filePrefix, subdirectoryValue, "CorrelationFile.rds", sep= "")
  }
  correlData = readRDS(correlationFileLocation)                            #Import the correlation data (non-permulated)
  
  # - Permulations - 
  if(usePermulations){
    if(useCategoricalPerms){
      correlData$permPValue = correlData$permP
    }else{
      if(!is.null(permulationPValOverride)){
        permulationFileLocation = permulationPValOverride
      }else{
        permulationFileLocation = paste(outputFolderName, filePrefix, permulationDefaultFilename, sep= "")
      }
      permData = readRDS(permulationFileLocation) 
      if(ncol(permData) ==1 ){
        correlData$permPValue = permData                       #Add a collumn to the data with the permulation p Values  
      }else{
        correlData$permPValue = permData$permpval                      #Add a collumn to the data with the permulation p Values  
      }
      
    }
  }
  
  # - Q values - 
  qValueNoperm = qvalue(p = correlData$p.adj)
  correlData$qValueNoperm = qValueNoperm$qvalues
  
  if(usePermulations){
    qValuePerm = qvalue(p = correlData$permPValue)
    correlData$qValuePerm = qValuePerm$qvalues
  }
  
  # - Negative and Positive Rho 
  correlDataPositive = correlData[which(correlData$Rho > 0),]
  correlDataNegative = correlData[which(correlData$Rho < 0),]
  correlDataZero = correlData[which(correlData$Rho == 0),]
  
  # - Enrichment results - 
  if(useGeneEnrichment){
    outputFolderFilenames = list.files(paste("./",outputFolderName, sep=""))
    enrichmentFileNames = outputFolderFilenames[grep("Enrichment-", outputFolderFilenames)]
    if(length(enrichmentFileNames) == 0){
      enrichmentRange = NA
      useGeneEnrichment = FALSE
      message("No enrichment files found. Deactivating enrichment analysis processing.")
    }else{
      enrichmentRange = length(enrichmentFileNames)
      enrichmentResultSets = NULL
      for(i in 1:enrichmentRange){
        EnrichmentSetNumber = paste("enrichment", i, sep="")
        EnrichmentData = readRDS(paste(outputFolderName, enrichmentFileNames[i], sep=""))
        enrichmentResultSets[i] = EnrichmentData
        names(enrichmentResultSets)[i] = names(EnrichmentData)
        assign(EnrichmentSetNumber, EnrichmentData)
      }
    }
  }
  
  
  # ------ Make plots ------ 
  
  # - p value histograms - 
  makePHistogram = function(data, column, titleVal){
    phist = ggplot(data, aes(x=data[[column]]))+
      geom_histogram(binwidth = 0.02, alpha=0.9, col="white", boundary=0)+
      scale_x_continuous(breaks = seq(0,1,0.1), lim = c(0,1))+
      labs(title = titleVal)+
      xlab("p Values")+
      theme(
            #axis.title.x = element_blank(),
            plot.title = element_text(size=18, hjust = 0.5),
            plot.title.position = "plot"
          )
  }
  
  postiveRhoNonpermHistogram = makePHistogram(correlDataPositive, "P", "Positive Rho Non-permulated")
  negativeRhoNonpermHistogram = makePHistogram(correlDataNegative, "P", "Negative Rho Non-permulated")
  
  if(usePermulations){
    postiveRhoPermHistogram = makePHistogram(correlDataPositive, "permPValue", "Positive Rho Permulated")
    negativeRhoPermHistogram = makePHistogram(correlDataNegative, "permPValue", "Negative Rho Permulated")
    
  }
  #I think I want to display positive and negative side-by-side, and permulated non-permulated above/below eachother.
  
  
  # - List plots - 
  makePvListPlot = function(data, column, length, titleVal){
    data = data[order(data[[column]]),]
    dataHead = data[1:length,]
    dataFront = dataHead[,1:2]
    dataBack = dataHead[,3:ncol(dataHead)]
    dataBack = format_table(dataBack, pretty_names = T, digits = "scientific3")
    dataMain = data.frame(rownames(dataHead))
    dataMain = append(dataMain, dataFront)
    names(dataMain)[1] = "Gene"
    dataMain = append(dataMain, dataBack)
    dataMain = as.data.frame(dataMain)
    listPlot = ggplot()+
      theme_void()+
      annotate(geom = "table",
               x=1,
               y=1,
               label = dataMain)+
      labs(title = titleVal)+
      theme(plot.title = element_text(size=18, hjust = 0.5, vjust=0))
  }
  
  topGenesPadj = makePvListPlot(correlData, "p.adj", 25, "Top genes by P-value non-permulated")
  topPositiveGenesPadj = makePvListPlot(correlDataPositive, "p.adj", 10,"Top Positive genes by P-value non-permulated")
  topNegativeGenesPadj = makePvListPlot(correlDataNegative, "p.adj", 10, "Top Negative genes by P-value non-permulated")
  
  topGenesNPQ = makePvListPlot(correlData, "qValueNoperm", 25, "Top genes by Q-Value non-permulated")
  topPositiveGenesNPQ = makePvListPlot(correlDataPositive,  "qValueNoperm",10, "Top Positive genes by Q-Value non-permulated")
  topNegativeGenesNPQ = makePvListPlot(correlDataNegative,  "qValueNoperm", 10,"Top Negative genes by Q-Value non-permulated")
  
  if(usePermulations){
    topGenesPerm = makePvListPlot(correlData, "permPValue", 25, "Top genes by P-value Permulated")
    topPositiveGenesPerm = makePvListPlot(correlDataPositive,  "permPValue", 10,"Top Positive genes by P-value Permulated")
    topNegativeGenesPerm = makePvListPlot(correlDataNegative,  "permPValue", 10,"Top Negative genes by P-value Permulated")
    
    topGenesPermQ = makePvListPlot(correlData, "qValuePerm", 25, "Top genes by Q-Value Permulated")
    topPositiveGenesPermQ = makePvListPlot(correlDataPositive,  "qValuePerm",10, "Top Positive genes by Q-Value Permulated")
    topNegativeGenesPermQ = makePvListPlot(correlDataNegative,  "qValuePerm",10, "Top Negative genes by Q-Value Permulated")
  }
  
  # - Gene Enrichment Plots - 
  if(useGeneEnrichment){
    enrichmentPlotSet = list()
    makeGeListPlot = function(data, column, length, titleVal, decreasing){
      setData = data[[1]]
      if(decreasing == "A"){
        setData = setData[order(abs(setData[[column]]), decreasing = T),]
      }else if(decreasing){
      setData = setData[order(abs(setData[[column]]), decreasing = T),]
      }else{
        setData = setData[order(abs(setData[[column]])),]
      }
      dataHead = setData[1:length,]
      dataHead$gene.vals= strsplit(dataHead$gene.vals, ",")
      for(i in 1:length){
        genesetlist = dataHead$gene.vals[[i]][1:7]
        genesetlistsingle = paste(genesetlist[1], genesetlist[2], genesetlist[3], genesetlist[4], genesetlist[5],genesetlist[6],genesetlist[7])
        dataHead$gene.vals[i] = genesetlistsingle
      }
      dataFront = dataHead[,c(1,4)]
      dataBack = dataHead[,c(2,3,5)]
      dataBack = format_table(dataBack, pretty_names = T, digits = "scientific3")
      dataMain = data.frame(substring(rownames(dataHead), 1, 40))
      dataMain = append(dataMain, dataFront)
      names(dataMain)[1] = "Geneset"
      dataMain = append(dataMain, dataBack)
      dataMain = as.data.frame(dataMain)
      listPlot = ggplot()+
        theme_void()+
        annotate(geom = "table",
                 x=1,
                 y=1,
                 label = dataMain)+
        labs(title = paste(names(data), titleVal))+
        theme(plot.title = element_text(size=18, hjust = 0.5, vjust=1))
    }
    for(i in 1:enrichmentRange){
      genesetPlotName = paste("genesetPlot", i, sep="")
      if(usePermulations){
        genesetPlot = makeGeListPlot(enrichmentResultSets[i], enrichmentSortColumn, 40, "Top pathways by permulation", enrichmentOrderDecrease)
      }else{
        genesetPlot = makeGeListPlot(enrichmentResultSets[i], enrichmentSortColumn, 40, "Top pathways by non-permulation", enrichmentOrderDecrease)
      }
      enrichmentPlotSet[[i]] = genesetPlot
      assign(genesetPlotName, genesetPlot)
      
    }
  }
  # ------- Make Plots ------- 
  
  # - Plot the top genes - 
  if(useBoth){
    headlineGenes = plot_grid(topGenesPerm, topGenesPermQ, topGenesPadj, topGenesNPQ, ncol = 2, nrow = 2)
    headlineRows = 2
  }else if(usePermulations){
    headlineGenes = plot_grid(topGenesPerm, topGenesPermQ, ncol = 2, nrow = 1)
    headlineRows = 1
  }else{
    headlineGenes = plot_grid(topGenesPadj, topGenesNPQ, ncol = 2, nrow = 1)
    headlineRows = 1
  }
  
  # - plot the histograms - 
  histograms = NULL
  if(useBoth){
    histograms = plot_grid(postiveRhoPermHistogram, negativeRhoPermHistogram, postiveRhoNonpermHistogram, negativeRhoNonpermHistogram, ncol = 2, nrow = 2)
    histogramRows = 2
  }else if(usePermulations){
    histograms = plot_grid(postiveRhoPermHistogram, negativeRhoPermHistogram, ncol = 2, nrow = 1)
    histogramRows = 1
  }else{
    histograms = plot_grid(postiveRhoNonpermHistogram, negativeRhoNonpermHistogram, ncol = 2, nrow = 1)
    histogramRows = 1
  }
  
  # - plot the signed genes - 
  if(useBoth){
    signedGenes = plot_grid(topPositiveGenesPerm, topNegativeGenesPerm, topPositiveGenesPadj, topNegativeGenesPadj, ncol = 2, nrow = 2)
    signedRows = 2
  }else if(usePermulations){
    signedGenes = plot_grid(topPositiveGenesPerm, topNegativeGenesPerm, ncol = 2, nrow = 1)
    signedRows = 1
  }else{
    signedGenes = plot_grid(topPositiveGenesPadj, topNegativeGenesPadj, ncol = 2, nrow = 1)
    signedRows = 1
  }
  
  # - plot enrichments - 
  enrichmentRows = 0 
  if(useGeneEnrichment){
    #enrichmentPlots = enrichmentPlotSet[1][1]
    #for(i in 2:enrichmentRange){
    #  enrichmentPlots = plot_grid(genesetPlot1, genesetPlot2, genesetPlot3, genesetPlot4, genesetPlot5, ncol = 1, nrow = 3)
    #}
    #enrichmentPlots= plot_grid(genesetPlot1, genesetPlot2, ncol = 1, nrow = 5)

    enrichmentPlots= plot_grid(genesetPlot3, genesetPlot2, genesetPlot1, genesetPlot4, genesetPlot5, genesetPlot6, genesetPlot7, ncol = 1, nrow = 7)

    enrichmentRows = pmax(length(enrichmentRange), enrichmentRange)
  }
  
  # ------ Output to pdf ------ 
  pdfRows = headlineRows+histogramRows+signedRows+enrichmentRows
  pdfLengthPerRow = 6.5
  pdfLength = pdfRows*pdfLengthPerRow
  
  outputPDFLocation = paste(outputFolderName, filePrefix, subdirectoryValue, "VisualizeOutput.pdf", sep= "") # this could be improved to be more dynamic
  
  pdf(file = outputPDFLocation, width = 15, height = pdfLength)
  
  mainplot = plot_grid(headlineGenes, histograms, signedGenes, ncol = 1, nrow = 3)
  print(mainplot)
  if(useGeneEnrichment){
    print(enrichmentPlots)
  }
 

  dev.off()
  
  
  outputPDFLocation = paste(outputFolderName, filePrefix, subdirectoryValue, "GOVisualization.pdf", sep= "") # this could be improved to be more dynamic
  pdf(file = outputPDFLocation, width = 15, height = pdfLength)
  if(useGeneEnrichment){
    print(enrichmentPlots)
  }
  dev.off()
  
  outputPDFLocation = paste(outputFolderName, filePrefix, subdirectoryValue, "GeneVisualization.pdf", sep= "") # this could be improved to be more dynamic
  pdfRows = headlineRows+histogramRows+signedRows
  pdfLengthPerRow = 6.5
  pdfLength = pdfRows*pdfLengthPerRow
  pdf(file = outputPDFLocation, width = 15, height = pdfLength)
  print(mainplot)
  dev.off()
  
    
  #if(useGeneEnrichment){
  #  for(i in 1:2){
  #    enrichmentPDFLocation = paste(outputFolderName, filePrefix, subdirectoryValue, "Geneset", names(enrichmentResultSets[i]), ".pdf", sep= "")
  #    ggsave(file = enrichmentPDFLocation, enrichmentPlotSet[[i]], width = 12, height = 10)
  #  }
  #}
}

#improvements for this script: 
  #Update output title to be more dynamic
  #fix justification of figure titles
  #Change row length to match length of content (shorten the signed row) 

