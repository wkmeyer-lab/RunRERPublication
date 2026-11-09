clusterRun = F
clusterRun = T
if(clusterRun){.libPaths("/share/ceph/wym219group/shared/libraries/R4")} #add path to custom libraries to searched locations

library(RERconverge)
library(tools)
source("Src/Reu/cmdArgImport.R")
library(data.table)
source("Src/Reu/trimGmtFilenameToSetname.R")

# -- Command arguments list
# r = filePrefix                                                               This is a prefix used to organize and separate files by analysis run. Always required. 
# v = <T or F>                                                                 This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# g = geneSet                                                                  This is the gmtset(as shown in the filename) of enrichment runs to be combined. 

#----------------

args = c('r=ComplexDietCentralAnalysisSimplify2', 's=g', 'g=gene', 'l=170')
args = c('r=ComplexDietCentralAnalysisSimplify2', 's=g', 'g=KeggReactome', 'l=170')





{
  # --- Standard start-up code ---
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
}



# --- Argument Imports ---
# Defaults
geneSets = NULL
usingGo = F
usingGene = T
validGeneArguments = c("NULL", "gene", "Gene")


#genesets
if(!all(is.na(cmdArgImport('g')))){
  geneSets = cmdArgImport('g')
  if(is.null(geneSets) | any(geneSets %in% validGeneArguments)){
    usingGene = T
    message("Null geneset specified or gene specified, performing gene correlation.")
  }else{
    usingGene = F
  }
  
  if(any(!geneSets %in% validGeneArguments)){
    usingGo = T
    geneSets = trimGmtFilenameToSetname(geneSets)
  }else{
    usingGo = F
  }
  
}else{
  message("No geneset specified or gene specified, performing gene correlation.")
}


# --- Main code ---

primaryOutputFolderName = outputFolderName
primaryFilePrefix = filePrefix
#---
if(usingGene){
  #load correlations
  PairwiseCorrelationsFilenames = list.files(paste0(outputFolderName, "Alternates/"), pattern = "PairwiseCorrelationFile.rds")
  PairwiseCorrelationsFilenames = PairwiseCorrelationsFilenames[-grep("Merged", PairwiseCorrelationsFilenames)]
  PairwiseCorrelationsFilenames = paste0(outputFolderName, "Alternates/", PairwiseCorrelationsFilenames)
  
  pairwiseCorrelationRuns = lapply(PairwiseCorrelationsFilenames, readRDS)
  
  
  
  overallCorrelationsFilenames = list.files(paste0(outputFolderName, "Alternates/Overall"), pattern = "CorrelationFile.rds")
  overallCorrelationsFilenames = overallCorrelationsFilenames[-grep("Merged", overallCorrelationsFilenames)]
  overallCorrelationsFilenames = paste0("Overall/", overallCorrelationsFilenames)
  overallCorrelationsFilenames = paste0(outputFolderName, "Alternates/", overallCorrelationsFilenames)
  
  overallCorrelationRuns = lapply(overallCorrelationsFilenames, readRDS)
  overallCorrelationRuns = lapply(overallCorrelationRuns, function(df) df[, -2, drop = FALSE])
  names(overallCorrelationRuns) = rep("Overall-Overall", length(overallCorrelationRuns))
  
  
  numberOfPairwiseComparisions = length(pairwiseCorrelationRuns[[1]])
  correlationRuns = pairwiseCorrelationRuns
  
  for (i in seq_along(correlationRuns)) {
    correlationRuns[[i]][(numberOfPairwiseComparisions+1)] = overallCorrelationRuns[1]
    names(correlationRuns[[i]])[(numberOfPairwiseComparisions+1)] = "Overall"
  }
  
  
  #----
  #Combine correlations
  #----
  
  
  combinedCorrelations <- lapply(names(correlationRuns[[1]]), function(nm) {
    
    # Extract list of data frames for this comparison
    dfs <- lapply(correlationRuns, function(x) x[[nm]])
    
    n <- nrow(dfs[[1]])
    
    data.frame(
      Rho   = I(lapply(seq_len(n), function(i) sapply(dfs, function(df) df$Rho[i]))),
      P     = I(lapply(seq_len(n), function(i) sapply(dfs, function(df) df$P[i]))),
      p.adj = I(lapply(seq_len(n), function(i) sapply(dfs, function(df) df$p.adj[i])))
    )
  })
  names(combinedCorrelations) <- names(correlationRuns[[1]])
  
  
  combinedCorrelations <- lapply(combinedCorrelations, function(df) {
    
    #Number of significant
    df$PNumSignificant<- sapply(df$P, function(x) sum(x < 0.05, na.rm = TRUE))
    df$PadjNumSignificant <- sapply(df$p.adj, function(x) sum(x < 0.05, na.rm = TRUE))
    
    #Max difference
    df$RhoMaxDiff <- sapply(df$Rho, function(x) {
      if (all(is.na(x))) return(NA_real_)
      max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
    })
    
    df$PMaxDiff <- sapply(df$P, function(x) {
      if (all(is.na(x))) return(NA_real_)
      max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
    })
    
    df$PadjMaxDiff <- sapply(df$p.adj, function(x) {
      if (all(is.na(x))) return(NA_real_)
      max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
    })
    
    
    # -------- Rho stats --------
    df$RhoMean   <- sapply(df$Rho, function(x) mean(x, na.rm = TRUE))
    df$RhoMedian <- sapply(df$Rho, function(x) median(x, na.rm = TRUE))
    df$RhoSD     <- sapply(df$Rho, function(x) sd(x, na.rm = TRUE))
    df$RhoCI95   <- sapply(df$Rho, function(x) {
      x <- x[!is.na(x)]
      if (length(x) < 2) return(NA_real_)
      1.96 * sd(x) / sqrt(length(x))
    })
    df$RhoIQR    <- sapply(df$Rho, function(x) IQR(x, na.rm = TRUE))
    
    # -------- P stats --------
    df$PMean   <- sapply(df$P, function(x) mean(x, na.rm = TRUE))
    df$PMedian <- sapply(df$P, function(x) median(x, na.rm = TRUE))
    df$PSD     <- sapply(df$P, function(x) sd(x, na.rm = TRUE))
    df$PCI95   <- sapply(df$P, function(x) {
      x <- x[!is.na(x)]
      if (length(x) < 2) return(NA_real_)
      1.96 * sd(x) / sqrt(length(x))
    })
    df$PIQR    <- sapply(df$P, function(x) IQR(x, na.rm = TRUE))
    
    # -------- p.adj stats --------
    df$PadjMean   <- sapply(df$p.adj, function(x) mean(x, na.rm = TRUE))
    df$PadjMedian <- sapply(df$p.adj, function(x) median(x, na.rm = TRUE))
    df$PadjSD     <- sapply(df$p.adj, function(x) sd(x, na.rm = TRUE))
    df$PadjCI95   <- sapply(df$p.adj, function(x) {
      x <- x[!is.na(x)]
      if (length(x) < 2) return(NA_real_)
      1.96 * sd(x) / sqrt(length(x))
    })
    df$PadjIQR    <- sapply(df$p.adj, function(x) IQR(x, na.rm = TRUE))
    
    df
    
    
    df
  })
  
  
  saveRDS(combinedCorrelations, paste0(primaryOutputFolderName, primaryFilePrefix, "AlternatesCombinedCorrelations.rds"))
  
  message("Completed gene combinations")
}

#----
#Combine enrichments
#----
for(j in 1:length(geneSets)){
  geneSet = geneSets[j]
  if(geneSet %in% validGeneArguments){
    isGo = F
    next()
  }else{
    isGo = T
  }

  alternateSets = readRDS(paste0(primaryOutputFolderName, primaryFilePrefix, "AlternatePruningSpecies.rds"))
  
  
  subdirectories = list.dirs(paste0(outputFolderName, "Alternates"))
  subdirectories = subdirectories[-which(subdirectories == paste0(outputFolderName, "Alternates"))]
  subdirectories = gsub(paste0(outputFolderName, "Alternates/"), "", subdirectories)
  
  pairwiseSets = subdirectories
  
  

  outputFolderName = paste0(outputFolderName, "Alternates/")
  
  

  for(j in 1:length(alternateSets)){
    message(j)
    currentAlternate = j
    filePrefix = paste0("Alternate", j, primaryFilePrefix)

      GOResults = list()
      for(i in 1:length(pairwiseSets)){
        currentSet = pairwiseSets[i]
        correlationSubsetName = gsub("-", " - ", currentSet)
        correlationPrefix = paste(substr(strsplit(currentSet, split = "-")[[1]],1,1), collapse = '')
        
        goFilename = paste0(outputFolderName, currentSet, "/", filePrefix, currentSet,"Enrichment-", geneSet, ".rds")
        currentGoData = readRDS(goFilename)[[1]]
        
        
        
        
        #names(currentGoData) = paste0(correlationPrefix, "-", names(currentGoData))
        
  
        GOResults[[i]] = currentGoData
        names(GOResults)[i] = currentSet
        
      }
      rm(currentGoData)
      
      GoCombinedResults = NA
      longestGOResult = which(sapply(GOResults, nrow) == max(sapply(GOResults, nrow)))[1]
      for(i in 1:length(GOResults)){
        if(all(rownames(GOResults[[i]]) == rownames(GOResults[[longestGOResult]]))){
          cat("Combining GO Data", i, "\n")
          GoCombinedResults = cbind(GoCombinedResults, GOResults[[i]])
        }else{
          tryCatch({
            cat("Mismatch in order; reordering current data to first data; Combining GO Data", i, "\n")
            
            which(!rownames(GOResults[[longestGOResult]]) %in% rownames(GOResults[[i]]))
            
            matchedOrder = match(rownames(GOResults[[longestGOResult]]), rownames(GOResults[[i]]))
            #matchedOrder = matchedOrder[-which(is.na(matchedOrder))]
            GOResults[[i]] = GOResults[[i]][matchedOrder,]
            #GoCombinedResults = cbind(GoCombinedResults, GOResults[[i]])
          },error = function(e){      
            stop("Rownames of GO results are not the same and cannot be reordered, there is an issue with the GO data.")
          }
          
          )

        }
      }
      #GoCombinedResults = GoCombinedResults[,-1]
      
     
      
      if(saveCombinedData){
        combinedGODataFilename = paste0(outputFolderName, filePrefix, "combinedGOResults-", geneSet)
        saveRDS(GOResults, paste0(combinedGODataFilename, ".rds"))
      }
  }
  
  message("Completed GO combinations for individual alternates")
  
  #---- 
  #combine alternates
  
  EnrichmentFilenames = list.files(paste0(primaryOutputFolderName, "Alternates/"), pattern = paste0("combinedGOResults-", geneSet, ".rds"))
  
  EnrichmentFilenames = paste0(primaryOutputFolderName, "Alternates/", EnrichmentFilenames)
  
  EnrichmentRuns = lapply(EnrichmentFilenames[1:length(EnrichmentFilenames)], readRDS)
  
  #Set the order of the enrichment runs correctly 
  


  
  #First, confirm all rownames are the same within the enrichment runs 
  errorSet = NA
  for(i in 1:100){
    mainRowNames = rownames(EnrichmentRuns[[i]][[1]])
    for(j in 2:length(EnrichmentRuns[[i]])){
      newRowNames = rownames(EnrichmentRuns[[i]][[j]]) 
      isEqual = all(newRowNames ==mainRowNames)
      if(!isEqual){
        errorCode = paste0(i, "-", j)
        names(errorCode) = i
        errorSet = append(errorSet, errorCode)
      }
    }
  }
  if(!(all(is.na(errorSet)))){
    message(paste("Problem in rownames of alternate:", unique(names(errorSet))))
    stop("THERE IS A DIFFERENCE OF ROWNAMES WITHIN AN ALTERNATE, ABORTING")
  }else{message("Rownames are matched WITHIN all alternates")}
  
  #Then, reorder all of the enrichments into that order
  rownameSets = list()
  for(i in 1:length(EnrichmentRuns)){
    currentRownames = rownames(EnrichmentRuns[[i]][[1]])
    rownameSets = append(rownameSets, list(currentRownames))
  }
  
  longestAlternateGoSetIndex = which(sapply(rownameSets, length) == max(sapply(rownameSets, length)))[1]
  longestAlternateGoSet = rownameSets[[longestAlternateGoSetIndex]]
  allIncludedGOSets = unique(unlist(rownameSets))
  
  #Check that the longest entry includes all of the contained GO sets 
  if(!all(allIncludedGOSets %in% longestAlternateGoSet)){
    stop("The longest GO set does not include all of the included GO sets, stoping combination.")
  }else{message("Longest GO set confirmed to include all GO sets processed in any alternate")}
  
  
  for(i in 1:length(EnrichmentRuns)){
    currentRownames = rownameSets[[i]]
    if(all(currentRownames == longestAlternateGoSet)){
      message(paste("Alternate ", i, " In correct order, continuing."))
    }else{
      message("Reordering Alternate set ", i, " to match longest GO set, confirmed to contain all processed GOs.")
      matchedOrder = match(longestAlternateGoSet, currentRownames)
      
      for(j in 1:length(EnrichmentRuns[[i]])){
        EnrichmentRuns[[i]][[j]] = EnrichmentRuns[[i]][[j]][matchedOrder, ]
      }
    }
  }
  
  
  
  
  combinedEnrichment <- lapply(seq_along(EnrichmentRuns[[1]]), function(g) {
    
    group_name <- names(EnrichmentRuns)[g]
    
    dfs <- lapply(EnrichmentRuns, function(run) run[[g]])
    
    n <- nrow(dfs[[1]])
    
    data.frame(
      stat = I(lapply(seq_len(n), function(i)
        sapply(dfs, function(df) df$stat[i])
      )),
      
      pval = I(lapply(seq_len(n), function(i)
        sapply(dfs, function(df) df$pval[i])
      )),
      
      p.adj = I(lapply(seq_len(n), function(i)
        sapply(dfs, function(df) df$p.adj[i])
      ))
    )
    
  })
  
  names(combinedEnrichment) <- names(EnrichmentRuns[[1]])
  
  
  combinedEnrichment <- lapply(combinedEnrichment, function(df) {
    
    df$PNumSignificant <- sapply(df$pval, function(x) sum(x < 0.05, na.rm = TRUE))
    df$PadjNumSignificant <- sapply(df$p.adj, function(x) sum(x < 0.05, na.rm = TRUE))
    
    
    df$statMaxDiff <- sapply(df$stat, function(x)
      if (all(is.na(x))) NA_real_ else diff(range(x, na.rm = TRUE))
    )
    
    df$pMaxDiff <- sapply(df$pval, function(x)
      if (all(is.na(x))) NA_real_ else diff(range(x, na.rm = TRUE))
    )
    
    df$padjMaxDiff <- sapply(df$p.adj, function(x)
      if (all(is.na(x))) NA_real_ else diff(range(x, na.rm = TRUE))
    )
    
    # -------- Rho stats --------
    df$statMean   <- sapply(df$stat, function(x) mean(x, na.rm = TRUE))
    df$statMedian <- sapply(df$stat, function(x) median(x, na.rm = TRUE))
    df$statSD     <- sapply(df$stat, function(x) sd(x, na.rm = TRUE))
    df$statCI95   <- sapply(df$stat, function(x) {
      x <- x[!is.na(x)]
      if (length(x) < 2) return(NA_real_)
      1.96 * sd(x) / sqrt(length(x))
    })
    df$statIQR    <- sapply(df$stat, function(x) IQR(x, na.rm = TRUE))
    
    # -------- P stats --------
    df$PMean   <- sapply(df$pval, function(x) mean(x, na.rm = TRUE))
    df$PMedian <- sapply(df$pval, function(x) median(x, na.rm = TRUE))
    df$PSD     <- sapply(df$pval, function(x) sd(x, na.rm = TRUE))
    df$PCI95   <- sapply(df$pval, function(x) {
      x <- x[!is.na(x)]
      if (length(x) < 2) return(NA_real_)
      1.96 * sd(x) / sqrt(length(x))
    })
    df$PIQR    <- sapply(df$pval, function(x) IQR(x, na.rm = TRUE))
    
    # -------- p.adj stats --------
    df$PadjMean   <- sapply(df$p.adj, function(x) mean(x, na.rm = TRUE))
    df$PadjMedian <- sapply(df$p.adj, function(x) median(x, na.rm = TRUE))
    df$PadjSD     <- sapply(df$p.adj, function(x) sd(x, na.rm = TRUE))
    df$PadjCI95   <- sapply(df$p.adj, function(x) {
      x <- x[!is.na(x)]
      if (length(x) < 2) return(NA_real_)
      1.96 * sd(x) / sqrt(length(x))
    })
    df$PadjIQR    <- sapply(df$p.adj, function(x) IQR(x, na.rm = TRUE))
    
    df
    
    df
  })
  
  for(i in 1:length(combinedEnrichment)){
    for(j in 1:length(combinedEnrichment[[i]])){
      names(combinedEnrichment[[i]][[j]]) = longestAlternateGoSet
    }
  }
  
  saveRDS(combinedEnrichment, paste0(primaryOutputFolderName, primaryFilePrefix, "AlternatesCombinedEnrichments-" , geneSet, ".rds"))
  
}
