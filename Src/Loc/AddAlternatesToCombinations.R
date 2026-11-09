clusterRun = F
clusterRun = T
if(clusterRun){.libPaths("/share/ceph/wym219group/shared/libraries/R4")} #add path to custom libraries to searched locations
library(RERconverge)
library(purrr)
source("Src/Reu/cmdArgImport.R")

args = c("r=ComplexDietCentralAnalysis", "p=NULL", "g=gene", "s=T", "l=F")
args = c("r=ComplexDietCentralAnalysis", "p=NULL", "g=KeggReactome", "s=T", "l=F")

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
saveData = T
usingGene = T
alternateCutoff = 50 

{ # Bracket used for collapsing purposes
  
  #geneset
  if(!is.na(cmdArgImport('g'))){
    geneSet = cmdArgImport('g')
    if(is.null(geneSet) | geneSet == "NULL" | geneSet == "gene" | geneSet == "Gene"){
      usingGene = T
      message("Null geneset specified or gene specified, performing gene correlation.")
    }else{
      usingGo = T
      usingGene = F
    }
  }else{
    message("No geneset specified or gene specified, performing gene correlation.")
  }
  
  #cutoff
  if(!is.na(cmdArgImport('c'))){
    saveData = cmdArgImport('c')
  }else{
    message("pvalue cuttoff not specified, using 0.05")
  }  
  
  
  #saveData 
  if(!is.na(cmdArgImport('s'))){
    saveData = as.logical(cmdArgImport('s'))
  }else{
    message("saveData not specified, using TRUE")
  }
  
  #Use Permualtions 
  if(!is.na(cmdArgImport('l'))){
    usePermulations = as.logical(cmdArgImport('l'))
  }else{
    message("use Permulations not specified, using FLASE")
  }
}




#----- 
#Main Gene code



if(usingGene){
  #Read in the files
  combinedDataFilename = paste0(outputFolderName, filePrefix, "combinedGeneResults.rds")
  combinedData = readRDS(combinedDataFilename)

  alternatesFilename = paste0(outputFolderName, filePrefix, "AlternatesCombinedCorrelations.rds")
  alternatesData = readRDS(alternatesFilename)
  
  
  
  
}else{
  combinedGODataFilename = paste0(outputFolderName, filePrefix, "combinedGOResults-", geneSet, ".rds")
  combinedData = readRDS(combinedGODataFilename)
  
  alternatesFilename = paste0(outputFolderName, filePrefix, "AlternatesCombinedEnrichments-" , geneSet, ".rds")
  alternatesData = readRDS(alternatesFilename)
  
}
  

#Reformat the alternates   
alternateSummaryColumns = map(alternatesData, ~ .x[, 4:ncol(.x), drop = FALSE])

alternatesSingleDF <- imap(alternateSummaryColumns, function(df, nm) {
  parts <- strsplit(nm, "-")[[1]]
  prefix <- paste0(substr(parts[1], 1, 1),
                   substr(parts[2], 1, 1))
  
  df %>% rename_with(~ paste0(prefix, "-", .x))
}) %>%
  bind_cols()


#Add the alternate data to the main dataframe 
combinedData = cbind(combinedData, alternatesSingleDF)




{  #Reorder the columns

  cols <- colnames(combinedData)
  
  # get all prefixes (e.g., CH, CO, HI, etc.)
  prefixes <- unique(sub("-.*", "", cols))
  
  # function to order columns for one prefix
  
  if(usingGene){
    order_prefix <- function(p) {
      c(
        paste0(p, "-significant"),
        paste0(p, "-Rho"),
        paste0(p, "-RhoMedian"),
        paste0(p, "-RhoMean"),
        paste0(p, "-RhoSD"),
        paste0(p, "-RhoCI95"),
        paste0(p, "-RhoIQR"),
        paste0(p, "-RhoMaxDiff"),
        paste0(p, "-P"),
        paste0(p, "-PNumSignificant"),
        paste0(p, "-PMedian"),
        paste0(p, "-PMean"),
        paste0(p, "-PSD"),
        paste0(p, "-PCI95"),
        paste0(p, "-PIQR"),
        paste0(p, "-PMaxDiff"),
        paste0(p, "-p.adj"),
        paste0(p, "-PadjNumSignificant"),
        paste0(p, "-PadjMedian"),
        paste0(p, "-PadjMean"),
        paste0(p, "-PadjSD"),
        paste0(p, "-PadjCI95"),
        paste0(p, "-PadjIQR"),
        paste0(p, "-PadjMaxDiff")
      )
    }
  }else{
    order_prefix <- function(p) {
      c(
        paste0(p, "-significant"),
        paste0(p, "-num.genes"),
        paste0(p, "-gene.vals"),
        paste0(p, "-stat"),
        paste0(p, "-statMedian"),
        paste0(p, "-statMean"),
        paste0(p, "-statSD"),
        paste0(p, "-statCI95"),
        paste0(p, "-statIQR"),
        paste0(p, "-statMaxDiff"),
        paste0(p, "-P"),
        paste0(p, "-PNumSignificant"),
        paste0(p, "-PMedian"),
        paste0(p, "-PMean"),
        paste0(p, "-PSD"),
        paste0(p, "-PCI95"),
        paste0(p, "-PIQR"),
        paste0(p, "-PMaxDiff"),
        paste0(p, "-p.adj"),
        paste0(p, "-PadjNumSignificant"),
        paste0(p, "-PadjMedian"),
        paste0(p, "-PadjMean"),
        paste0(p, "-PadjSD"),
        paste0(p, "-PadjCI95"),
        paste0(p, "-PadjIQR"),
        paste0(p, "-PadjMaxDiff")
      )
    }
  }

  
  # build desired order for all prefixes
  ordered_main <- unlist(lapply(prefixes, order_prefix))
  
  # keep only columns that actually exist (important!)
  ordered_main <- ordered_main[ordered_main %in% cols]
  
  # everything else (overlaps, deltas, etc.)
  remaining <- setdiff(cols, ordered_main)
  
  # final order
  #combinedData <- combinedData[, c(ordered_main, remaining)]
  combinedData <- combinedData[, c(ordered_main)]
}


# --- Add column which includes the main analysis as an alternate

prefixes <- unique(sub("-.*", "", names(combinedData)))

for (prefix in prefixes) {
  
  p_col <- paste0(prefix, "-P")
  padj_col <- paste0(prefix, "-p.adj")
  
  pnum_col <- paste0(prefix, "-PNumSignificant")
  padjnum_col <- paste0(prefix, "-PadjNumSignificant")
  
  # skip if columns don't exist
  if (!all(c(p_col, padj_col, pnum_col, padjnum_col) %in% names(combinedData))) next
  
  # new column names
  pnum_new <- paste0(prefix, "-PNumSignificantInclusive")
  padjnum_new <- paste0(prefix, "-PadjNumSignificantInclusive")
  
  # create new values
  combinedData[[pnum_new]] <- combinedData[[pnum_col]] + 
    ifelse(combinedData[[p_col]] < 0.05, 1, 0)
  
  combinedData[[padjnum_new]] <- combinedData[[padjnum_col]] + 
    ifelse(combinedData[[padj_col]] < 0.05, 1, 0)
  
  # move columns to correct position (after originals)
  p_index <- match(pnum_col, names(combinedData))
  padj_index <- match(padjnum_col, names(combinedData))
  
  # reorder for PNum
  combinedData <- combinedData[, append(
    names(combinedData)[-which(names(combinedData) == pnum_new)],
    pnum_new,
    after = p_index
  )]
  
  # reorder for PadjNum
  combinedData <- combinedData[, append(
    names(combinedData)[-which(names(combinedData) == padjnum_new)],
    padjnum_new,
    after = padj_index
  )]
}



for (prefix in prefixes) {
  
  sigCol <- paste0(prefix, "-significant")

  padjnum_col <- paste0(prefix, "-PadjNumSignificant")
  
  # skip if columns don't exist
  if (!all(c(sigCol,  padjnum_col) %in% names(combinedData))) next
  
  # new column names
  newCol <- paste0(prefix, "-significantRobust")

  
  # create new values
  combinedData[[newCol]] = combinedData[[sigCol]] & combinedData[[padjnum_col]] > alternateCutoff
    
  
  # move columns to correct position (after originals)
  sigCol_index <- match(sigCol, names(combinedData))
  
  # reorder for PNum
  combinedData <- combinedData[, append(
    names(combinedData)[-which(names(combinedData) == newCol)],
    newCol,
    after = sigCol_index
  )]
  
}



if(saveData){
  
  if(usingGene){
    combinedDataFilename = paste0(outputFolderName, filePrefix, "combinedGeneResultsWithAlternates")
    write.csv(combinedData, paste0(combinedDataFilename, ".csv"))
    saveRDS(combinedData, paste0(combinedDataFilename, ".rds"))    
  }else{
    combinedGODataFilename = paste0(outputFolderName, filePrefix, "combinedGOResultsWithAlternates-", geneSet)
    write.csv(combinedData, paste0(combinedGODataFilename, ".csv"))
    saveRDS(combinedData, paste0(combinedGODataFilename, ".rds"))
  }

}






