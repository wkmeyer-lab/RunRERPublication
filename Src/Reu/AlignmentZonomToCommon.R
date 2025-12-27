#This function converts Zoonomia names to common names on a FASTA file, using the manual annotations spreadsheet. 
#If provided a file prefix, it will filter the alignment to the species used by that prefix, and sort the species by foreground/background
alignmentZonomToCommon = function(inFile, outfile, prefix = NULL, binary = T, Format = "fasta"){
  library(RERconverge)
  Alignment = read.phyDat(inFile, format = Format)
  if(!is.null(prefix)){
    outputFolderName = paste("Output/",prefix,"/", sep = "")
    speciesFilter = readRDS(paste(outputFolderName, prefix, "SpeciesFilter.rds", sep=""))
    foregroundList = readRDS(paste(outputFolderName, prefix, "BinaryTreeForegroundSpecies.rds", sep=""))
    AlignFiltered = Alignment[which(names(Alignment) %in% speciesFilter),]
    Alignment = AlignFiltered
    if(binary){
      AlignForeground = AlignFiltered[which(names(AlignFiltered) %in% foregroundList),]
      AlignBackground = AlignFiltered[-which(names(AlignFiltered) %in% foregroundList),]
      AlignSorted = c(AlignForeground, AlignBackground)
      Alignment = AlignSorted
    }
  }

  speciesNames = names(Alignment)
  source("ZonomNameConvertVector.R")
  commonNames = ZonomNameConvertVectorCommon(speciesNames)
  names(Alignment) = commonNames
  Alignment
  write.phyDat(Alignment, file=outfile, format = Format)
}