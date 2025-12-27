mainPrefix = "NewHiller4Phen"
subdirectory = "_Omnivore-Carnivore"
subdirectory = "_Omnivore-Insectivore"
subdirectory = "_Omnivore-Herbivore"
filetag = paste(mainPrefix, subdirectory, sep="")
outputFolder = paste("Output", mainPrefix, subdirectory, "", sep="/")

mainPrefix = "CVHRemake"
subdirectory = ""
filetag = paste(mainPrefix, sep="")
outputFolder = paste("Output", mainPrefix, "", sep="/")


correlationsfile = paste(outputFolder, mainPrefix, subdirectory, "CorrelationFile.rds", sep="")
correlations = readRDS(correlationsfile)
rownames(correlations)
correlations$Gene = rownames(correlations)

sorteationsSorted = correlations[order(correlations$Rho, decreasing = T),]
sorteationsSorted$Rank = nrow(sorteationsSorted):1
rnkFile = sorteationsSorted[,c(4,1)]

rnkFileName = correlationsfile = paste(outputFolder, mainPrefix, subdirectory, "RankFile.rnk", sep="")
write.table(rnkFile, rnkFileName, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

outputFolder = paste("Output", mainPrefix, subdirectory, "", sep="/")
enrichmentFile = paste(outputFolder, mainPrefix, subdirectory, "Enrichment-GO_Biological_Process_2023.rds", sep="")
enrichment = readRDS(enrichmentFile)

enrichmentReorder = enrichment$GO_Biological_Process_2023
enrichmentReorder$phenoType = sign(enrichmentReorder$stat)
enrichmentReorder$GO.ID = rownames(enrichmentReorder)
enrichmentReorder$Description = rownames(enrichmentReorder)
enrichmentReorder$Description = sub("\\(.*", "",enrichmentReorder$Description)
enrichmentReorder$Empty = rep("")

enrichmentReorderClean = enrichmentReorder[,c(7,8,3,9,6)]
names(enrichmentReorderClean)[3] = "p.Val"

cytoEnrichFile = paste(outputFolder, mainPrefix, subdirectory, "CytoscapeEnrichmentsFile.tsv", sep="")
write.table(enrichmentReorderClean, cytoEnrichFile, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

enrichmentReorderCleanSignificant = enrichmentReorderClean[enrichmentReorderClean$p.Val <0.1,]
enrichmentReorderCleanSignificant


makeCytoEnrichTable = function( subdirectory, mainprefix = mainPrefix, save = F){
  outputFolder = paste("Output", mainprefix, subdirectory, "", sep="/")
  enrichmentFile = paste(outputFolder, mainprefix, subdirectory, "Enrichment-GO_Biological_Process_2023.rds", sep="")
  enrichment = readRDS(enrichmentFile)
  
  enrichmentReorder = enrichment$GO_Biological_Process_2023
  enrichmentReorder$phenoType = sign(enrichmentReorder$stat)
  enrichmentReorder$GO.ID = rownames(enrichmentReorder)
  enrichmentReorder$Description = rownames(enrichmentReorder)
  enrichmentReorder$Description = sub("\\(.*", "",enrichmentReorder$Description)
  enrichmentReorder$Empty = rep("")
  
  enrichmentReorderClean = enrichmentReorder[,c(7,8,3,9,6)]
  names(enrichmentReorderClean)[3] = "p.Val"
  enrichmentReorderClean
  
  if(save){
    cytoEnrichFile = paste(outputFolder, mainPrefix, subdirectory, "CytoscapeEnrichmentsFile.tsv", sep="")
    write.table(enrichmentReorderClean, cytoEnrichFile, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    
    enrichmentReorderCleanSignificant = enrichmentReorderClean[enrichmentReorderClean$p.Val <0.1,]
    enrichmentReorderCleanSignificant
    
    cytoEnrichSignificantFile = paste(outputFolder, mainPrefix, subdirectory, "CytoscapeSignificantEnrichmentsFile.tsv", sep="")
    write.table(enrichmentReorderCleanSignificant, cytoEnrichSignificantFile, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    
    
  }
  enrichmentReorderClean
}


OvCEnrichments = makeCytoEnrichTable("_Omnivore-Carnivore", save = T)
OvCEnrichmentsSignificant = OvCEnrichments[OvCEnrichments$p.Val <0.1,]

OvIEnrichments = makeCytoEnrichTable("_Omnivore-Insectivore", save = T)
OvIEnrichmentsSignificant = OvIEnrichments[OvIEnrichments$p.Val <0.1,]

OvHEnrichments = makeCytoEnrichTable("_Omnivore-Herbivore", save = T)
OvHEnrichmentsSignificant = OvHEnrichments[OvHEnrichments$p.Val <0.1,]



OvCEnrichmentsSignificant$phenoType[OvCEnrichmentsSignificant$phenoType == 1]= "1"
OvCEnrichmentsSignificant$phenoType[OvCEnrichmentsSignificant$phenoType == -1]= "-1"


OvIEnrichmentsSignificant$phenoType[OvIEnrichmentsSignificant$phenoType == 1]= "2"
OvIEnrichmentsSignificant$phenoType[OvIEnrichmentsSignificant$phenoType == -1]= "-2"


OvHEnrichmentsSignificant$phenoType[OvHEnrichmentsSignificant$phenoType == 1]= "1"
OvHEnrichmentsSignificant$phenoType[OvHEnrichmentsSignificant$phenoType == -1]= "-1"

allOmnivoreEnrichments = rbind(OvCEnrichmentsSignificant, OvIEnrichmentsSignificant)
allOmnivoreEnrichments = rbind(allOmnivoreEnrichments, OvHEnrichmentsSignificant)

AllOmnicytoEnrichFile = paste("Output/NewHiller4Phen/Overall/AllOmnivoreCytoscapeEnrichmentsFile.tsv", sep="")
write.table(allOmnivoreEnrichments, AllOmnicytoEnrichFile, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

