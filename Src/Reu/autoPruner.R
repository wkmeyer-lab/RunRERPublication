source("Src/Reu/ZoonomTreeNameToCommon.R")



autopruner= function(masterTree, dropPercent = NA, dropValue = 0.01, tipsToKeep = NA, returnEdgeTable = F, procedurePlot = F, nameConversionColumn = NA, nameConversionData = "Data/mergedData.csv", preDroppedTips = NA, originalTree = NA, mainTreesObject = mainTrees){
  if(all(is.na(originalTree))){originalTree = masterTree}
  
  report= mainTreesObject$report
  speciesGeneNumber = colSums(report)
  
  # -- determine the branch length cutoff -- 
  if(is.character(dropPercent)){dropPercent = as.numeric(dropPercent)}
  if(is.character(dropValue)){dropValue = as.numeric(dropValue)}
  if(!is.na(dropPercent)){
    if(dropPercent > 1){dropPercent = dropPercent/100}
    lengthCutoff = quantile(masterTree$edge.length, dropPercent)  
  }else{
    lengthCutoff = dropValue
  }
  print(lengthCutoff)
  
  # -- Make a table of the terminal edges-- 
  trimmedTree = masterTree 
  terminalEdges = which(trimmedTree$edge[,2] <= length(trimmedTree$tip.label))
  terminalEdges = which(trimmedTree$edge[,2] <= length(trimmedTree$tip.label))
  terminalEdgeTable = data.frame(edgeNumber = terminalEdges, edgeLength = trimmedTree$edge.length[terminalEdges], edgeStart = trimmedTree$edge[terminalEdges,1], edgeEnd = trimmedTree$edge[terminalEdges,2])
  terminalEdgeTable$edgeTip = trimmedTree$tip.label[terminalEdgeTable$edgeEnd]
  terminalEdgeTable = terminalEdgeTable[order(terminalEdgeTable$edgeLength),]
  if(returnEdgeTable){
    return(terminalEdgeTable)
    stop()
  }
  
  #message("Species in order of shortest tips:")
  #message(paste(terminalEdgeTable$edgeTip, collapse="\n"))
  
  droppedTips = vector()
  skippedTips = vector()
  
  # --- drop tips until the branch lengths are within the cutoff --- 
  while(min(terminalEdgeTable$edgeLength) <= lengthCutoff){
    {
    selectingBranch = T; i = 1
    while(selectingBranch){
      shortestTerminalBranch = terminalEdgeTable[i,]
      if(shortestTerminalBranch$edgeTip %in% tipsToKeep){
        if(!shortestTerminalBranch$edgeTip %in% skippedTips){
          message(paste("Skipping", shortestTerminalBranch$edgeTip))
          skippedTips = append(skippedTips, shortestTerminalBranch$edgeTip)
        }
        i = i+1
      }else{
        message(paste("Pruning", shortestTerminalBranch$edgeTip, " -- length:", shortestTerminalBranch$edgeLength))
        i = 1
        tipDropped = shortestTerminalBranch$edgeTip
        names(tipDropped) = shortestTerminalBranch$edgeLength
        droppedTips = append(droppedTips, tipDropped)
        selectingBranch = F
      }
    }
    trimmedTree = drop.tip(trimmedTree, shortestTerminalBranch$edgeTip)
    
    terminalEdges = which(trimmedTree$edge[,2] <= length(trimmedTree$tip.label))
    terminalEdgeTable = data.frame(edgeNumber = terminalEdges, edgeLength = trimmedTree$edge.length[terminalEdges], edgeStart = trimmedTree$edge[terminalEdges,1], edgeEnd = trimmedTree$edge[terminalEdges,2])
    terminalEdgeTable$edgeTip = trimmedTree$tip.label[terminalEdgeTable$edgeEnd]
    if(!length(which(terminalEdgeTable$edgeTip %in% skippedTips))==0){
      terminalEdgeTable = terminalEdgeTable[-which(terminalEdgeTable$edgeTip %in% skippedTips),] #This prevents the while loop from continuing due to a skipped tip 
    }
    terminalEdgeTable = terminalEdgeTable[order(terminalEdgeTable$edgeLength),]
    
    if(procedurePlot){
      tipDropPlot(nameConversionColumn, nameConversionData, originalTree, trimmedTree, droppedTips, skippedTips, preDroppedTips, message = F)
    }
    } 
    min(trimmedTree$edge.length[terminalEdges])
  }

  
  message(paste("Skipped tips:"))
  message(paste(skippedTips, collaspe = ", ", sep=""))
  message(paste("Dropped tips:"))
  message((paste(droppedTips, collaspe = ", ", sep="")))  
  
  if(length(droppedTips > 0)){
    if(!is.na(nameConversionColumn)){
      droppedCommonTips = ZonomNameConvertVectorCommon(droppedTips, annotationLocation = nameConversionData, tipColumn = nameConversionColumn)
      droppedTable = data.frame(droppedTips, droppedCommonTips, names(droppedTips), speciesGeneNumber[match(droppedTips,  names(speciesGeneNumber))], match(droppedTips,  originalTree$tip.label))
      names(droppedTable) = c("Dropped Tip", "Common Name", "Branch Length", "Num Genes", "Tip Position")
      tipDropPlot(nameConversionColumn, nameConversionData, originalTree, trimmedTree, droppedTips, skippedTips, preDroppedTips, message = T)
    }else{
      droppedTable = data.frame(droppedTips, names(droppedTips), speciesGeneNumber[match(droppedTips,  names(speciesGeneNumber))], match(droppedTips,  originalTree$tip.label))
      names(droppedTable) = c("Dropped Tip", "Branch Length", "Num Genes", "Tip Position")
    }
    droppedTable = droppedTable[order(droppedTable$`Tip Position`),]
    print(droppedTable)
    
    par(mfrow = c(1,2))
    plotTreeHighlightBranches(originalTree, hlspecies = droppedTips, hlcols = "red")
    plot.phylo(trimmedTree)
    par(mfrow = c(1,1))
  }
  
  droppedTips = append(droppedTips, preDroppedTips)
  droppedTips <<- droppedTips
  trimmedTree
}

tipDropPlot = function(nameConversionColumn, nameConversionData, masterTree, trimmedTree, droppedTips, skippedTips, preDroppedTips, message = F){
  
  if(!all(is.na(preDroppedTips))){
    droppedTips = append(droppedTips, preDroppedTips)
  }
  
  if(!is.na(nameConversionColumn)){
    plotMasterTree = ZoonomTreeNameToCommon(masterTree, manualAnnotLocation = nameConversionData, tipCol = nameConversionColumn, plot = F)
    plotTrimmedTree = ZoonomTreeNameToCommon(trimmedTree, manualAnnotLocation = nameConversionData, tipCol = nameConversionColumn, plot = F)
    plotDroppedTips = ZonomNameConvertVectorCommon(droppedTips, annotationLocation = nameConversionData, tipCol = nameConversionColumn)
    plotSkippedTips = ZonomNameConvertVectorCommon(skippedTips, annotationLocation = nameConversionData, tipCol = nameConversionColumn)
  }else{
    plotMasterTree = masterTree
    plotTrimmedTree = trimmedTree
    plotDroppedTips = droppedTips
    plotSkippedTips = skippedTips
  }
  
  par(mfrow = c(1,2))
  plotTreeHighlightBranches(plotMasterTree, hlspecies = plotDroppedTips, hlcols = "red")
  plot.phylo(plotTrimmedTree)
  par(mfrow = c(1,1))
  
  if(message){
    message(paste("Skipped tips:"))
    message(paste(plotSkippedTips, collaspe = ", ", sep=""))
    message(paste("Dropped tips:"))
    message((paste(plotDroppedTips, collaspe = ", ", sep="")))
  }
}

dropFewGeneSpecies = function(mainTrees, masterTree = NA, cutoff = -3, nameConversionColumn = NA, nameConversionData = "Data/mergedData.csv"){
  reportTable= mainTrees$report
  
  speciesGeneNumber = colSums(mainTrees$report)
  
  zScores = scale(speciesGeneNumber)[,1]
  fewGeneSp = zScores[which(zScores < -3)]
  fewGeneSpNumber = speciesGeneNumber[names(speciesGeneNumber) %in% names(fewGeneSp)]
  commonFewGeneSp = ZonomNameConvertVectorCommon(names(fewGeneSp), annotationLocation = nameConversionData, tipCol = nameConversionColumn)
  
  fewGeneTable = data.frame(commonFewGeneSp, names(fewGeneSp), fewGeneSp, fewGeneSpNumber)
  names(fewGeneTable) = c("CommonName", "TipName", "Z-score", "Num_of_Genes")
  
  droppedTips = fewGeneTable
  if(!all(is.na(masterTree))){
    droppedTips = droppedTips[droppedTips$TipName %in% masterTree$tip.label,]
  }
  message("Dropping these tips due to few genes from species:")
  message(paste(names(droppedTips), collapse = "   "))
  for(i in 1:nrow(droppedTips)){
    message(paste(droppedTips[i,], collapse = "   "))
  }
  tipsToDrop = droppedTips$TipName
  
  if(length(tipsToDrop) == 0){tipsToDrop = NA }
  return(tipsToDrop)
}

