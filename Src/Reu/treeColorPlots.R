treeColorPlot = function(trgene, rer.cex = 0.7, tip.cex = 0.7, nalab = "NA"){
  par(mar = c(1, 1, 1, 0))
  edgcols <- rep("black", nrow(trgene$edge))
  edgwds <- rep(1, nrow(trgene$edge))
  edgcols <- rep("gray", nrow(trgene$edge))
  edgwds <- rep(1, nrow(trgene$edge))
  for(i in 1:length(unique(trgene$edge.length))){
    edgcols[trgene$edge.length == i] <- palette()[i]
  }
  plot.phylo(trgene, font = 2, edge.color = edgcols, edge.width = edgwds, 
             cex = tip.cex)
  rerlab <- round(trgene$edge.length, 3)
  rerlab[is.na(rerlab)] <- nalab
  edgelabels(rerlab, bg = NULL, adj = c(0.5, 0.9), col = edgcols, 
             frame = "none", cex = rer.cex, font = 2)
}

treeColorByLabel = function(phenMasterTree){
  startingPalette = palette()
  bufferedPalette = append("white", startingPalette)
  palette(bufferedPalette)
  par(mar = c(1, 1, 1, 0))
  edgcols <- rep("gray", nrow(phenMasterTree$edge))
  tipcols <- rep("gray", length(phenMasterTree$tip.label))
  tipEnds = substr(phenMasterTree$tip.label, nchar(phenMasterTree$tip.label)-2, nchar(phenMasterTree$tip.label))
  for(i in 2:length(unique(phenMasterTree$node.label))){
    nodeLabel = unique(phenMasterTree$node.label)[i]
    #print(nodeLabel)
    parentNodes = which(phenMasterTree$node.label == nodeLabel)
    parentNodes = parentNodes + length(phenMasterTree$tip.label)
    #print(parentNodes)
    parentTips = which(tipEnds == nodeLabel)
    tipcols[parentTips] <- palette()[i]
    #print(parentTips)
    parentPostions = append(parentNodes, parentTips)
    parentEdges = which(phenMasterTree$edge[,2] %in% parentPostions)
    #print(parentEdges)
    edgcols[parentEdges] <- palette()[i]
    #print(edgcols)
  }
  palette(startingPalette)
  plot = plot.phylo(phenMasterTree, font = 2, edge.color = edgcols, cex = 0.7, tip.color = tipcols)
  
}
