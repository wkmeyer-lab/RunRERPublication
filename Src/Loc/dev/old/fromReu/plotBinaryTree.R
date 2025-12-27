#This script takes an input multiphylo and gene name, and outputs 1) a tree based on the master-tree branch lengths with the gene tree topology and 2) the gene tree. 
#If an RER object is provided, it will instead trim to only tips which have RER values. 
#If provided a foreground vector, it can color the foreground and background different colors. 
#This script is no longer dependent on RERConverge
#This script will by default convert names from zoonomia names to common names. This requires "Data/manualAnnotationsSheet.csv". This can be toggled off using convertNames = F=
#This script can also make a plot of overall genome length vs gene length by toggling correlationPlot = T.

source("Src/Reu/ZoonomTreeNameToCommon.R")

plotBinaryTree = function(mainTrees, binaryTree, foregroundVector = NULL, fgcols = "blue", bgcolor = "black", convertNames = T, mainTitle = NULL, tipColumn = "tipName"){
  masterTree = mainTrees$masterTree
  namesToKeep = binaryTree$tip.label
  prunedMaster = drop.tip(masterTree, which(!masterTree$tip.label %in% namesToKeep))
  
  if(convertNames){
    if(file.exists("Src/Reu/ZonomNameConvertVector.R")){source("Src/Reu/ZonomNameConvertVector.R")}
    if(file.exists("Src/Reu/ZoonomTreeNameToCommon.R")){source("Src/Reu/ZoonomTreeNameToCommon.R")}
    commonMaster = ZoonomTreeNameToCommon(prunedMaster, plot = F, tipCol = tipColumn)
    plotMaster = commonMaster
  }else{
    plotMaster = prunedMaster
  }
  par(mfrow = c(1,1))
  
  binaryForegroundEdges = which(binaryTree$edge.length ==1 )
  binaryForegroundTips = NA
  for(i in binaryForegroundEdges){
    endNode = binaryTree$edge[i,][[2]]
    if(endNode <= length(binaryTree$tip.label)){
      binaryForegroundTips = append(binaryForegroundTips, endNode)
      }
  }

  edgeColors = rep(bgcolor, length(plotMaster$edge))
  edgeColors[binaryForegroundEdges] = fgcols
  
  tipColors = rep(bgcolor, length(plotMaster$tip.label))
  tipColors[binaryForegroundTips ] = fgcols
  
  plotobj = plot.phylo(plotMaster, main = mainTitle, edge.color = edgeColors, tip.color = tipColors, edge.width = 2, cex = 0.8)
}




# this script outputs the foreground edges for an "all" clades foreground given a tree and a foreground species vector. 

getForegroundEdges = function(inputTree, foregroundVector, plot = F){
  startnodes = min(inputTree$edge[,1]):max(inputTree$edge[,1])
  edges = inputTree$edge
  foregroundTips = which(inputTree$tip.label %in% foregroundVector)
  backgroundTips = which(!1:length(inputTree$tip.label) %in% foregroundTips)
  unsureNodes = startnodes
  for(i in 1:10){
    for(i in unsureNodes){
      endNodes = edges[which(edges[,1] == i),2]
      if(all(endNodes %in% foregroundTips)){
        foregroundTips = append(foregroundTips, i)
        unsureNodes = unsureNodes[unsureNodes != i]
      }else if(any(endNodes %in% backgroundTips)){
        backgroundTips = append(backgroundTips, i)
        unsureNodes = unsureNodes[unsureNodes != i]
      }else{
        unsureNodes = append(unsureNodes, i)
      }
      #message(endNodes)
    }
  }
  foregroundEdges = which(edges[,2] %in% foregroundTips)
  if(plot){
    plotTreeHighlightBranches(inputTree, hlspecies = foregroundEdges, hlcols = "blue")
  }
  return(foregroundEdges)
}