rerViolinPlot = function(mainTrees, RERObject, pathsObject, phenotypeSet, geneOfInterest, colorScale = c("red", "blue", "black", "green", "yellow"), correlationFile = NULL, legend = F){
  source("Src/Reu/RERConvergeFunctions.R")
  rerTree = returnRersAsTree(mainTrees, RERObject, geneOfInterest, pathsObject, plot = F)
  relativeRate = rerTree$edge.length
  geneTree = mainTrees$trees[[geneOfInterest]]
  
  edgeIndexMap = edgeIndexRelativeMaster(geneTree, mainTrees$masterTree)
  relevantBranches = mainTrees$matIndex[edgeIndexMap[,c(2,1)]]
  relevantPath = pathsObject[relevantBranches]
  
  phenotypeVector = c(rep(NA, length(relativeRate)))
  for(i in unique(relevantPath)){
    if(!is.na(i)){
      phenotypeVector[relevantPath == i] = phenotypeSet[i]
    }
  }
  
  #Adding Variables with different caps for better Aesthetics 
  Phenotype = phenotypeVector
  RelativeRate = relativeRate
  
  rateWithPhenotype = data.frame(RelativeRate, Phenotype)
  rateWithPhenotype = rateWithPhenotype[!is.na(rateWithPhenotype$Phenotype),]
  
  plot = ggplot(rateWithPhenotype, aes(x=Phenotype, y=RelativeRate, col=Phenotype)) + 
    geom_violin(adjust=1/3, show.legend = legend) +
    geom_jitter(position=position_jitter(0.2), show.legend = legend) +
    theme_classic() +
    scale_color_manual(values=colorScale) +
    theme(text = element_text(size = 20))+
    ggtitle(geneOfInterest)
  
  if(!is.null(correlationFile)){
    rowNumber = grep(paste("\\<",geneOfInterest,"\\>", sep=""), rownames(correlationFile))
    if("permPValue" %in% colnames(correlationFile)){
      plot = plot + ggtitle(geneOfInterest, subtitle = paste("p.adj = ", signif(correlationFile$p.adj[rowNumber], 4), "    Permulated P = ", signif(correlationFile$permPValue[rowNumber],4))) 
    }else{
      plot = plot + ggtitle(geneOfInterest, subtitle = paste("p.adj = ", signif(correlationFile$p.adj[rowNumber], 4))) 
    }
  }
  plot
  return(plot)
}

rerViolinPlotOld = function(mainTrees, RERObject, phenotypeTree, foregroundSpecies, geneOfInterest, foregroundName = "Foreground", backgroundName = "Background", foregroundColor = "lightsalmon", backgroundColor= "darkgreen", correlationFile = NULL){
  source("Src/Reu/RERConvergeFunctions.R")
  rerTree = returnRersAsTree(mainTrees, RERObject, geneOfInterest, foregroundSpecies, plot = F)
  relativeRate = rerTree$edge.length
  geneTree = mainTrees$trees[[geneOfInterest]]
  phenotypePath = tree2Paths(phenotypeTree, mainTrees)
  
  edgeIndexMap = edgeIndexRelativeMaster(geneTree, mainTrees$masterTree)
  relevantBranches = mainTrees$matIndex[edgeIndexMap[,c(2,1)]]
  relevantPath = phenotypePath[relevantBranches]
  
  phenotypeVector = c(rep(NA, length(relativeRate)))
  phenotypeVector[relevantPath == 1] = foregroundName
  phenotypeVector[relevantPath == 0] = backgroundName
  
  #Adding Variables with different caps for better Aesthetics 
  Phenotype = phenotypeVector
  RelativeRate = relativeRate
  
  rateWithPhenotype = data.frame(RelativeRate, Phenotype)
  rateWithPhenotype = rateWithPhenotype[!is.na(rateWithPhenotype$Phenotype),]
  
  plot = ggplot(rateWithPhenotype, aes(x=Phenotype, y=RelativeRate, col=Phenotype)) + 
    geom_violin(adjust=1/3) +
    geom_jitter(position=position_jitter(0.2)) +
    theme_classic() +
    scale_color_manual(values=c(foregroundColor,backgroundColor)) +
    theme(text = element_text(size = 20))+
    ggtitle(geneOfInterest)
  
  if(!is.null(correlationFile)){
    rowNumber = grep(paste("\\<",geneOfInterest,"\\>", sep=""), rownames(correlationFile))
    if("permPValue" %in% colnames(correlationFile)){
      plot = plot + ggtitle(geneOfInterest, subtitle = paste("p.adj = ", signif(correlationFile$p.adj[rowNumber], 4), "    Permulated P = ", signif(correlationFile$permPValue[rowNumber],4))) 
    }else{
      plot = plot + ggtitle(geneOfInterest, subtitle = paste("p.adj = ", signif(correlationFile$p.adj[rowNumber], 4))) 
    }
  }
  plot
  return(plot)
}


