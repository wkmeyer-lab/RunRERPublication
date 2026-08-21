library(scales)
source("Src/Reu/ZonomNameConvertVectorCommon.R")
displayCategoricalRERTree = function(treesObj, rermat, index, phenv = NULL, subsetTree = T, equalLengths = T, minWidth = 1, maxWidth = 6, tipCol = "tipColumn", annotLocation = "Data/mergedData.csv"){
  treesObj$trees[[index]]$tip.label = ZonomNameConvertVectorCommon(treesObj$trees[[index]]$tip.label, tipColumn = tipCol, annotationLocation = annotLocation)
  treesObj$masterTree$tip.label = ZonomNameConvertVectorCommon(treesObj$masterTree$tip.label, tipColumn = tipCol, annotationLocation = annotLocation)
  colnames(rermat) = ZonomNameConvertVectorCommon(colnames(rermat), tipColumn = tipCol)
  returnRersAsTreeNew(treesObj, rermat, index, phenv, 0.7, 0.7, 'NA', T, subsetTree, equalLengths, minWidth, maxWidth)
}



returnRersAsTreeNew <- function(treesObj, rermat, index, phenv = NULL, rer.cex = 0.7,
                                tip.cex = 0.7, nalab = 'NA', plot = T, subsetTree = F, equalLengths=T, minWidth=1, maxWidth=6){
  trgene <- treesObj$trees[[index]]
  if (subsetTree == TRUE) {
    #Use phenv to identify non-NA paths and then subset gene tree by species
    subsp <- unique(colnames(rermat)[which(!is.na(phenv))])
    trgene <- pruneTree(trgene, subsp) #drops only species not in subsp
  }
  if(equalLengths){trgene$edge.length <- rep(2,nrow(trgene$edge))}
  ee=edgeIndexRelativeMaster(trgene, treesObj$masterTree)
  ii= treesObj$matIndex[ee[, c(2,1)]]
  rertree=rermat[index,ii]
  rertree[is.nan(rertree)]=NA #replace NaNs from C functions
  relativeRER = abs(scale(rertree))
  if(any(relativeRER > 3)){
    message("WARNING: Potential misannotated gene!")
    for(i in which(relativeRER > 3)){
      if(!is.na(rownames(relativeRER)[i])){
        message(paste("Potential Misannotated gene at branch:", i ))
        message(paste("Species of branch=:", rownames(relativeRER)[i]))
        if(plot){
          message("Branch highlighted in PINK.")
        }
      }
    }
  }
  
  if (plot) {
    clampedRelativeRER = relativeRER
    clampedRelativeRER[order(relativeRER, decreasing = TRUE)[1:3]] <- sort(relativeRER, decreasing = TRUE)[4] #Set the top three high RERs to the fourth highest RER; this matches the calculation ignoring top values and prevents mis-IDed genes from throwing off the scale. 
    rerWidth = rescale(clampedRelativeRER, c(minWidth,maxWidth))
    rerWidth[is.na(rerWidth)] = 1
    
    par(mar = c(1,0,1,0))
    edgcols <- rep('black', nrow(trgene$edge))
    edgwds <- rep(1, nrow(trgene$edge))
    if(!is.null(phenv)){
      edgcols <- rep('black', nrow(trgene$edge))
      edgwds <- rerWidth
      if(length(unique(pathsObject) < length(palette()))){ # add a catch for continuous phenotypes and not run it in that case
        for(j in unique(pathsObject)[!is.na(unique(pathsObject))]){
          edgcols[phenv[ii]==j] <- palette()[j]
        }      
      }
    }
    if(any(relativeRER > 3)){
      for(i in which(relativeRER > 3)){
        if(!is.na(rownames(relativeRER)[i])){
          edgcols[i] = 'hotpink'
        }
      }
    }
    plot.phylo(trgene, font = 2, edge.color = edgcols, edge.width = edgwds, cex = tip.cex)
    rerlab <- round(rertree,3)
    rerlab[is.na(rerlab)] <- nalab
    if(!is.null(phenv)){ #reset the color so that the RER labels still have the correct palette color even with the pink warning
      for(j in unique(pathsObject)[!is.na(unique(pathsObject))]){
        edgcols[phenv[ii]==j] <- palette()[j]
      }
    }
    edgelabels(rerlab, bg = NULL, adj = c(0.5,0.9), col = edgcols, frame = 'none',cex = rer.cex, font =2)
  }
  
  
  trgene$edge.length <- rertree
  return(trgene)
}

{
  edgeIndexRelativeMaster= function(tree, masterTree){
    map=matchAllNodes(tree,masterTree)
    newedge=tree$edge
    newedge[,1]=map[newedge[,1],2]
    newedge[,2]=map[newedge[,2],2]
    newedge
  }
  matchAllNodes = function(tree1, tree2){
    map=matchNodesInject(tree1,tree2)
    map=map[order(map[,1]),]
    map
  }
  matchNodesInject = function (tr1, tr2){
    if(length(tmpsp<-setdiff(tr1$tip.label, tr2$tip.label))>0){
      #stop(paste(paste(tmpsp, ","), "in tree1 do not exist in tree2"))
      stop(c("The following species in tree1 do not exist in tree2: ",paste(tmpsp, ", ")))
    }
    commontiplabels <- intersect(tr1$tip,tr2$tip)
    if(RF.dist(pruneTree(tr1,commontiplabels),pruneTree(tr2,commontiplabels))>0){
      stop("Discordant tree topology detected - gene/trait tree and treesObj$masterTree have irreconcilable topologies")
    }
    #if(RF.dist(tr1,tr2)>0){
    #  stop("Discordant tree topology detected - trait tree and treesObj$masterTree have irreconcilable topologies")
    #}
    
    toRm=setdiff(tr2$tip.label, tr1$tip.label)
    desc.tr1 <- lapply(1:tr1$Nnode + length(tr1$tip), function(x) extract.clade(tr1,
                                                                                x)$tip.label)
    names(desc.tr1) <- 1:tr1$Nnode + length(tr1$tip)
    desc.tr2 <- lapply(1:tr2$Nnode + length(tr2$tip), function(x) extract.clade(tr2,
                                                                                x)$tip.label)
    names(desc.tr2) <- 1:tr2$Nnode + length(tr2$tip)
    Nodes <- matrix(NA, length(desc.tr1), 2, dimnames = list(NULL,
                                                             c("tr1", "tr2")))
    for (i in 1:length(desc.tr1)) {
      Nodes[i, 1] <- as.numeric(names(desc.tr1)[i])
      for (j in 1:length(desc.tr2)) if (all(desc.tr1[[i]] %in%
                                            desc.tr2[[j]]))
        Nodes[i, 2] <- as.numeric(names(desc.tr2)[j])
    }
    
    iim=match(tr1$tip.label, tr2$tip.label)
    Nodes=rbind(cbind(1:length(tr1$tip.label),iim),Nodes)
    if(any(table(Nodes[,2])>1)){
      stop("Incorrect pseudorooting detected - use fixPseudoroot() function to correct trait tree topology")
    }
    
    Nodes
  }
  
}