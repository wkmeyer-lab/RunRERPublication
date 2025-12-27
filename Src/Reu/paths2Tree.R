paths2Tree =function (treesObj, phenv, index, 
          rer.cex = 0.7, tip.cex = 0.7, nalab = "NA", plot = T) 
{
  trgene <- treesObj$trees[[index]]
  trgene$edge.length <- rep(2, nrow(trgene$edge))
  ee = edgeIndexRelativeMaster(trgene, treesObj$masterTree)
  ii = treesObj$matIndex[ee[, c(2, 1)]]
  rertree = phenv[ii]
  rertree[is.nan(rertree)] = NA
  if (plot) {
    par(mar = c(1, 1, 1, 0))
    edgcols <- rep("black", nrow(trgene$edge))
    edgwds <- rep(1, nrow(trgene$edge))
    if (!is.null(phenv)) {
      edgcols <- rep("gray", nrow(trgene$edge))
      edgwds <- rep(1, nrow(trgene$edge))
      for(i in 1:length(unique(phenv[ii]))){
        edgcols[phenv[ii] == i] <- palette()[i]
      }
    }
    plot.phylo(trgene, font = 2, edge.color = edgcols, edge.width = edgwds, 
               cex = tip.cex)
    rerlab <- round(rertree, 3)
    rerlab[is.na(rerlab)] <- nalab
    edgelabels(rerlab, bg = NULL, adj = c(0.5, 0.9), col = edgcols, 
               frame = "none", cex = rer.cex, font = 2)
  }
  trgene$edge.length <- rertree
  return(trgene)
}


#' @keywords  internal
edgeIndexRelativeMaster=function(tree, masterTree){
  map=matchAllNodes(tree,masterTree)
  newedge=tree$edge
  newedge[,1]=map[newedge[,1],2]
  newedge[,2]=map[newedge[,2],2]
  newedge
}
#' @keywords  internal
matchAllNodes=function(tree1, tree2){
  map=matchNodesInject(tree1,tree2)
  map=map[order(map[,1]),]
  map
}
#' @keywords  internal
matchNodesInject=function (tr1, tr2){
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