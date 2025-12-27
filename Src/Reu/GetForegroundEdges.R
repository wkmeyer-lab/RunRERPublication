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