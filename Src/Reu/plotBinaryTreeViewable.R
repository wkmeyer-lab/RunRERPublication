
plotBinaryTreeViewable = function(inTree){
  treeDisplayable = inTree
  treeDisplayable$edge.length = replace(treeDisplayable$edge.length, treeDisplayable$edge.length==0, 0.5)
  treeDisplayable$edge.length = replace(treeDisplayable$edge.length, treeDisplayable$edge.length==1, 4)
  
  plotTreeHighlightBranches(treeDisplayable, hlspecies=which(inTree$edge.length>=1), hlcols="blue",)
  return()
}
