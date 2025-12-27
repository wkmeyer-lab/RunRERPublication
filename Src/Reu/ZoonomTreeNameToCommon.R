#USAGE: 
#Input a tree you want tip labels changed to common name 
#Chose if you want a plot 
#Specify if foreground tree (only used for plotting)
#Will output tree of same shape with tips renamed

#Will expect the manual annotations spreadsheet to be placed in "Data/manualAnnotationsSheet.csv" by default, if not there, specify location 

#For manual use: 
#treeToConvertLocation = "Data/CVHRemakeBinaryForegroundTree.rds"
#inputTree = readRDS(treeToConvertLocation)

source("Src/Reu/ZonomNameConvertVectorCommon.R")

ZoonomTreeNameToCommon = function(tree, plot = T, isForegroundTree = T, manualAnnotLocation = "Data/mergedData.csv", hlcol = "blue", bgcol = "black", fontSize = 0.8, scientific = F, scientificCol = "ScientificName", commonCol = "CommonName", tipCol = "tipName"){
  
  inputTree = tree
  tipNames = inputTree$tip.label
  
  tipNames = ZonomNameConvertVectorCommon(tipNames, annotationLocation = manualAnnotLocation, toScientific = scientific, scientificColumn = scientificCol, commonColumn = commonCol, tipColumn = tipCol)
  
  inputTree$tip.label = tipNames
  
  if(plot){
    if(isForegroundTree){
      readableTree = inputTree
      readableTree$edge.length[readableTree$edge.length == 0] = 1
      plotTreeHighlightBranches2(readableTree, hlspecies = which(inputTree$edge.length == 1), hlcols = hlcol, bgcol = bgcol, fontSize = fontSize)
    }else{
      plotTree(inputTree)
    }
  }
  return(inputTree)
}


plotTreeHighlightBranches2 =function (tree, outgroup = NULL, hlspecies, hlcols = NULL, main = "", 
                                      
                                      
                                      useGG = FALSE, bgcol = "black", fontSize =0.8) 
  
  
{
  
  
  if (is.null(hlcols)) {
    
    
    hlcols <- c(2:(length(hlspecies) + 1))
    
    
  }
  
  
  if (length(hlcols) < length(hlspecies)) {
    
    
    hlcols <- rep_len(hlcols, length(hlspecies))
    
    
  }
  
  
  if (!is.null(outgroup)) {
    
    
    outgroup <- outgroup[outgroup %in% tree$tip.label]
    
    
    if (length(outgroup) > 0) {
      
      
      if (is.numeric(hlspecies)) {
        
        
        dummytree <- tree
        
        
        dummytree$edge.length <- c(rep(1, nrow(dummytree$edge)))
        
        
        for (i in c(1:length(hlspecies))) {
          
          
          dummytree$edge.length[hlspecies] <- i + 1
          
          
        }
        
        
        dummyrooted <- root(dummytree, outgroup)
        
        
      }
      
      
      rooted <- root(tree, outgroup)
      
      
    }
    
    
    else {
      
      
      print("No members of requested outgroup found in tree; keeping unrooted.")
      
      
      rooted <- tree
      
      
      outgroup <- NULL
      
      
    }
    
    
  }
  
  
  else {
    
    
    rooted <- tree
    
    
  }
  
  
  if (useGG) {
    
    
    if (is.numeric(hlcols)) {
      
      
      hlcols <- rep_len("#0000ff", length(hlspecies))
      
      
    }
    
    
    hlspecies_named <- vector(mode = "character")
    
    
    if (is.numeric(hlspecies)) {
      
      
      for (i in 1:length(hlspecies)) {
        
        
        hlspecies_named[i] <- tree$tip.label[hlspecies[i]]
        
        
      }
      
      
    }
    
    
    else hlspecies_named <- hlspecies
    
    
    rooted2 <- rooted
    
    
    mm <- min(rooted2$edge.length[rooted2$edge.length > 0])
    
    
    rooted2$edge.length[rooted2$edge.length == 0] <- max(0.02, 
                                                         
                                                         
                                                         mm/20)
    
    
    tipCols <- vector(mode = "character")
    
    
    x <- 1
    
    
    for (i in 1:length(tree$tip.label)) {
      
      
      if (tree$tip.label[i] %in% hlspecies_named) {
        
        
        tipCols[i] <- hlcols[x]
        
        
        x <- x + 1
        
        
      }
      
      
      else tipCols[i] <- bgcol
      
      
    }
    
    
    nlabel <- rooted2$Nnode + length(rooted2$tip.label)
    
    
    edgeCols <- vector(mode = "character", length = nlabel)
    
    
    x <- 1
    
    
    for (i in 1:nlabel) {
      
      
      if (i < length(rooted2$tip.label)) {
        
        
        if (rooted2$tip.label[i] %in% hlspecies_named) {
          
          
          edgeCols[i] <- hlcols[x]
          
          
          x <- x + 1
          
          
        }
        
        
        else edgeCols[i] <- bgcol
        
        
      }
      
      
      else edgeCols[i] <- bgcol
      
      
    }
    
    
    plotobj = ggtree(rooted2, color = edgeCols)
    
    
    plotobj = plotobj + geom_tiplab(color = tipCols, geom = "text", 
                                    
                                    
                                    cex = 3) + labs(title = main)
    
    
    return(plotobj)
    
    
  }
  
  
  else {
    
    
    colMaster <- c(rep(bgcol, nrow(rooted$edge)))
    
    
    if (is.numeric(hlspecies)) {
      
      
      if (!is.null(outgroup)) {
        
        
        hlcols <- c(bgcol, hlcols)
        
        
        colMaster <- hlcols[dummyrooted$edge.length]
        
        
      }
      
      
      else {
        
        
        for (i in 1:length(hlspecies)) {
          
          
          colMaster[hlspecies[i]] <- hlcols[i]
          
          
        }
        
        
      }
      
      
    }
    
    
    else {
      
      
      wspmr <- rooted$tip.label[rooted$edge[, 2]]
      
      
      for (i in 1:length(hlspecies)) {
        
        
        colMaster[which(wspmr == hlspecies[i])] <- hlcols[i]
        
        
      }
      
      
    }
    
    
    termedge <- order(rooted$edge[, 2])[1:length(rooted$tip.label)]
    
    
    colMasterTip <- colMaster[termedge]
    
    
    rooted2 = rooted
    
    
    mm = min(rooted2$edge.length[rooted2$edge.length > 0])
    
    
    rooted2$edge.length[rooted2$edge.length == 0] = max(0.02, 
                                                        
                                                        
                                                        mm/20)
    
    
    plotobj = plot.phylo(rooted2, main = main, edge.color = colMaster, 
                         
                         
                         tip.color = colMasterTip, edge.width = 2, cex = fontSize)
    
    
    return(plotobj)
    
    
  }
  
  
}

