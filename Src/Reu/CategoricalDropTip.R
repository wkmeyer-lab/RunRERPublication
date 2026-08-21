categoricalDropTip = function (phy, tip, trim.internal = TRUE, root.edge = 0, 
          rooted = is.rooted(phy), collapse.singles = TRUE, interactive = FALSE, 
          ...) 
{
  Ntip <- length(phy$tip.label)
  if (is.character(tip)) 
    tip <- which(phy$tip.label %in% tip)
  out.of.range <- tip > Ntip
  if (any(out.of.range)) {
    warning("some tip numbers were larger than the number of tips: they were ignored")
    tip <- tip[!out.of.range]
  }
  
  if (!length(tip)) 
    return(phy)
  if (length(tip) == Ntip) {
    if (Nnode(phy) < 3 || trim.internal) {
      warning("drop all tips of the tree: returning NULL")
      return(NULL)
    }
  }
  
  
  wbl <- !is.null(phy$edge.length)
  
  if (!rooted) {
    phy$root.edge <- NULL
  }
  phy <- reorder(phy) 
  NEWROOT <- ROOT <- Ntip + 1 #Figure out the first internal node, by finding the node number one higher than the last terminal node (tip)
  Nnode <- phy$Nnode
  Nedge <- dim(phy$edge)[1]
  
  edge1 <- phy$edge[, 1] #store the parent nodes
  edge2 <- phy$edge[, 2] #store the daughter nodes
  keep <- !logical(Nedge) # make a logical of the length of the number of edges
  keep[match(tip, edge2)] <- FALSE #Mark all edges with a child node in the droplist to be removed 
  
  droppedInternals = integer()
  
  if (trim.internal) {
    ints <- edge2 > Ntip #make a list of the internal branches
    repeat {
      sel <- !(edge2 %in% edge1[keep]) & ints & keep #select branches which are 1)have a daughter node whose child node is NOT being kept; is an internal branch; and is currently being kept
      if (!sum(sel)) # If none were selected in this repeat
        break #end the loop
      droppedInternals = append(droppedInternals, which(sel))
      keep[sel] <- FALSE #if some were selected, stop keeping them
    }
    
    
    if (root.edge && wbl) {
      degree <- tabulate(edge1[keep])
      if (degree[ROOT] == 1) {
        j <- integer(0)
        repeat {
          i <- which(edge1 == NEWROOT & keep)
          j <- c(i, j)
          NEWROOT <- edge2[i]
          if (degree[NEWROOT] > 1) 
            break
        }
        keep[j] <- FALSE
        j <- j[1:root.edge]
        NewRootEdge <- sum(phy$edge.length[j])
        if (length(j) < root.edge && !is.null(phy$root.edge)) 
          NewRootEdge <- NewRootEdge + phy$root.edge
        phy$root.edge <- NewRootEdge
      }
    }
  }
  
  
  phy$edge <- phy$edge[keep, ] # remove all of the edges not being kept 
  
  ogEdgelengths = phy$edge.length
  
  if (wbl) 
    phy$edge.length <- phy$edge.length[keep] #remove all of the edge lengths of the removed edges from the edgelengths
  
  TERMS <- !(phy$edge[, 2] %in% phy$edge[, 1])  #Get a list of the edges with a daughter node NOT in the list of parent nodes
  oldNo.ofNewTips <- phy$edge[TERMS, 2] #This is a list of the nodes which are now terminal (tips)
  
  n <- length(oldNo.ofNewTips) #This is the number of tips in the new tree
  phy$edge[TERMS, 2] <- rank(phy$edge[TERMS, 2]) #Change the node number of the terminal branches daughter nodes to values starting at 1; because all of the branches already using those daugther numbers are included in this list, there are no repeats
  if (length(tip)) 
    phy$tip.label <- phy$tip.label[-tip]
  
  phy$Nnode <- dim(phy$edge)[1] - n + 1L #Change the new internal node number to be correct; the total number of edges, minus the number of terminal branches, plus 1 
  newNb <- integer(Ntip + Nnode) #Make an integer vector of length equal to the number of terminal branches and the number of internal nodes in the ORIGINAL tree; this is used to map old node values to new node values; where the index in this vector is the old number, and the value in this vector is the new number 
  newNb[NEWROOT] <- n + 1L #Set the value at the position of the root in the original tree to be the new first internal node 
  sndcol <- phy$edge[, 2] > n # get a list of edges with a daughter node that is not a terminal branch (internal branches)
  newNb[sort(phy$edge[sndcol, 2])] <- (n + 2):(n + phy$Nnode) #for all of the internal branches, sort them by their daughter node (Because they would otherwise be ordered by parent node). Then, in the integer vector, set a new value made up of internal node numbers starting at one past the new root node's number. 
        #Effectively, this makes a vector that makes new internal node numbers in acceding order; by using only the indexes from the branches, it also skips any internal nodes that no longer are used. 
  phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2]] #Set the daughter node values of the internal branches to the now decreased node values
  phy$edge[, 1] <- newNb[phy$edge[, 1]] #Set the parent nodes values to the new node values using the conversion table
  storage.mode(phy$edge) <- "integer"
  if (!is.null(phy$node.label))  #IF there are node labels 
    phy$node.label <- phy$node.label[which(newNb > 0) - Ntip] #Reduce the number of node labels by dropping any labels with a conversion table VALUE (aka new node number) of 0 (not in the new tree)
  
  
  #Alright, at this point, all of the nodes have been re-assigned. The problem is that the branch lengths are still wrong -- 
  phyUnmerged = phy
  if (collapse.singles)
  {
    
    tree = phy
    
    n <- length(tree$tip.label)
    if (n == 0) {  #if there are no tip labels
      return(tree) #exit early 
    }
    
    tree <- reorder(tree) #reorder the tree (no effect) 
    e1 <- tree$edge[, 1] #save the parent nodes
    e2 <- tree$edge[, 2] #save the daughter nodes
    tab <- tabulate(e1) #make a table of how many times each number between 1 and the largest parent node is a parent node. 
    
    
    if (all(tab[-c(1:n)] > 1)) #if each internal node (numbers larger than the number of tips) is only a parent once
      return(tree) #exit early
    if (is.null(tree$edge.length)) {
      wbl <- FALSE
    }else {
      wbl <- TRUE
      el <- tree$edge.length #save the edgelengths 
    }
    
    
    ROOT <- n + 1L   #Get the root node, the first internal node 
    while (tab[ROOT] == 1) { #code for if the root node is parent only once 
      i <- which(e1 == ROOT) #Remove nodes and branches until the root is parent more than once 
      ROOT <- e2[i]
      if (wbl) {
        el <- el[-i]
      }
      e1 <- e1[-i]
      e2 <- e2[-i]
    }
    
    edgeFrame = data.frame(1:length(e1), e1,e2,el)
    names(edgeFrame) = c("index", "ogE1", "ogE2", "ogEl")
    
    singles <- which(tabulate(e1) == 1) #get the nodes which are parents only once 
    if (length(singles) > 0) {
      ii <- sort(match(singles, e1), decreasing = TRUE) #get the indexes of branches which are singles 
      jj <- match(e1[ii], e2) #Get the indexes of branches with a daughter node that is the parent to a single 
      for (i in 1:length(singles)) {
        e2[jj[i]] <- e2[ii[i]] #set the daughter of 'the branch whose daughter is a single' to the daughter of that single 
        if (wbl) 
          el[jj[i]] <- el[ii[i]] #set the branch length of that branch to the daughter's length (CHANGED FROM SUMMING)
      }
      e1 <- e1[-ii] #remove the single from the parent list
      e2 <- e2[-ii] #remove the single from the daughter list
      if (wbl) 
        el <- el[-ii] #remove the single from the branch length list 
      edgeFrame = edgeFrame[-ii,]
    }
    
    
    Nnode <- length(e1) - n + 1L
    oldnodes <- unique(e1)
    if (!is.null(tree$node.label)) 
      tree$node.label <- tree$node.label[oldnodes - n]
    newNb <- integer(max(oldnodes))
    newNb[ROOT] <- n + 1L
    sndcol <- e2 > n
    e2[sndcol] <- newNb[e2[sndcol]] <- n + 2:Nnode
    e1 <- newNb[e1]
    tree$edge <- cbind(e1, e2, deparse.level = 0)
    tree$Nnode <- Nnode
    edgeFrame$e1 = e1
    edgeFrame$e2 = e2
    edgeFrame$el = el
    if (wbl) {
      tree$edge.length <- el
    }
    tree
    
    phy = tree
  }
    

  phy
}
