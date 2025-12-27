# This file is scartch sapce for permualtions updates being made on the main RER repo. 
a = b #prevent full runs
library(RERconverge)
library(tools)

source("../RERConverge/R/RERfuncs.R")
source("../RERConverge/R/PermulationFuncs.R")

# -----------------------------------------------------------------
# --- Needed internal functions ---
# -----------------------------------------------------------------
{
  #' @keywords internal
  getForegroundsFromBinaryTree=function(tree){
    nameEdgesPerms.tree = nameEdgesPerms(tree)
    edge.length = as.logical(tree$edge.length)
    foregrounds = nameEdgesPerms.tree[edge.length]
    ind.tip = which(foregrounds != "")
    foregrounds = foregrounds[ind.tip]
    return(foregrounds)
  }
  
  #' @keywords internal
  nameEdgesPerms=function(tree){
    if (is.null(tree$tip.label)) {
      nn = NULL
    } else {
      nn=character(nrow(tree$edge))
      iim=match(1:length(tree$tip.label), tree$edge[,2])
      nn[iim]=tree$tip.label
    }
    nn
  }
}

# ------------------------------------------------------------------
# --- Updated function --- 
# ------------------------------------------------------------------

getPermsBinary=function(numperms, fg_vec, sisters_list, root_sp, RERmat, trees, mastertree, permmode="cc", method="k", min.pos=2, trees_list=NULL, calculateenrich=F, annotlist=NULL){
  pathvec = foreground2Paths(fg_vec, trees, clade="all",plotTree=F)
  col_labels = colnames(trees$paths)
  names(pathvec) = col_labels
  
  message("As of RERConverge [X.xx], permulation functions have been updated. Old versions have been moved to ccLegacy and ssmLegacy.")
  if(permmode=="cc"){
    print("Running CC permulation. sisters_list is required only for enrichments, otherwise sisters_list = NA is sufficient.")
    
    print("Generating permulated trees")
    
    # --- new code since legacy method; switching to categorical function to infer phenotype tree --
    #covert fg_vec to a categorical phenotypeVector
    phenotypeVector = rep(0, length(trees$masterTree$tip.label))
    names(phenotypeVector) = trees$masterTree$tip.label
    phenotypeVector[names(phenotypeVector) %in% fg_vec] = 1
    
    
    permulationData = categoricalPermulations(trees, phenotypeVector, rm = "ER", rp = "auto", ntrees = numperms)
    message("cat perm done")
    
    permulatedTrees = lapply(permulationData$trees, function(x) {
      tr = trees$masterTree
      tr$edge.length = c(x$tips, x$nodes)[tr$edge[,2]]
      tr$edge.length = tr$edge.length-1
      names(tr$edge.length) = NULL
      #tree2Paths(tr, treesObj, categorical = TRUE, useSpecies = names(phenvals))
      tr
    })
    permulated.binphens = list(permulatedTrees)
    message(permulated.binphens)
    message("cat conversion done")
    #----
    #permulated.binphens = generatePermulatedBinPhen(trees$masterTree, numperms, trees, root_sp, fg_vec, sisters_list, pathvec, permmode="cc")
    permulated.fg = mapply(getForegroundsFromBinaryTree, permulated.binphens[[1]])
    permulated.fg.list = as.list(data.frame(permulated.fg))
    phenvec.table = mapply(foreground2Paths,permulated.fg.list,MoreArgs=list(treesObj=trees,clade="all"))
    phenvec.list = lapply(seq_len(ncol(phenvec.table)), function(i) phenvec.table[,i])
    
    print("Calculating correlations")
    corMatList = lapply(phenvec.list, correlateWithBinaryPhenotype, RERmat=RERmat)
    
    #make enrich list/matrices to fill
    permPvals=data.frame(matrix(ncol=numperms, nrow=nrow(RERmat)))
    rownames(permPvals)=rownames(RERmat)
    permRhovals=data.frame(matrix(ncol=numperms, nrow=nrow(RERmat)))
    rownames(permRhovals)=rownames(RERmat)
    permStatvals=data.frame(matrix(ncol=numperms, nrow=nrow(RERmat)))
    rownames(permStatvals)=rownames(RERmat)
    
    for (i in 1:length(corMatList)){
      permPvals[,i] = corMatList[[i]]$P
      permRhovals[,i] = corMatList[[i]]$Rho
      permStatvals[,i] = sign(corMatList[[i]]$Rho)*-log10(corMatList[[i]]$P)
    }
    
  }
  else if (permmode=="ccLegacy"){
    print("Running CC Legacy permulation")
    
    print("Generating permulated trees")
    permulated.binphens = generatePermulatedBinPhen(trees$masterTree, numperms, trees, root_sp, fg_vec, sisters_list, pathvec, permmode="cc")
    permulated.fg = mapply(getForegroundsFromBinaryTree, permulated.binphens[[1]])
    permulated.fg.list = as.list(data.frame(permulated.fg))
    phenvec.table = mapply(foreground2Paths,permulated.fg.list,MoreArgs=list(treesObj=trees,clade="all"))
    phenvec.list = lapply(seq_len(ncol(phenvec.table)), function(i) phenvec.table[,i])
    
    print("Calculating correlations")
    corMatList = lapply(phenvec.list, correlateWithBinaryPhenotype, RERmat=RERmat)
    
    #make enrich list/matrices to fill
    permPvals=data.frame(matrix(ncol=numperms, nrow=nrow(RERmat)))
    rownames(permPvals)=rownames(RERmat)
    permRhovals=data.frame(matrix(ncol=numperms, nrow=nrow(RERmat)))
    rownames(permRhovals)=rownames(RERmat)
    permStatvals=data.frame(matrix(ncol=numperms, nrow=nrow(RERmat)))
    rownames(permStatvals)=rownames(RERmat)
    
    for (i in 1:length(corMatList)){
      permPvals[,i] = corMatList[[i]]$P
      permRhovals[,i] = corMatList[[i]]$Rho
      permStatvals[,i] = sign(corMatList[[i]]$Rho)*-log10(corMatList[[i]]$P)
    }
    
  } else if (permmode=="ssmLegacy"){
    print("Running SSM Legacy permulation")
    
    if (is.null(trees_list)){
      trees_list = trees$trees
    }
    
    RERmat = RERmat[match(names(trees_list), rownames(RERmat)),]
    
    print("Generating permulated trees")
    permulated.binphens = generatePermulatedBinPhenSSMBatched(trees_list,numperms,trees,root_sp,fg_vec,sisters_list,pathvec)
    
    # Get species membership of the trees
    df.list = lapply(trees_list,getSpeciesMembershipStats,masterTree=mastertree,foregrounds=fg_vec)
    df.converted = data.frame(matrix(unlist(df.list), nrow=length(df.list), byrow=T),stringsAsFactors=FALSE)
    attr = attributes(df.list[[1]])
    col_names = attr$names
    attr2 = attributes(df.list)
    row_names = attr2$names
    
    colnames(df.converted) = col_names
    rownames(df.converted) = row_names
    
    df.converted$num.fg = as.integer(df.converted$num.fg)
    df.converted$num.spec = as.integer(df.converted$num.spec)
    
    spec.members = df.converted$spec.members
    
    # Group gene trees based on the similarity of their species membership
    grouped.trees = groupTrees(spec.members)
    ind.unique.trees = grouped.trees$ind.unique.trees
    ind.unique.trees = unlist(ind.unique.trees)
    ind.tree.groups = grouped.trees$ind.tree.groups
    
    # For each unique tree, produce a permuted tree. We already have this function, but we need a list of trees to feed in.
    unique.trees = trees_list[ind.unique.trees]
    
    # precompute clade mapping for each unique tree
    unique.map.list = mapply(matchAllNodesClades,unique.trees,MoreArgs=list(treesObj=trees))
    
    # calculate paths for each permulation
    unique.permulated.binphens = permulated.binphens[ind.unique.trees]
    unique.permulated.paths = calculatePermulatedPaths_apply(unique.permulated.binphens,unique.map.list,trees)
    
    permulated.paths = vector("list", length = length(trees_list))
    for (j in 1:length(permulated.paths)){
      permulated.paths[[j]] = vector("list",length=numperms)
    }
    for (i in 1:length(unique.permulated.paths)){
      ind.unique.tree = ind.unique.trees[i]
      ind.tree.group = ind.tree.groups[[i]]
      unique.path = unique.permulated.paths[[i]]
      for (k in 1:length(ind.tree.group)){
        permulated.paths[[ind.tree.group[k]]] = unique.path
      }
    }
    attributes(permulated.paths)$names = row_names
    
    print("Calculating correlations")
    RERmat.list = lapply(seq_len(nrow(RERmat[])), function(i) RERmat[i,])
    corMatList = mapply(calculateCorPermuted,permulated.paths,RERmat.list)
    permPvals = extractCorResults(corMatList,numperms,mode="P")
    rownames(permPvals) = names(trees_list)
    permRhovals = extractCorResults(corMatList,numperms,mode="Rho")
    rownames(permRhovals) = names(trees_list)
    permStatvals = sign(permRhovals)*-log10(permPvals)
    rownames(permStatvals) = names(trees_list)
    
  } else {
    stop("Invalid binary permulation mode.")
  }
  
  if (calculateenrich){
    realFgtree = foreground2TreeClades(fg_vec, sisters_list, trees, plotTree=F)
    realpaths = tree2PathsClades(realFgtree, trees)
    realresults = getAllCor(RERmat, realpaths, method=method, min.pos=min.pos)
    realstat =sign(realresults$Rho)*-log10(realresults$P)
    names(realstat) = rownames(RERmat)
    realenrich = fastwilcoxGMTall(na.omit(realstat), annotlist, outputGeneVals=F)
    
    #sort real enrichments
    groups=length(realenrich)
    c=1
    while(c<=groups){
      current=realenrich[[c]]
      realenrich[[c]]=current[order(rownames(current)),]
      c=c+1
    }
    #make matrices to fill
    permenrichP=vector("list", length(realenrich))
    permenrichStat=vector("list", length(realenrich))
    c=1
    while(c<=length(realenrich)){
      newdf=data.frame(matrix(ncol=numperms, nrow=nrow(realenrich[[c]])))
      rownames(newdf)=rownames(realenrich[[c]])
      permenrichP[[c]]=newdf
      permenrichStat[[c]]=newdf
      c=c+1
    }
    
    counter=1;
    while (counter <= numperms){
      stat = permStatvals[,counter]
      names(stat) = rownames(RERmat)
      enrich=fastwilcoxGMTall(na.omit(stat), annotlist, outputGeneVals=F)
      #sort and store enrichment results
      groups=length(enrich)
      c=1
      while(c<=groups){
        current=enrich[[c]]
        enrich[[c]]=current[order(rownames(current)),]
        enrich[[c]]=enrich[[c]][match(rownames(permenrichP[[c]]), rownames(enrich[[c]])),]
        permenrichP[[c]][,counter]=enrich[[c]]$pval
        permenrichStat[[c]][,counter]=enrich[[c]]$stat
        c=c+1
      }
      counter = counter+1
    }
  }
  
  if(calculateenrich){
    data=vector("list", 5)
    data[[1]]=permPvals
    data[[2]]=permRhovals
    data[[3]]=permStatvals
    data[[4]]=permenrichP
    data[[5]]=permenrichStat
    names(data)=c("corP", "corRho", "corStat", "enrichP", "enrichStat")
  } else {
    data=vector("list", 3)
    data[[1]]=permPvals
    data[[2]]=permRhovals
    data[[3]]=permStatvals
    names(data)=c("corP", "corRho", "corStat")
  }
  data
}


# ------------------------------------------------------------------
# --- Testing run  --- 
# ------------------------------------------------------------------
rerpath = find.package('RERconverge')
toytreefile = "subsetMammalGeneTrees.txt" 
toyTrees=readTrees(paste(rerpath,"/extdata/",toytreefile,sep=""), max.read = 200)
marineFg = c("Killer_whale", "Dolphin", "Walrus", "Seal", "Manatee")
sisters_marine = list("clade1"=c("Killer_whale", "Dolphin"))
marineFgTree = foreground2TreeClades(marineFg,sisters_marine,toyTrees,plotTree=F)

# Calculating paths from the foreground tree
pathvec = tree2PathsClades(marineFgTree, toyTrees)
# Calculate RERs
mamRERw = getAllResiduals(toyTrees, transform="sqrt", weighted=T, scale=T)
# Calculate correlation 
res = correlateWithBinaryPhenotype(mamRERw, pathvec, min.sp=10, min.pos=2,
                                   weighted="auto")
#define the root species
root_sp = "Human"


permCC = getPermsBinary(10, marineFg, sisters_marine, root_sp, mamRERw, toyTrees,
                        masterTree, permmode="cc")

permCC = getPermsBinary(10, marineFg, sisters_marine, root_sp, mamRERw, toyTrees,
                        masterTree, permmode="ccLegacy")

# ------------------------------------------------------------------
# --- Looking into things for the permualtions update  ----- 
# ------------------------------------------------------------------
outputFolderName = "Output/CategoricalInsVertivoreTree/"
filePrefix = "CategoricalInsVertivoreTree"

phenotypeVectorFilename = paste(outputFolderName, filePrefix, "CategoricalPhenotypeVector.rds",sep="") #make a filename based on the prefix
phenotypeVector = readRDS(phenotypeVectorFilename)                              #Load the phenotype vector 


?generatePermulatedBinPhen()
categoricalPermulations
correlateWithBinaryPhenotype

# -- 
# toy example 

library(RERconverge)
rerpath = find.package('RERconverge')

#read trees
toytreefile = "subsetMammalGeneTrees.txt" 
toyTrees=readTrees(paste(rerpath,"/extdata/",toytreefile,sep=""), max.read = 200)

marineFg = c("Killer_whale", "Dolphin", "Walrus", "Seal", "Manatee")
sisters_marine = list("clade1"=c("Killer_whale", "Dolphin"))

marineFgTree = foreground2TreeClades(marineFg,sisters_marine,toyTrees,plotTree=F)
marineplot1 = plotTreeHighlightBranches(marineFgTree,
                                        hlspecies=which(marineFgTree$edge.length==1),
                                        hlcols="blue", main="Marine mammals trait tree")


# Calculating paths from the foreground tree
pathvec = tree2PathsClades(marineFgTree, toyTrees)

# Calculate RERs
mamRERw = getAllResiduals(toyTrees, transform="sqrt", weighted=T, scale=T)

# Calculate correlation 
res = correlateWithBinaryPhenotype(mamRERw, pathvec, min.sp=10, min.pos=2,
                                   weighted="auto")


#define the root species
root_sp = "Human"

masterTree = toyTrees$masterTree


numperms = 10
fg_vec = marineFg
sisters_list = sisters_marine
RERmat = mamRERw
trees = toyTrees

phenotypeVector = rep(0, length(trees$masterTree$tip.label))
names(phenotypeVector) = trees$masterTree$tip.label
phenotypeVector[names(phenotypeVector) %in% fg_vec] = 1


permulated.binphens = generatePermulatedBinPhen(trees$masterTree, numperms, trees, root_sp, fg_vec, sisters_list, pathvec, permmode="cc")

test = permulated.binphens[[1]]

test[[1]][[1]]

all.equal(test[[1]][[1]], test[[10]][[1]])
all.equal(test[[1]][[3]], test[[10]][[3]])
all.equal(test[[1]][[4]], test[[10]][[4]])
# This has confirmed that the trees are all identical except for the branch lengths. 


permulationData = categoricalPermulations(trees, phenotypeVector, rm = "ER", rp = "auto", ntrees = numperms)
#Alright that produces something, but not whole trees -- which is what the funciton is otherwise expecting. See if I can convert that into full trees. 
length(test[[1]][[2]]) # that's one less than the number of tips plus the number of nodes, so the tips+nodes don't directly translate into edges, but they are close. 

demoTree = permulationData$trees[[1]]

permulatedTrees = lapply(permulationData$trees, function(x) {
  tr = trees$masterTree
  tr$edge.length = c(x$tips, x$nodes)[tr$edge[,2]]
  tr$edge.length = tr$edge.length-1
  names(tr$edge.length) = NULL
  #tree2Paths(tr, treesObj, categorical = TRUE, useSpecies = names(phenvals))
  tr
})

permulatedTrees[[1]]$edge.length

tr = trees$masterTree
tr$edge.length = c(demoTree$tips, demoTree$nodes)[tr$edge[,2]]






function (realCors, nullPhens, phenvals, treesObj, RERmat, method = "kw", 
          min.sp = 10, min.pos = 2, winsorizeRER = NULL, winsorizetrait = NULL, 
          weighted = F, extantOnly = FALSE, report=F) 
{
  tree = treesObj$masterTree
  keep = intersect(names(phenvals), tree$tip.label)
  tree = pruneTree(tree, keep)
  if (is.rooted(tree)) {
    tree = unroot(tree)
  }
  if(report){pathStartTime = Sys.time()}
  message("Generating null paths")
  nullPaths = lapply(nullPhens, function(x) {
    if(report){message("One path complete")}
    tr = tree
    tr$edge.length = c(x$tips, x$nodes)[tr$edge[,2]]
    tree2Paths(tr, treesObj, categorical = TRUE, useSpecies = names(phenvals))
  })
  if(report){pathsEndTime = Sys.time(); pathsDuration = pathsEndTime - pathStartTime; message(paste("Completed paths;","Duration", pathsDuration, attr(pathsDuration, "units")))}
  
  message("Calculating correlation statistics")
  corsMatPvals = matrix(nrow = nrow(RERmat), ncol = length(nullPhens), dimnames = list(rownames(RERmat), NULL))
  corsMatEffSize = matrix(nrow = nrow(RERmat), ncol = length(nullPhens), dimnames = list(rownames(RERmat), NULL))
  if(report){message("Matrixes")}
  Ppvals = lapply(1:length(realCors[[2]]), matrix, data = NA, nrow = nrow(RERmat), ncol = length(nullPhens), dimnames = list(rownames(RERmat), NULL))
  names(Ppvals) = names(realCors[[2]])
  Peffsize = lapply(1:length(realCors[[2]]), matrix, data = NA, nrow = nrow(RERmat), ncol = length(nullPhens), dimnames = list(rownames(RERmat), NULL))
  names(Peffsize) = names(realCors[[2]])
  if(report){message("pVals")}
  for (i in 1:length(nullPaths)) {
    if(report){corStartTime = Sys.time()}
    cors = getAllCor(RERmat, nullPaths[[i]], method = method, 
                     min.sp = min.sp, min.pos = min.pos, winsorizeRER = winsorizeRER, 
                     winsorizetrait = winsorizetrait, weighted = weighted)
    if(report){corEndTime = Sys.time(); corDuration = corEndTime - corStartTime; message(paste("Completed Correlation", i, "Duration", corDuration, attr(corDuration, "units")))}
    corsMatPvals[, i] = cors[[1]]$P
    corsMatEffSize[, i] = cors[[1]]$Rho
    for (j in 1:length(cors[[2]])) {
      Ppvals[[names(cors[[2]])[j]]][, i] = cors[[2]][[j]]$P
      Peffsize[[names(cors[[2]])[j]]][, i] = cors[[2]][[j]]$Rho
    }
    #if(report){message(paste("compelted", i))}
    gc()
  }
  output = list(corsMatEffSize, Peffsize, corsMatPvals, Ppvals)
  names(output) = c("corsMatEffSize", "Peffsize", "corsMatPvals", "Ppvals")
  return(output)
}



#perform binary CC permulation
permCC = getPermsBinary(10, marineFg, sisters_marine, root_sp, mamRERw, toyTrees,
                        masterTree, permmode="cc")

#calculate permulation p-values
permpvalCC = permpvalcor(res,permCC)