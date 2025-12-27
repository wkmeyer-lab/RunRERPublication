library(RERconverge) #load RERconverge package
library(RERconverge)
library("tools")
library(fastmatch)
library(bench)
source("Src/Reu/cmdArgImport.R")
source("Src/Reu/convertLogiToNumeric.R")

function (numperms, fg_vec, sisters_list, root_sp, RERmat, trees, 
          mastertree, permmode = "cc", method = "k", min.pos = 2, trees_list = NULL, 
          calculateenrich = F, annotlist = NULL) 
fg_vec = foregroundString
sisters_list = sistersList  
trees = mainTrees  
mastertree = masterTree
numperms = 3
root_sp = rootSpecies
RERmat = RERObject
{
  pathvec = foreground2Paths(fg_vec, trees, clade = "all", 
                             plotTree = F)
  col_labels = colnames(trees$paths)
  names(pathvec) = col_labels
  if (permmode == "cc") {
    print("Running CC permulation")
    print("Generating permulated trees")
    permulated.binphens = generatePermulatedBinPhen(trees$masterTree, 
                                                    numperms, trees, root_sp, fg_vec, sisters_list, pathvec, 
                                                    permmode = "cc")
    permulated.fg = mapply(getForegroundsFromBinaryTree, 
                           permulated.binphens[[1]])
    permulated.fg = mapply(getForegroundsFromBinaryTree, 
                           permulated.binphens[[1]])
    permulated.fg.list = as.list(data.frame(permulated.fg))
    phenvec.table = mapply(foreground2Paths, permulated.fg.list, 
                           MoreArgs = list(treesObj = trees, clade = "all"))
    phenvec.list = lapply(seq_len(ncol(phenvec.table)), function(i) phenvec.table[, 
                                                                                  i])
    print("Calculating correlations")
    corMatList = lapply(phenvec.list, correlateWithBinaryPhenotype, 
                        RERmat = RERmat)
    permPvals = data.frame(matrix(ncol = numperms, nrow = nrow(RERmat)))
    rownames(permPvals) = rownames(RERmat)
    permRhovals = data.frame(matrix(ncol = numperms, nrow = nrow(RERmat)))
    rownames(permRhovals) = rownames(RERmat)
    permStatvals = data.frame(matrix(ncol = numperms, nrow = nrow(RERmat)))
    rownames(permStatvals) = rownames(RERmat)
    for (i in 1:length(corMatList)) {
      permPvals[, i] = corMatList[[i]]$P
      permRhovals[, i] = corMatList[[i]]$Rho
      permStatvals[, i] = sign(corMatList[[i]]$Rho) * -log10(corMatList[[i]]$P)
    }
  }
  
  
  generatePermulatedBinPhen
  function (tree, numperms, trees, root_sp, fg_vec, sisters_list, 
            pathvec, permmode = "cc") 
  {
    if (permmode == "cc") {
      tree_rep = lapply(1:numperms, rep_tree, tree = trees)
      permulated.binphens = lapply(tree_rep, simBinPhenoCC, 
                                   mastertree = trees$masterTree, root_sp = root_sp, 
                                   fg_vec = fg_vec, sisters_list = sisters_list, pathvec = pathvec, 
                                   plotTreeBool = F)
    }
    
    
    simBinPhenoCC
    function (trees, mastertree, root_sp, fg_vec, sisters_list = NULL, 
              pathvec, plotTreeBool = F) 
    {
      tip.labels = mastertree$tip.label
      res = getForegroundInfoClades(fg_vec, sisters_list, trees, 
                                    plotTree = F, useSpecies = tip.labels)
      fg_tree = res$tree
      fg.table = res$fg.sisters.table
      fgnum = length(which(fg_tree$edge.length == 1))
      internal = nrow(fg.table)
      tips = fgnum - internal
      num.tip.sisters.real = length(which(as.vector(fg.table) <= 
                                            length(tip.labels)))
      top = NA
      num.tip.sisters.fake = 10000
      while (num.tip.sisters.fake != num.tip.sisters.real) {
        blsum = 0
        while (blsum != fgnum) {
          t = root.phylo(trees$masterTree, root_sp, resolve.root = T)
          rm = ratematrix(t, pathvec)
          sims = sim.char(t, rm, nsim = 1)
          nam = rownames(sims)
          s = as.data.frame(sims)
          simulatedvec = s[, 1]
          names(simulatedvec) = nam
          top = names(sort(simulatedvec, decreasing = TRUE))[1:tips]
          t = foreground2Tree(top, trees, clade = "all", plotTree = F)
          blsum = sum(t$edge.length)
          t.table = findPairs(t)
          num.tip.sisters.fake = length(which(as.vector(t.table) <= 
                                                length(tip.labels)))
        }
      }
      if (plotTreeBool) {
        plot(t)
      }
      return(t)
    }
    
    
    
    
    else if (permmode == "ssm") {
      tree_rep = lapply(1:numperms, rep_tree, tree = tree)
      permulated.binphens = lapply(tree_rep, simBinPhenoSSM, 
                                   trees = trees, root_sp = root_sp, fg_vec = fg_vec, 
                                   sisters_list = sisters_list, pathvec = pathvec)
    }
    else {
      stop("Invalid binary permulation mode.")
    }
    output.list <- list()
    output.list[[1]] <- permulated.binphens
    return(output.list)
  }

  else {
    data = vector("list", 3)
    data[[1]] = permPvals
    data[[2]] = permRhovals
    data[[3]] = permStatvals
    names(data) = c("corP", "corRho", "corStat")
  }
  data
}


getPermulatedFg = function (numperms, fg_vec, sisters_list, root_sp, RERmat, trees, 
          mastertree, permmode = "cc", method = "k", min.pos = 2, trees_list = NULL, 
          calculateenrich = F, annotlist = NULL) 
{
  environment() = environment(getPermsBinary) 
  pathvec = foreground2Paths(fg_vec, trees, clade = "all", 
                             plotTree = F)
  col_labels = colnames(trees$paths)
  names(pathvec) = col_labels
  if (permmode == "cc") {
    print("Running CC permulation")
    print("Generating permulated trees")
    permulated.binphens = generatePermulatedBinPhen(trees$masterTree, 
                                                    numperms, trees, root_sp, fg_vec, sisters_list, pathvec, 
                                                    permmode = "cc")
    permulated.fg = mapply(getForegroundsFromBinaryTree, 
                           permulated.binphens[[1]])
    permulated.fg.list = as.list(data.frame(permulated.fg))
    phenvec.table = mapply(foreground2Paths, permulated.fg.list, 
                           MoreArgs = list(treesObj = trees, clade = "all"))
    phenvec.list = lapply(seq_len(ncol(phenvec.table)), function(i) phenvec.table[, 
                                                                                  i])
    print("Calculating correlations")
    return(phenvec.list)
    corMatList = lapply(phenvec.list, correlateWithBinaryPhenotype, 
                        RERmat = RERmat)
    permPvals = data.frame(matrix(ncol = numperms, nrow = nrow(RERmat)))
    rownames(permPvals) = rownames(RERmat)
    permRhovals = data.frame(matrix(ncol = numperms, nrow = nrow(RERmat)))
    rownames(permRhovals) = rownames(RERmat)
    permStatvals = data.frame(matrix(ncol = numperms, nrow = nrow(RERmat)))
    rownames(permStatvals) = rownames(RERmat)
    for (i in 1:length(corMatList)) {
      permPvals[, i] = corMatList[[i]]$P
      permRhovals[, i] = corMatList[[i]]$Rho
      permStatvals[, i] = sign(corMatList[[i]]$Rho) * -log10(corMatList[[i]]$P)
    }
  }
  else if (permmode == "ssm") {
    print("Running SSM permulation")
    if (is.null(trees_list)) {
      trees_list = trees$trees
    }
    RERmat = RERmat[match(names(trees_list), rownames(RERmat)), 
    ]
    print("Generating permulated trees")
    permulated.binphens = generatePermulatedBinPhenSSMBatched(trees_list, 
                                                              numperms, trees, root_sp, fg_vec, sisters_list, pathvec)
    df.list = lapply(trees_list, getSpeciesMembershipStats, 
                     masterTree = masterTree, foregrounds = fg_vec)
    df.converted = data.frame(matrix(unlist(df.list), nrow = length(df.list), 
                                     byrow = T), stringsAsFactors = FALSE)
    attr = attributes(df.list[[1]])
    col_names = attr$names
    attr2 = attributes(df.list)
    row_names = attr2$names
    colnames(df.converted) = col_names
    rownames(df.converted) = row_names
    df.converted$num.fg = as.integer(df.converted$num.fg)
    df.converted$num.spec = as.integer(df.converted$num.spec)
    spec.members = df.converted$spec.members
    grouped.trees = groupTrees(spec.members)
    ind.unique.trees = grouped.trees$ind.unique.trees
    ind.unique.trees = unlist(ind.unique.trees)
    ind.tree.groups = grouped.trees$ind.tree.groups
    unique.trees = trees_list[ind.unique.trees]
    unique.map.list = mapply(matchAllNodesClades, unique.trees, 
                             MoreArgs = list(treesObj = trees))
    unique.permulated.binphens = permulated.binphens[ind.unique.trees]
    unique.permulated.paths = calculatePermulatedPaths_apply(unique.permulated.binphens, 
                                                             unique.map.list, trees)
    permulated.paths = vector("list", length = length(trees_list))
    for (j in 1:length(permulated.paths)) {
      permulated.paths[[j]] = vector("list", length = numperms)
    }
    for (i in 1:length(unique.permulated.paths)) {
      ind.unique.tree = ind.unique.trees[i]
      ind.tree.group = ind.tree.groups[[i]]
      unique.path = unique.permulated.paths[[i]]
      for (k in 1:length(ind.tree.group)) {
        permulated.paths[[ind.tree.group[k]]] = unique.path
      }
    }
    attributes(permulated.paths)$names = row_names
    print("Calculating correlations")
    RERmat.list = lapply(seq_len(nrow(RERmat[])), function(i) RERmat[i, 
    ])
    corMatList = mapply(calculateCorPermuted, permulated.paths, 
                        RERmat.list)
    permPvals = extractCorResults(corMatList, numperms, mode = "P")
    rownames(permPvals) = names(trees_list)
    permRhovals = extractCorResults(corMatList, numperms, 
                                    mode = "Rho")
    rownames(permRhovals) = names(trees_list)
    permStatvals = sign(permRhovals) * -log10(permPvals)
    rownames(permStatvals) = names(trees_list)
  }
  else {
    stop("Invalid binary permulation mode.")
  }
  if (calculateenrich) {
    realFgtree = foreground2TreeClades(fg_vec, sisters_list, 
                                       trees, plotTree = F)
    realpaths = tree2PathsClades(realFgtree, trees)
    realresults = getAllCor(RERmat, realpaths, method = method, 
                            min.pos = min.pos)
    realstat = sign(realresults$Rho) * -log10(realresults$P)
    names(realstat) = rownames(RERmat)
    realenrich = fastwilcoxGMTall(na.omit(realstat), annotlist, 
                                  outputGeneVals = F)
    groups = length(realenrich)
    c = 1
    while (c <= groups) {
      current = realenrich[[c]]
      realenrich[[c]] = current[order(rownames(current)), 
      ]
      c = c + 1
    }
    permenrichP = vector("list", length(realenrich))
    permenrichStat = vector("list", length(realenrich))
    c = 1
    while (c <= length(realenrich)) {
      newdf = data.frame(matrix(ncol = numperms, nrow = nrow(realenrich[[c]])))
      rownames(newdf) = rownames(realenrich[[c]])
      permenrichP[[c]] = newdf
      permenrichStat[[c]] = newdf
      c = c + 1
    }
    counter = 1
    while (counter <= numperms) {
      stat = permStatvals[, counter]
      names(stat) = rownames(RERmat)
      enrich = fastwilcoxGMTall(na.omit(stat), annotlist, 
                                outputGeneVals = F)
      groups = length(enrich)
      c = 1
      while (c <= groups) {
        current = enrich[[c]]
        enrich[[c]] = current[order(rownames(current)), 
        ]
        enrich[[c]] = enrich[[c]][match(rownames(permenrichP[[c]]), 
                                        rownames(enrich[[c]])), ]
        permenrichP[[c]][, counter] = enrich[[c]]$pval
        permenrichStat[[c]][, counter] = enrich[[c]]$stat
        c = c + 1
      }
      counter = counter + 1
    }
  }
  if (calculateenrich) {
    data = vector("list", 5)
    data[[1]] = permPvals
    data[[2]] = permRhovals
    data[[3]] = permStatvals
    data[[4]] = permenrichP
    data[[5]] = permenrichStat
    names(data) = c("corP", "corRho", "corStat", "enrichP", 
                    "enrichStat")
  }
  else {
    data = vector("list", 3)
    data[[1]] = permPvals
    data[[2]] = permRhovals
    data[[3]] = permStatvals
    names(data) = c("corP", "corRho", "corStat")
  }
  data
}





getForegroundsFromBinaryTree=function(tree){
  nameEdgesPerms.tree = nameEdgesPerms(tree)
  edge.length = as.logical(tree$edge.length)
  foregrounds = nameEdgesPerms.tree[edge.length]
  ind.tip = which(foregrounds != "")
  foregrounds = foregrounds[ind.tip]
  return(foregrounds)
}

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
rm(corList2)
corList2 = list()
corList2 = append(corList2, list(correl))
corList2 = append(corList2, correl)
permulationsfg = getPermulatedFg(permulationNumber, foregroundString, sistersList, rootSpecies, RERObject, mainTrees, mastertree =  masterTree, permmode = "cc")

testListFunction = function(n, mList= list()){
  mList = append(mList, n)
  mList = append(mList, "a")
}
testList = testListFunction(n=1)
leafBitMaps = NULL
phenvec = NULL; rm=NULL; leafBitMaps=NULL
tree=rootedMasterTree; phenvec=phenotypeVector; internal=internalNumber
function(tree, tips=0, internal, phenvec = NULL, rm=NULL, leafBitMaps=NULL){
  insum=0
  top = character(0)
  if (!is.null(phenvec)) {
    print("phenvec present")
    tips = sum(phenvec)
    if(is.null(rm)) {
      print("Making rate matrix")
      rm = ratematrix(tree, phenvec)
    }
  } else if (is.null(rm)){
    print("No rate matrix supplied & cannot compute one")
    return()
  }
  
  if(is.null(leafBitMaps)){
    leafBitMaps = makeLeafMap(tree)
  }
  
  repeat { #implements a do-while loop
    t=tree 
    #root.phylo(tree, root, resolve.root = T) #For now, let's use the ancestral root
    sims=sim.char(t, rm, nsim = 1)
    nam=rownames(sims)
    s=as.data.frame(sims)
    simulatedvec=s[,1]
    names(simulatedvec)=nam
    top=names(sort(simulatedvec, decreasing = TRUE))[1:tips]
    insum=countInternal(tree, leafBitMaps, top)
    insum
    if (insum==internal) {break}
  }
  # plot(t)
  return(top)
}


permulationCorList = correlationList
convertPermulationFormat = function(permulationCorList, RERObject = RERObject, permulationNumber = permulationNumber){
  permulationCorList
  permPvals = data.frame(matrix(ncol = permulationNumber, nrow = nrow(RERObject)))
  rownames(permPvals) = rownames(RERObject)
  permRhovals = data.frame(matrix(ncol = permulationNumber, nrow = nrow(RERObject)))
  rownames(permRhovals) = rownames(RERObject)
  permStatvals = data.frame(matrix(ncol = permulationNumber, nrow = nrow(RERObject)))
  rownames(permStatvals) = rownames(RERObject)
  for (i in 1:length(permulationCorList)) {
    permPvals[, i] = permulationCorList[[i]]$P
    permRhovals[, i] = permulationCorList[[i]]$Rho
    permStatvals[, i] = sign(permulationCorList[[i]]$Rho) * -log10(permulationCorList[[i]]$P)
  }
  output = vector("list", 3)
  output[[1]] = permPvals
  output[[2]] = permRhovals
  output[[3]] = permStatvals
  names(output) = c("corP", "corRho", "corStat")
  output
}


fastSimBinPhenoVecReport=function(tree, tips=0, internal, phenvec = NULL, rm=NULL, leafBitMaps=NULL){
  insum=0
  top = character(0)
  if (!is.null(phenvec)) {
    tips = sum(phenvec)
    if(is.null(rm)) {
      rm = ratematrix(tree, phenvec)
      message(rm)
    }
  } else if (is.null(rm)){
    print("No rate matrix supplied & cannot compute one")
    return()
  }
  
  if(is.null(leafBitMaps)){
    leafBitMaps = makeLeafMap(tree)
  }
  loop=0
  insumList = vector()
  repeat { #implements a do-while loop
    t=tree 
    #root.phylo(tree, root, resolve.root = T) #For now, let's use the ancestral root
    sims=sim.char(t, rm, nsim = 1)
    nam=rownames(sims)
    s=as.data.frame(sims)
    simulatedvec=s[,1]
    names(simulatedvec)=nam
    top=names(sort(simulatedvec, decreasing = TRUE))[1:tips]
    insum=countInternal(tree, leafBitMaps, top)
    loop = loop+1
    message("loop: ", loop)
    
    message("insum: ",insum)
    #insumList = append(insumList, insum)
    
    if (insum==internal) {break}
  }
  # plot(t)
  return(top)
}
phenotypeVectorSaver = phenotypeVector
phenotypeVector = phenotypeVectorSaver
phenotypeVector = 1

oldLauraPhenv = readRDS("Data/lauraphenvec.rds")
names(phenotypeVector) %in% names(oldLauraPhenv)
newForeground = names(phenotypeVector[phenotypeVector ==1])
oldForeground = names(oldLauraPhenv[oldLauraPhenv==1])

newBackground = names(phenotypeVector[phenotypeVector ==0])
oldBackground = names(oldLauraPhenv[oldLauraPhenv==0])

lengthDifference = length(phenotypeVector) - length(oldLauraPhenv)

fgOverlap = newForeground %in% oldForeground
all.equal(newForeground[order(newForeground)], oldForeground[order(oldForeground)])

bgOverlap = newBackground %in% oldBackground
bgOverlap
bgNewOnly = which(!newBackground %in% oldBackground)
newOnly = which(!names(phenotypeVector) %in% names(oldLauraPhenv))
length(newOnly)

trimmedNewPhenotypeVector = phenotypeVector[-newOnly]
length(trimmedNewPhenotypeVector)
length(oldLauraPhenv)
sortedTrimmedPhenVec = trimmedNewPhenotypeVector[order(trimmedNewPhenotypeVector)]

oldNames = names(oldLauraPhenv)
newNames = names(trimmedNewPhenotypeVector)

reOrder = match(newNames, oldNames)
reorderedNames = newNames[order(reOrder)]


orderedTrimmedPhenVec = trimmedNewPhenotypeVector[order(reOrder)]


phenotypeVector = trimmedNewPhenotypeVector

phenotypeVector = oldLauraPhenv

phenotypeVector = sortedTrimmedPhenVec

phenotypeVector = orderedTrimmedPhenVec





all.equal(orderedTrimmedPhenVec, oldLauraPhenv)
all.equal(phenotypeVector, oldLauraPhenv)

oldLauraTree = readRDS("Data/oldLauraTree.rds")
oldLauraMasterTree = readRDS("Data/oldLauraMasterTree.rds")
rootedMasterTree = oldLauraMasterTree

prunedMasterTree = pruneTree(rootedMasterTree, names(phenotypeVector))
doublePrunedMasterTree = pruneTree(rootedMasterTree, names(orderedTrimmedPhenVec))
all.equal(prunedMasterTree, oldLauraMasterTree)
rootedMasterTree = prunedMasterTree
internalNumber = 29




dev.new()
par(mfrow=c(1,3))

plot(oldLauraMasterTree, main = "oldLaura")

plot(doublePrunedMasterTree, main = "DoublePrune")
plot(prunedMasterTree, main = "SinglePrune")
inputTree = oldLauraMasterTree
testplot2 = plotTreeHighlightBranches(inputTree,hlspecies=which(inputTree$edge.length== 3),hlcols="blue", main="oldLauraMasterTree"); edgelabels(cex = 0.7, frame="none", font=2, adj=c(0,-0.2), col="blue"); nodelabels(cex = 0.7, frame="none", font=2, adj=c(-0.2,0.3), col="dark green"); tiplabels(cex = 0.8, frame="none", font=2, adj=c(0.2,0), col="dark red")
inputTree = initalRootedMasterTree
testplot2 = plotTreeHighlightBranches(inputTree,hlspecies=which(inputTree$edge.length== 3),hlcols="blue", main="initalRootedMasterTree"); edgelabels(cex = 0.7, frame="none", font=2, adj=c(0,-0.2), col="blue"); nodelabels(cex = 0.7, frame="none", font=2, adj=c(-0.2,0.3), col="dark green"); tiplabels(cex = 0.8, frame="none", font=2, adj=c(0.2,0), col="dark red")
inputTree = prunedMasterTree
testplot2 = plotTreeHighlightBranches(inputTree,hlspecies=which(inputTree$edge.length== 3),hlcols="blue", main="SinglePrune"); edgelabels(cex = 0.7, frame="none", font=2, adj=c(0,-0.2), col="blue"); nodelabels(cex = 0.7, frame="none", font=2, adj=c(-0.2,0.3), col="dark green"); tiplabels(cex = 0.8, frame="none", font=2, adj=c(0.2,0), col="dark red")

countInternalReport=function(t, bitMaps, fg){
  nInternal = t$Nnode
  fgNums = which(t$tip.label %in% fg)
  fgBitMap = makeBitMap(fgNums, nInternal + 1)
  count = 0
  internalList = vector()
  for (i in 1:nInternal) {
    if (compare(bitMaps[,i], fgBitMap)) {
      count = count +1
      internalList = append(internalList, paste(i, "---"))
    }
  }
  message(internalList)
  count
}

internalReport = countInternalReport(rootedMasterTree, bitMap,fg=names(phenotypeVector)[which(phenotypeVector==1)])
dev.new()
par(mfrow=c(1,3))
inputTree = rootedBinaryTree
testplot2 = plotTreeHighlightBranches(inputTree,hlspecies=which(rootedBinaryTree$edge.length== 1),hlcols="blue", main="RootedBinary"); edgelabels(cex = 0.7, frame="none", font=2, adj=c(0,-0.2), col="blue"); nodelabels(cex = 0.7, frame="none", font=2, adj=c(-0.2,0.3), col="dark green"); tiplabels(cex = 0.8, frame="none", font=2, adj=c(0.2,0), col="dark red")
inputTree = binaryTree
testplot2 = plotTreeHighlightBranches(inputTree,hlspecies=which(binaryTree$edge.length== 1),hlcols="blue", main="Binary"); edgelabels(cex = 0.7, frame="none", font=2, adj=c(0,-0.2), col="blue"); nodelabels(cex = 0.7, frame="none", font=2, adj=c(-0.2,0.3), col="dark green"); tiplabels(cex = 0.8, frame="none", font=2, adj=c(0.2,0), col="dark red")
inputTree = prunedBinaryTree
testplot2 = plotTreeHighlightBranches(inputTree,hlspecies=which(prunedBinaryTree$edge.length== 1),hlcols="blue", main="PrunedBinary"); edgelabels(cex = 0.7, frame="none", font=2, adj=c(0,-0.2), col="blue"); nodelabels(cex = 0.7, frame="none", font=2, adj=c(-0.2,0.3), col="dark green"); tiplabels(cex = 0.8, frame="none", font=2, adj=c(0.2,0), col="dark red")


inputTree$edge.length

binaryTree = inputTree

prunedBinaryTree = pruneTree(binaryTree, names(phenotypeVector))

testTree = mainTrees$masterTree

rootedBinaryTree = multi2di(binaryTree)

testBitMap = makeLeafMap(rootedBinaryTree)
internalReport = countInternalReport(rootedBinaryTree, testBitMap,fg=names(phenotypeVector)[which(phenotypeVector==1)])

permualationsDataFileName = paste(outputFolderName, filePrefix, "PermulationsData", runInstanceValue, ".rds", sep= "")
correctFormatPermulation = readRDS("Data/allInsectivoryPermulationsData1.rds")
correctFormatPermulation$corP[,1:4]
convertedPermulations$corP
permulationPValues = permpvalcor(cladesCorrelation, convertedPermulations)





# ------------ INSEPCT permPvalCor ------------------------------

realcor = cladesCorrelation
permvals = combinedPermulationsData

function (realcor, permvals) 
{
  permcor = permvals$corRho
  realstat = realcor$Rho
  names(realstat) = rownames(realcor)
  permcor = permcor[match(names(realstat), rownames(permcor)), 
  ]
  permpvals = vector(length = length(realstat))
  names(permpvals) = names(realstat)
  count = 1
  while (count <= length(realstat)) {
    if (is.na(realstat[count])) {
      permpvals[count] = NA
    }
    else {
      num = sum(abs(permcor[count, ]) > abs(realstat[count]), 
                na.rm = T)
      denom = sum(!is.na(permcor[count, ]))
      permpvals[count] = num/denom
    }
    count = count + 1
  }
  permpvals
}


# Comment 

function (realcor, permvals) 
{
  permcor = permvals$corRho                                                     #Rho values from the permulations
  realstat = realcor$Rho                                                        #Rho values from the real data 
  names(realstat) = rownames(realcor)                                           #Name the values after the row names 
  permcor = permcor[match(names(realstat), rownames(permcor)),]                 #trim the permulations' Rho values to only ones in the real data 
  permpvals = vector(length = length(realstat))                                 #Make a vector to store the pValues of a length equal to the number of realstat genes
  names(permpvals) = names(realstat)                                            #Copy the names from the real Rho values onto that vector
  count = 1                                                                     #Set a count equal to one 
  
  while (count <= length(realstat)) {                                           #While the count hasn't done all of the entries                                         
    if (is.na(realstat[count])) {                                               #if the real correlation is NA
      permpvals[count] = NA                                                     #Make the permulation value NA
    }
    else {                                                                      #If not 
      num = sum(abs(permcor[count, ]) > abs(realstat[count]),                   #Make a numerator which is the sum of the permulated values greater than the actual values; after removing NA entires. 
                na.rm = T)
      denom = sum(!is.na(permcor[count, ]))                                     #Make a denominator which is the sum of the non-NA permulation values
      permpvals[count] = num/denom                                              #p-value = numerator/denominator
    }
    count = count + 1                                                           #increase count by one
  }
  permpvals                                                                     #Return the vector of pValues 
}

# Rewrite

slowmatchPVal = function (realcor, permvals) 
{
  permcor = permvals$corRho                                                     #Rho values from the permulations
  realstat = realcor$Rho                                                        #Rho values from the real data 
  names(realstat) = rownames(realcor)                                           #Name the values after the row names 
  permcor = permcor[match(names(realstat), rownames(permcor)),]                 #trim the permulations' Rho values to only ones in the real data 
  permpvals = vector(length = length(realstat))                                 #Make a vector to store the pValues of a length equal to the number of realstat genes
  names(permpvals) = names(realstat)                                            #Copy the names from the real Rho values onto that vector
  count = 1                                                                     #Set a count equal to one 
  
  while (count <= length(realstat)) {                                           #While the count hasn't done all of the entries                                         
    if (is.na(realstat[count])) {                                               #if the real correlation is NA
      permpvals[count] = NA                                                     #Make the permulation value NA
    }
    else {                                                                      #If not 
      num = sum(abs(permcor[count, ]) > abs(realstat[count]),                   #Make a numerator which is the sum of the permulated values greater than the actual values; after removing NA entires. 
                na.rm = T)
      denom = sum(!is.na(permcor[count, ]))                                     #Make a denominator which is the sum of the non-NA permulation values
      permpvals[count] = num/denom                                              #p-value = numerator/denominator
    }
    count = count + 1                                                           #increase count by one
  }
  permpvals                                                                     #Return the vector of pValues 
}

permPValCorReport = function (realcor, permvals) {
  timeStart = Sys.time()
  permcor = permvals$corRho                                                     #Rho values from the permulations
  timePermExtractEnd = Sys.time()
  message("Time to extract permulation values: ", timePermExtractEnd - timeStart, attr(timePermExtractEnd - timeStart, "units"))

  realstat = realcor$Rho                                                        #Rho values from the real data 
  timeStatExtractEnd = Sys.time()
  message("Time to extract real values: ", timeStatExtractEnd - timePermExtractEnd, attr(timeStatExtractEnd - timePermExtractEnd, "units"))
  names(realstat) = rownames(realcor)                                           #Name the values after the row names 
  
  permcor = permcor[fmatch(names(realstat), rownames(permcor)),]                 #trim the permulations' Rho values to only ones in the real data 
  timeMatchEnd = Sys.time()
  message("Time to match names: ", timeMatchEnd - timeStatExtractEnd, attr(timeMatchEnd - timeStatExtractEnd, "units"))
  
  permpvals = vector(length = length(realstat))                                 #Make a vector to store the pValues of a length equal to the number of realstat genes
  names(permpvals) = names(realstat)                                            #Copy the names from the real Rho values onto that vector
  count = 1                                                                     #Set a count equal to one 
  
  timeBeforeLoopStart = Sys.time()
  message("Total time before loop start: ", timeBeforeLoopStart - timeStart, attr(timeBeforeLoopStart - timeStart, "units"))
  
  while (count <= length(realstat)) {                                           #While the count hasn't done all of the entries                                         
    timeLoopBegin = Sys.time()
    message("Gene number: ", count)
    if (is.na(realstat[count])) {                                               #if the real correlation is NA
      permpvals[count] = NA                                                     #Make the permulation value NA
      message("is NA")
    }
    else {                                                                      #If not 
      num = sum(abs(permcor[count, ]) > abs(realstat[count]),                   #Make a numerator which is the sum of the permulated values greater than the actual values; after removing NA entires. 
                na.rm = T)
      timeNumeratorCalc = Sys.time()
      message("Numerator Calculate time: ", timeNumeratorCalc - timeLoopBegin, attr(timeNumeratorCalc - timeLoopBegin, "units"))
     
       denom = sum(!is.na(permcor[count, ]))                                     #Make a denominator which is the sum of the non-NA permulation values
      timeDenomSum = Sys.time()
      message("Denominator sum time: ", timeDenomSum - timeNumeratorCalc, attr(timeDenomSum - timeNumeratorCalc, "units"))
       permpvals[count] = num/denom                                              #p-value = numerator/denominator
    }
    timeLoopEnd = Sys.time()
    message("Total loop time: ", timeLoopEnd - timeLoopBegin, attr(timeLoopEnd - timeLoopBegin, "units"))
    count = count + 1                                                           #increase count by one
  }
  message("calculation compelte.")
  timeEnd = Sys.time()
  message("Total calculation time: ", timeEnd - timeStart, attr(timeEnd - timeStart, "units"))
  permpvals                                                                     #Return the vector of pValues 
}













mark(
  slowmatchPVal(realcor, permvals),
  fastmatchPVal(realcor, permvals), iterations = 10
)[c("expression", "min", "median", "itr/sec", "n_gc")]

testTime1 = Sys.time()
slowmatchPVal(realcor, permvals)
testTime2 = Sys.time()
timeSpent = testTime2 - testTime1
timeSpent

testTime1 = Sys.time()
fastmatchPVal(realcor, permvals)
testTime2 = Sys.time()
timeSpent = testTime2 - testTime1
timeSpent

# Convert to apply 

fastmatchPVal = function (realcor, permvals) 
{
  permcor = permvals$corRho                                                     #Rho values from the permulations
  realstat = realcor$Rho                                                        #Rho values from the real data 
  names(realstat) = rownames(realcor)                                           #Name the values after the row names 
  permcor = permcor[fmatch(names(realstat), rownames(permcor)),]                 #trim the permulations' Rho values to only ones in the real data 
  permpvals = vector(length = length(realstat))                                 #Make a vector to store the pValues of a length equal to the number of realstat genes
  names(permpvals) = names(realstat)                                            #Copy the names from the real Rho values onto that vector
  count = 1                                                                     #Set a count equal to one 
  
  calcP = function(realValue){
    currentName = names(realValue)
    currentName
    
  }
  test = vapply(realstat, calcP)
  
  while (count <= length(realstat)) {                                           #While the count hasn't done all of the entries                                         
    if (is.na(realstat[count])) {                                               #if the real correlation is NA
      permpvals[count] = NA                                                     #Make the permulation value NA
    }
    else {                                                                      #If not 
      num = sum(abs(permcor[count, ]) > abs(realstat[count]),                   #Make a numerator which is the sum of the permulated values greater than the actual values; after removing NA entires. 
                na.rm = T)
      denom = sum(!is.na(permcor[count, ]))                                     #Make a denominator which is the sum of the non-NA permulation values
      permpvals[count] = num/denom                                              #p-value = numerator/denominator
    }
    count = count + 1                                                           #increase count by one
  }
  permpvals                                                                     #Return the vector of pValues 
}

getP1sided <- function(corrow, permtib, geneind, rhoind) {
  #compute 1-sided p-value for a row of the correlation table
  #geneind = index of column containing gene name; rhoind = index of column containing Rho value
  ptcol <- permtib[,corrow[geneind]]
  if (is.na(as.numeric(corrow[rhoind])) == FALSE && sum(is.na(ptcol)==FALSE) > 0) {
    if (as.numeric(corrow[rhoind]) > 0) {
      p1s <- sum(ptcol >= as.numeric(corrow[rhoind]),na.rm=TRUE)/sum(is.na(ptcol)==FALSE)
    } else if (as.numeric(corrow[rhoind]) < 0) {
      p1s <- sum(ptcol <= as.numeric(corrow[rhoind]),na.rm=TRUE)/sum(is.na(ptcol)==FALSE)
    } else {
      p1s <- sum(ptcol != 0, na.rm=TRUE)/sum(is.na(ptcol)==FALSE)
    }
  } else {
    p1s <- NA
  }
  return(p1s)
}




fullPermData = readRDS("Data/allInsectivoryPermulationsData1.rds")
fullPermData2 = combinePermData(fullPermData, fullPermData, enrich = F)
fullPermData4 = combinePermData(fullPermData2, fullPermData2, enrich = F)
fullPermData8 = combinePermData(fullPermData4, fullPermData4, enrich = F)
fullPermData16 = combinePermData(fullPermData8, fullPermData8, enrich = F)

rm(fullPermData8)
rm(fullPermData4)
rm(fullPermData2)

realcor = cladesCorrelation
permvals = combinedPermulationsData
permvals = fullPermData
permvals = fullPermData2
permvals = fullPermData4
permvals = fullPermData8
permvals = fullPermData16

#permPValCorReport = function (realcor, permvals) 
{
  timeStart = Sys.time()
  permcor = fullPermData16$corRho                                                     #Rho values from the permulations
  timePermExtractEnd = Sys.time()
  message("Time to extract permulation values: ", timePermExtractEnd - timeStart, attr(timePermExtractEnd - timeStart, "units"))
  
  realstat = realcor$Rho                                                        #Rho values from the real data 
  timeStatExtractEnd = Sys.time()
  message("Time to extract real values: ", timeStatExtractEnd - timePermExtractEnd, attr(timeStatExtractEnd - timePermExtractEnd, "units"))
  names(realstat) = rownames(realcor)                                           #Name the values after the row names 
  
  permcor = permcor[match(names(realstat), rownames(permcor)),]                 #trim the permulations' Rho values to only ones in the real data 
  timeMatchEnd = Sys.time()
  message("Time to match names: ", timeMatchEnd - timeStatExtractEnd, attr(timeMatchEnd - timeStatExtractEnd, "units"))
  
  permpvals = vector(length = length(realstat))                                 #Make a vector to store the pValues of a length equal to the number of realstat genes
  names(permpvals) = names(realstat)                                            #Copy the names from the real Rho values onto that vector
  count = 2                                                                     #Set a count equal to one 
  
  timeBeforeLoopStart = Sys.time()
  message("Total time before loop start: ", timeBeforeLoopStart - timeStart, attr(timeBeforeLoopStart - timeStart, "units"))
  
  #while (count <= length(realstat)) {                                           #While the count hasn't done all of the entries                                         
    timeLoopBegin = Sys.time()
    message("Gene number: ", count)
    if (is.na(realstat[count])) {                                               #if the real correlation is NA
      permpvals[count] = NA                                                     #Make the permulation value NA
      message("is NA")
    }else {                                                                      #If not 
      permRow = permcor[count, ]
      permRow = as.vector(permcor[count, ])
      permRow = t(permRow)
      permRowMakeTime = Sys.time()
      message("Time to make permulation row: ", permRowMakeTime - timeLoopBegin, attr(permRowMakeTime - timeLoopBegin, "units"))
      
      realRow = realstat[count]
      realRowMakeTime = Sys.time()
      message("Time to make real row: ", realRowMakeTime - permRowMakeTime, attr(realRowMakeTime - permRowMakeTime, "units"))
      
      num = sum(abs(permRow) > abs(realRow),                   #Make a numerator which is the sum of the permulated values greater than the actual values; after removing NA entires. 
                na.rm = T)
      timeNumeratorCalc = Sys.time()
      message("Numerator Calculate time: ", timeNumeratorCalc - realRowMakeTime, attr(timeNumeratorCalc - realRowMakeTime, "units"))
      
      denom = sum(!is.na(permRow))                                     #Make a denominator which is the sum of the non-NA permulation values
      timeDenomSum = Sys.time()
      message("Denominator sum time: ", timeDenomSum - timeNumeratorCalc, attr(timeDenomSum - timeNumeratorCalc, "units"))
      permpvals[count] = num/denom                                              #p-value = numerator/denominator
    }
    timeLoopEnd = Sys.time()
    message("Total loop time: ", timeLoopEnd - timeLoopBegin, attr(timeLoopEnd - timeLoopBegin, "units"))
    count = count + 1                                                           #increase count by one
  #}
  message("calculation compelte.")
  timeEnd = Sys.time()
  message("Total calculation time: ", timeEnd - timeStart, attr(timeEnd - timeStart, "units"))
  #permpvals                                                                     #Return the vector of pValues 
}


permRow = permcor[count, ]
permRowVec = as.vector(permcor[count, ])
permRowCollumn = t(permRow)
permRowVec2 = as.vector(permRowCollumn)
realRow = realstat[count]

numCalcVec = function(permRowType, realVal = realRow){
  num = sum(abs(permRowType) > abs(realVal),                   #Make a numerator which is the sum of the permulated values greater than the actual values; after removing NA entires. 
            na.rm = T)
}

numCalcCol = function(permRowType, realVal = realRow){
  num = sum(abs(permRowType) > abs(realVal),                   #Make a numerator which is the sum of the permulated values greater than the actual values; after removing NA entires. 
            na.rm = T)
}
numCalcRow = function(permRowType, realVal = realRow){
  num = sum(abs(permRowType) > abs(realVal),                   #Make a numerator which is the sum of the permulated values greater than the actual values; after removing NA entires. 
            na.rm = T)
}

mark(
  numCalcVec(permRowVec2, realRow),
  numCalcCol(permRowCollumn, realRow), 
  numCalcRow(permRow, realRow), iterations = 10
  
)[c("expression", "min", "median", "itr/sec", "n_gc")]

mark(
  numCalcRow(permRow, realRow),
  numCalcCol(permRowCollumn, realRow),  iterations = 10
  
)[c("expression", "min", "median", "itr/sec", "n_gc")]

mark(
  numCalcVec(permRowVec2, realRow),
  numCalcCol(permRowCollumn, realRow), iterations = 10
  
)[c("expression", "min", "median", "itr/sec", "n_gc")]


num = sum(abs(permRowCollumn) > abs(realRow),                   #Make a numerator which is the sum of the permulated values greater than the actual values; after removing NA entires. 
          na.rm = T)


testOutputFormat = permPValCorReport()
oldOutputFormat = permulationPValues

all.equal(oldOutputFormat, permulationPValues)


#Time to make permulation row: 0.857425212860107secs
#Time to make real row: 0.000335931777954102secs
#Numerator Calculate time: 0.00166702270507812secs
#Denominator sum time: 0.000598907470703125secs

#Time to make permulation row: 1.34229707717896secs
#Time to make real row: 0.000356912612915039secs
#Numerator Calculate time: 0.0234830379486084secs
#Denominator sum time: 0.000855922698974609secs
pathCladesObjectUnnamed = unname(pathCladesObject)
length(pathsObject[pathsObject ==1])
length(pathCladesObject[pathCladesObject ==1])



correl2 = correl

all.equal(correl2, correl)

RERObject2 = RERObject
saveRDS(RERObject2, file = "Results/carnvHerbRERFileNewGeneration.rds")


speciesFilterOld = speciesFilter
speciesFilter == speciesFilter[speciesFilter %in% speciesFilterOld]


all.equal(RERObject, RERObject2)
length(RERObject)
length(which(is.na(RERObject)))/length(RERObject)
length(which(is.na(RERObject2)))/length(RERObject2)

REROldSpeciesList = readRDS("Output/carnvHerbs/carnvHerbsRERFile.rds")
all.equal(REROldSpeciesList, RERObject)
length(which(is.na(REROldSpeciesList)))
RERNewSpeciesList = RERObject2
length(which(is.na(RERObject2)))

RERObject = RERNewSpeciesList

all.equal(RERObject, REROldSpeciesList)

length(which(is.na(correl)))

rm(foreground2Tree)
foreground2Tree =function (foreground, treesObj, plotTree = T, clade = c("ancestral", 
                                                                         "terminal", "all"), weighted = F, transition = "unidirectional", 
                           useSpecies = NULL) 
{
  clade <- match.arg(clade)
  res = treesObj$masterTree
  if (!is.null(useSpecies)) {
    sp.miss = setdiff(res$tip.label, useSpecies)
    if (length(sp.miss) > 0) {
      #message(paste0("Species from master tree not present in useSpecies: ", 
      #               paste(sp.miss, collapse = ",")))
    }
    useSpecies = intersect(useSpecies, res$tip.label)
    res = pruneTree(res, useSpecies)
  }
  else {
    useSpecies = res$tip.label
  }
  foreground = intersect(foreground, useSpecies)
  res$edge.length <- rep(0, length(res$edge.length))
  if (clade == "terminal") {
    res$edge.length[nameEdges(res) %in% foreground] = 1
    names(res$edge.length) = nameEdges(res)
  }
  else if (clade == "ancestral") {
    weighted = F
    if (transition == "bidirectional") {
      res <- inferBidirectionalForegroundClades(res, foreground, 
                                                ancestralOnly = T)
    }
    else {
      res <- inferUnidirectionalForegroundClades(res, foreground, 
                                                 ancestralOnly = T)
    }
  }
  else {
    if (transition == "bidirectional") {
      res <- inferBidirectionalForegroundClades(res, foreground, 
                                                ancestralOnly = F)
    }
    else {
      res <- inferUnidirectionalForegroundClades(res, foreground, 
                                                 ancestralOnly = F)
    }
  }
  if (weighted) {
    if (clade == "all") {
      tobeweighted <- rep(TRUE, length(res$edge.length))
      tobeweighted[res$edge.length == 0] <- FALSE
      while (sum(tobeweighted) > 0) {
        edgetodo <- which(tobeweighted == T)[1]
        clade.down.edges = getAllCladeEdges(res, edgetodo)
        if (length(clade.down.edges) > 1) {
          clade.edges = c(clade.down.edges, edgetodo)
          clade.edges.toweight <- clade.edges[res$edge.length[clade.edges] == 
                                                1]
          res$edge.length[clade.edges.toweight] <- 1/(length(clade.edges.toweight))
          tobeweighted[clade.edges] <- FALSE
        }
        else {
          tobeweighted[clade.down.edges] <- FALSE
        }
      }
    }
    else if (clade == "terminal") {
      tobeweightededgeterminalnode <- unique(res$edge[(res$edge[, 
                                                                2] %in% c(1:length(res$tip.label))), 1])
      tobeweighted <- setdiff(match(tobeweightededgeterminalnode, 
                                    res$edge[, 2]), NA)
      for (edgetodo in tobeweighted) {
        clade.down.edges = getAllCladeEdges(res, edgetodo)
        if (all(res$edge.length[clade.down.edges] == 
                1)) {
          res$edge.length[clade.down.edges] <- 0.5
        }
      }
    }
  }
  if (plotTree) {
    res2 = res
    mm = min(res2$edge.length[res2$edge.length > 0])
    res2$edge.length[res2$edge.length == 0] = max(0.02, mm/20)
    plot(res2, main = paste0("Clade: ", clade, "\nTransition: ", 
                             transition, "\nWeighted: ", weighted), cex = 0.5)
    if (weighted) {
      labs <- round(res$edge.length, 3)
      labs[labs == 0] <- NA
      edgelabels(labs, col = "black", bg = "transparent", 
                 adj = c(0.5, -0.5), cex = 0.4, frame = "n")
    }
  }
  res
}

debugTruncatedFG2Tree = function (foreground, treesObj, plotTree = T, clade = c("ancestral", "terminal", "all"), weighted = F, transition = "unidirectional", useSpecies = NULL) 
{
  clade <- match.arg(clade)
  res = treesObj$masterTree
  if (!is.null(useSpecies)) {
    sp.miss = setdiff(res$tip.label, useSpecies)
    if (length(sp.miss) > 0) {
     # message(paste0("Species from master tree not present in useSpecies: ", 
    #               paste(sp.miss, collapse = ",")))
    }
    useSpecies = intersect(useSpecies, res$tip.label)
    res = pruneTree(res, useSpecies)
  }
  foreground = intersect(foreground, useSpecies)
  #message("fed in foreground: ", foreground)
  if(length(foreground) == 0){message("No foreground in species filter")}
  res$edge.length <- rep(0, length(res$edge.length))
  {
    if (transition == "bidirectional") {
      res <- debugInferBidirectionalForegroundClades(res, foreground, 
                                                ancestralOnly = F)
    }
  }
  res
}
####
useSpecies = newSpeciesFilter
foreground = brokenPermulatedForegroundExample 


####
treesObj = mainTrees
res = treesObj$masterTree
treeinput  = res
foreground = permulatedForeground
ancestralOnly = F

debugInferBidirectionalForegroundClades <- function(treeinput, foreground = NULL, ancestralOnly = F){
  tree <- treeinput
  tip.vals=rep(0, length(tree$tip.label))
  names(tip.vals)=tree$tip.label
  tip.vals[foreground]=1
  tmp=cbind(as.character(tip.vals))
  rownames(tmp)=names(tip.vals)
  tip.vals=tmp
  #Add option to function for "type" within ancestral.pars
  ancres=ancestral.pars(tree, df<-as.phyDat(tip.vals, type="USER", levels=unique(as.character(tip.vals))),type="ACCTRAN" )
  ancres=unlist(lapply(ancres, function(x){x[2]}))
  internalVals=ancres
  #evals=matrix(nrow=nrow(treesObj$masterTree$edge), ncol=2)
  evals=matrix(nrow=nrow(tree$edge), ncol=2)
  eres=ancres
  #evals[,1]=eres[treesObj$masterTree$edge[,1]]
  evals[,1]=eres[tree$edge[,1]]
  #evals[,2]=eres[treesObj$masterTree$edge[,2]]
  evals[,2]=eres[tree$edge[,2]]
  tree$edge.length=evals[,2]-evals[,1]
  #res$edge.length[res$edge.length<1]=0
  if(!ancestralOnly){
    edgeIndex=which(tree$edge.length>0)
    edgeIndexNeg=which(tree$edge.length<0)
    edgeIndexAll = c(edgeIndex,edgeIndexNeg)
    edgeDirection = c(rep(1, length(edgeIndex)),rep(-1, length(edgeIndexNeg)))
    edgedf = data.frame(edgeIndexAll,edgeDirection)
    edgedf = edgedf[order(edgedf$edgeIndexAll),]
    clade.edges=NA
    clade.lengths=NA
    cladedf = data.frame(clade.edges,clade.lengths)
    for(i in 1:nrow(edgedf)) { #Does this go from ancestral to terminal?
      #save the clade until the edges no longer overlap
      #message("ancEdge Input: ", edgedf$edgeIndexAll[i])
      clade.edges=debugGetAllCladeEdges(tree, edgedf$edgeIndexAll[i])
      
      clade.edges=unique(c(edgedf$edgeIndexAll[i], clade.edges))
      if (any(clade.edges %in% cladedf$clade.edges)==F) {
        tree$edge.length[cladedf$clade.edges[which(cladedf$clade.lengths==1)]]=1
        if (edgedf$edgeDirection[i] == 1) {
          clade.lengths = c(rep(1,length(clade.edges)))
        } else {
          clade.lengths = c(rep(0,length(clade.edges)))
        }
        cladedf = data.frame(clade.edges,clade.lengths)
      } else {
        #update df lengths
        if (edgedf$edgeDirection[i] == 1) {
          cladedf$clade.lengths[which(cladedf$clade.edges %in% clade.edges)] = 1
        } else {
          cladedf$clade.lengths[which(cladedf$clade.edges %in% clade.edges)] = 0
        }
      }
    }
    #update edge lengths from the final clade
    tree$edge.length[cladedf$clade.edges[which(cladedf$clade.lengths==1)]]=1
    tree$edge.length[tree$edge.length<0]=0
  }
  tree$edge.length[tree$edge.length<0]=0
  tree
}




debugGetAllCladeEdges=function(tree, AncEdge){
  node=tree$edge[AncEdge,2]
  #message(tree$edge)
  #message("AncEdge: ", AncEdge)
  #message("Tree$edge[AncEdge]: ", tree$edge[AncEdge,])
  #message("Tree$edge[AncEdge,2]: ", tree$edge[AncEdge,2])
  #message("getClades node: ",node)
  #get descendants
  iid=debugGetDescendants(tree, node)
  #find their edges
  iim=match(iid, tree$edge[,2])
  iim
}

debugGetDescendants =function (tree, node, curr = NULL) 
{
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of class \"phylo\".")
  if (is.null(curr)) 
    curr <- vector()
  daughters <- tree$edge[which(tree$edge[, 1] == node), 2]
  curr <- c(curr, daughters)
  #message(length(curr) == 0)
  #message(node <= Ntip(tree))
  if(is.na(node <= Ntip(tree))){
    message("node: ",node)
    message("nTipTree: ", Ntip(tree))
  }
  if (length(curr) == 0 && node <= Ntip(tree)) 
    curr <- node
  w <- which(daughters > Ntip(tree))
  if (length(w) > 0) 
    for (i in 1:length(w)) curr <- getDescendants(tree, daughters[w[i]], 
                                                  curr)
  return(curr)
}



for(i in 1){
  permulatedTree = debugTruncatedFG2Tree(brokenPermulatedForegroundExample, mainTrees, plotTree=F, clade="all", transition="bidirectional", useSpecies=speciesFilter) #generate a tree using that foregound
  message(i)
}

for(i in 1){
  permulatedTree = debugTruncatedFG2Tree(brokenPermulatedForegroundExample, mainTrees, plotTree=F, clade="all", transition="bidirectional", useSpecies=oldSpeciesFilter) #generate a tree using that foregound
  message(i)
}

#(length(curr) == 0 && node <= Ntip(tree)) curr <- node : 


#IMport old speciesfilter: 
speciesFilter = readRDS("Data/carnvHerbsSpeciesFilterFullTreeWithNonlaura.rds")
oldSpeciesFilter = readRDS("Data/carnvHerbsSpeciesFilterFullTreeWithNonlaura.rds")

#import new speciesfilter: 
speciesFilter = readRDS(speciesFilterFileName)
newSpeciesFilter = readRDS(speciesFilterFileName)


all.equal(oldSpeciesFilter, newSpeciesFilter)

oldSpeciesFilterTrimmed = oldSpeciesFilter[which(oldSpeciesFilter %in% newSpeciesFilter)]
length(oldSpeciesFilterTrimmed)

all.equal(oldSpeciesFilterTrimmed, speciesFilter)

oldSpeciesFilterTrimmedSorted = oldSpeciesFilterTrimmed[order(oldSpeciesFilterTrimmed)]
newSpeciesFilterSorted = newSpeciesFilter[order(newSpeciesFilter)]

all.equal(oldSpeciesFilterTrimmedSorted, newSpeciesFilterSorted)

oldSpeciesFilterUnTrimmedSorted = oldSpeciesFilter[order(oldSpeciesFilter)]

oldSpeciesFilterCut = oldSpeciesFilter[1:71]

oldSpeciesFilterRemoveFew = oldSpeciesFilter[1:90]

oldSpeciesFilter100 = oldSpeciesFilter[1:100]
oldSpeciesFilter110= oldSpeciesFilter[1:110]

for(i in 1:100){
  singlePermCorrelation = computeCorrelationOnePermulationDebug(rootedMasterTree, phenotypeVector, mainTrees, RERObject)
  message(i)
}

for(i in 1:100){
  permulatedForeground = fastSimBinPhenoVec(tree=rootedMasterTree, phenvec=phenotypeVector, internal=internalNumber)
  inferBidirectionalForegroundClades(treeinput, permulatedForeground, F)
  message(i)
}

for(i in 1:100){
  permulatedForeground = fastSimBinPhenoVec(tree=rootedMasterTree, phenvec=phenotypeVector, internal=internalNumber)
  permulatedTree = debugTruncatedFG2Tree(permulatedForeground, mainTrees, plotTree=F, clade="all", transition="bidirectional") #generate a tree using that foregound
  message(i)
}

for(i in 1:100){
  permulatedForeground = fastSimBinPhenoVec(tree=rootedMasterTree, phenvec=phenotypeVector, internal=internalNumber)
  permulatedTree = debugTruncatedFG2Tree(permulatedForeground, mainTrees, plotTree=F, clade="all", transition="bidirectional", useSpecies=oldSpeciesFilter) #generate a tree using that foregound
  message(i)
}

for(i in 1:100){
  permulatedForeground = fastSimBinPhenoVec(tree=rootedMasterTree, phenvec=phenotypeVector, internal=internalNumber)
  permulatedTree = debugTruncatedFG2Tree(permulatedForeground, mainTrees, plotTree=F, clade="all", transition="bidirectional", useSpecies=oldSpeciesFilterUnTrimmedSorted) #generate a tree using that foregound
  message(i)
}

for(i in 1:100){
  permulatedForeground = fastSimBinPhenoVec(tree=rootedMasterTree, phenvec=phenotypeVector, internal=internalNumber)
  permulatedTree = debugTruncatedFG2Tree(permulatedForeground, mainTrees, plotTree=F, clade="all", transition="bidirectional", useSpecies=oldSpeciesFilterCut) #generate a tree using that foregound
  message(i)
}

for(i in 1:100){
  permulatedForeground = fastSimBinPhenoVec(tree=rootedMasterTree, phenvec=phenotypeVector, internal=internalNumber)
  permulatedTree = debugTruncatedFG2Tree(permulatedForeground, mainTrees, plotTree=F, clade="all", transition="bidirectional", useSpecies=oldSpeciesFilterRemoveFew) #generate a tree using that foregound
  message(i)
}

for(i in 1:100){
  permulatedForeground = fastSimBinPhenoVec(tree=rootedMasterTree, phenvec=phenotypeVector, internal=internalNumber)
  permulatedTree = debugTruncatedFG2Tree(permulatedForeground, mainTrees, plotTree=F, clade="all", transition="bidirectional", useSpecies=oldSpeciesFilter100) #generate a tree using that foregound
  message(i)
}

for(i in 1:100){
  permulatedTree = debugTruncatedFG2Tree(permulatedForeground, mainTrees, plotTree=F, clade="all", transition="bidirectional", useSpecies=newSpeciesFilter) #generate a tree using that foregound
  message(i)
}



for(i in 1:1000){
  permulatedForeground = fastSimBinPhenoVec(tree=rootedMasterTree, phenvec=phenotypeVector, internal=internalNumber)
  permulatedTree = debugTruncatedFG2Tree(permulatedForeground, mainTrees, plotTree=F, clade="all", transition="bidirectional", useSpecies=oldSpeciesFilter110) #generate a tree using that foregound
  message(i)
}

for(i in 1:100){
  permulatedForeground = fastSimBinPhenoVec(tree=rootedMasterTree, phenvec=phenotypeVector, internal=internalNumber)
  permulatedTree = debugTruncatedFG2Tree(permulatedForeground, mainTrees, plotTree=F, clade="all", transition="bidirectional", useSpecies=speciesFilter) #generate a tree using that foregound
  message(i)
}

for(i in 1:1000){
  permulatedForeground = fastSimBinPhenoVec(tree=rootedMasterTree, phenvec=phenotypeVector, internal=internalNumber)
  permulatedTree = debugTruncatedFG2Tree(permulatedForeground, mainTrees, plotTree=F, clade="all", transition="bidirectional", useSpecies=speciesFilter) #generate a tree using that foregound
  message(i)
}

for(i in 1:100){
  permulatedTree = debugTruncatedFG2Tree(permulatedForeground, mainTrees, plotTree=F, clade="all", transition="bidirectional", useSpecies=speciesFilter) #generate a tree using that foregound
  message(i)
}
#saveRDS(permulatedForeground, "brokenPermualtedForegroundExample.rds")
brokenPermulatedForegroundExample = readRDS("brokenPermualtedForegroundExample.rds")

for(i in 1:100){
  permulatedTree = debugTruncatedFG2Tree(brokenPermulatedForegroundExample, mainTrees, plotTree=F, clade="all", transition="bidirectional", useSpecies=speciesFilter) #generate a tree using that foregound
  message(i)
}

for(i in 1:100){
  permulatedTree = debugTruncatedFG2Tree(brokenPermulatedForegroundExample, mainTrees, plotTree=F, clade="all", transition="bidirectional", useSpecies=oldSpeciesFilter) #generate a tree using that foregound
  message(i)
}
brokenPermulatedForegroundExample %in% speciesFilter
brokenPermulatedForegroundExample %in% oldSpeciesFilter

for(i in 1:100){
  permulatedForeground = fastSimBinPhenoVec(tree=rootedMasterTree, phenvec=phenotypeVector, internal=internalNumber)
  message(permulatedForeground %in% speciesFilter)
  message(paste(permulatedForeground %in% oldSpeciesFilter))
  permulatedTree = debugTruncatedFG2Tree(permulatedForeground, mainTrees, plotTree=F, clade="all", transition="bidirectional", useSpecies=oldSpeciesFilter) #generate a tree using that foregound
  message(i)
}
permulatedForeground %in% rootedMasterTree$tip.label

for(i in 1:50){
  permulatedForeground = fastSimBinPhenoVec(tree=rootedMasterTree, phenvec=phenotypeVector, internal=internalNumber)
  fgInFilter = permulatedForeground[(permulatedForeground %in% oldSpeciesFilter)]
  fgInRemovedFilter = fgInFilter[!fgInFilter %in% newSpeciesFilter]
  fgInKeptFilter = fgInFilter[fgInFilter %in% newSpeciesFilter]
  message("fg length: ", length(permulatedForeground))
  message("fg in kept portion: ", length(fgInKeptFilter))
  message("fg only in removed portion: ", length(fgInRemovedFilter))
  permulatedTree = debugTruncatedFG2Tree(permulatedForeground, mainTrees, plotTree=F, clade="all", transition="bidirectional", useSpecies=oldSpeciesFilter) #generate a tree using that foregound
  message(i)
}

fgInFilter = permulatedForeground[(permulatedForeground %in% oldSpeciesFilter)]
fgInRemovedFilter = fgInFilter[!fgInFilter %in% newSpeciesFilter]
message("fg only in removed portion: ", length(fgInRemovedFilter))

fgInRemovedFilter %in% newSpeciesFilter
fgInRemovedFilter %in% oldSpeciesFilter

computeCorrelationOnePermulationDebug

computeCorrelationOnePermulationDebug = function(rootedMasterTree, phenotypeVector, mainTrees, RERObject, min.sp =35, internalNumber = internalNumberValue){
  #singlePermStartTime = Sys.time()
  message("perm")
  permulatedForeground = fastSimBinPhenoVec(tree=rootedMasterTree, phenvec=phenotypeVector, internal=internalNumber)                                     #generate a null foreground via permulation
  
  message("tree")
  #permulatedTree = foreground2Tree(permulatedForeground, mainTrees, plotTree=F, clade="all", transition="bidirectional", useSpecies=speciesFilter) #generate a tree using that foregound
  permulatedTree = debugTruncatedFG2Tree(permulatedForeground, mainTrees, plotTree=F, clade="all", transition="bidirectional", useSpecies=speciesFilter) #generate a tree using that foregound
  
  message("paths")
  #permulatedPaths = tree2Paths(permulatedTree, mainTrees, binarize=T, useSpecies=speciesFilter)                                                    #generate a path from that tree
  
  #singlePermulationEndTime = Sys.time()
  #permulationDuration = singlePermulationEndTime - singlePermStartTime
  #message("Permulation time: ", permulationDuration)
  #permulatedCorrelations = correlateWithBinaryPhenotype(RERObject, permulatedPaths, min.sp=min.sp)                                                 #Use that path to get a coreelation of the null phenotype to genes (this is the outbut of a Get PermsBinary run)
  
  #correlationEndTime = Sys.time()
  #correlationDuration = correlationEndTime - singlePermulationEndTime
  #message("Correlation Duration: ", correlationDuration)
  #permulatedCorrelations
}

oldSpeciesFilter %in% newSpeciesFilter


names(phenotypeVector) %in% speciesFilter

permulationPValues = permPValCorReport(cladesCorrelation, combinedPermulationsData, startNumber = startValue, geneNumber = geneNumberValue, plusOne = plusOneValue)

realcor = cladesCorrelation; permvals = combinedPermulationsData; startNumber=1; geneNumber = NA; report = TRUE; plusOne = FALSE


# ---- Compare clades vs non-clades phenotype trees ---

toytreefile = "subsetMammalGeneTrees.txt" 
mainTrees

#Calculate RERs
RERObject 

#Clades-based correlation results:

CladesFgTree = foreground2TreeClades(foregroundSpecies,sistersListExport,mainTrees,plotTree=T)

plotTree(CladesFgTree)

#Non-Clades-Based correlation results:
noCladesFgTree = foreground2Tree(foregroundSpecies,mainTrees,plotTree=T, clade = "all")
filteredNoCladesFgTree = foreground2Tree(foregroundSpecies,mainTrees,plotTree=T, clade = "all", useSpecies = speciesFilter)

length(NoCladesFgTree$tip.label)
#clades tree modified
FilteredCladesFgTree = foreground2TreeClades(foregroundSpecies,sistersListExport,mainTrees,plotTree=T, useSpecies = speciesFilter)
prunedCladesTree = pruneTree(CladesFgTree, speciesFilter)

CladesFgTree$tip.label
foregroundSpecies %in% speciesFilter
foregroundSpecies %in% CladesFgTree$tip.label
length(which(CladesFgTree$tip.label %in% speciesFilter))
length(speciesFilter)

par(mfrow = c(1,4))
plotTree(noCladesFgTree)
plotTree(CladesFgTree)
plotTree(prunedCladesTree) 
title("clades Based")
plotTree(filteredNoCladesFgTree)
title("non-clades Based")

prunedCladesTreeViewable = prunedCladesTree
prunedCladesTreeViewable$edge.length[1:length(prunedCladesTreeViewable$edge.length)] = 0.5
prunedCladesTreeViewable$edge.length[prunedCladesTree$edge.length == 1] = 3

filteredNoCladesFgTreeViewable = filteredNoCladesFgTree
filteredNoCladesFgTreeViewable$edge.length[1:length(filteredNoCladesFgTreeViewable$edge.length)] = 0.5
filteredNoCladesFgTreeViewable$edge.length[filteredNoCladesFgTree$edge.length == 1] = 3

par(mfrow = c(1,2))
plotTreeHighlightBranches(prunedCladesTreeViewable, hlspecies=which(prunedCladesTree$edge.length==1), hlcols="blue") 
title("clades Based")
plotTreeHighlightBranches(filteredNoCladesFgTreeViewable, hlspecies=which(filteredNoCladesFgTree$edge.length==1), hlcols="blue")
title("non-clades Based")


all.equal(noCladesFgTree, CladesFgTree)

plotTree(binaryPhenotypeTree)
plotTree(FilteredCladesFgTree)





# --- toy examples of clades vs non clades 
library(RERconverge)
rerpath = find.package('RERconverge')
par(mfrow = c(1,2))
#read trees
toytreefile = "subsetMammalGeneTrees.txt" 
toyTrees=readTrees(paste(rerpath,"/extdata/",toytreefile,sep=""), max.read = 200)

#Calculate RERs
mamRERw = getAllResiduals(toyTrees, transform="sqrt", weighted=T, scale=T)

#Clades-based correlation results:
marineFg = c("Killer_whale", "Dolphin", "Walrus", "Seal", "Manatee")
sisters_marine = list("clade1"=c("Killer_whale", "Dolphin"))
useSpeciesList = toyTrees$masterTree$tip.label[20:62]
marineFgTree = foreground2TreeClades(marineFg,sisters_marine,toyTrees,plotTree=T)
fileteredMarineFgTree = foreground2TreeClades(marineFg,sisters_marine,toyTrees,plotTree=T, useSpecies = useSpeciesList)
pathvec = tree2PathsClades(marineFgTree, toyTrees)
res = correlateWithBinaryPhenotype(mamRERw, pathvec, min.sp=10, min.pos=2,
                                   weighted="auto")

#Non-Clades-Based correlation results:
marineNoCladesFgTree = foreground2Tree(marineFg,toyTrees,plotTree=T, clade = "all")
nonCladePath = tree2Paths(marineNoCladesFgTree, toyTrees)
noCladeRes = correlateWithBinaryPhenotype(mamRERw, nonCladePath, min.sp=10, min.pos=2,
                                          weighted="auto")
head(noCladeRes)
head(res)




library(RERconverge)
rerpath = find.package('RERconverge')

#read trees
toytreefile = "subsetMammalGeneTrees.txt" 
toyTrees=readTrees(paste(rerpath,"/extdata/",toytreefile,sep=""), max.read = 200)

#Calculate RERs
mamRERw = getAllResiduals(toyTrees, transform="sqrt", weighted=T, scale=T)

#Clades-based correlation results:
marineFg = c("Killer_whale", "Dolphin", "Walrus", "Seal", "Manatee")
sisters_marine = list("clade1"=c("Killer_whale", "Dolphin"))
marineFgTree = foreground2TreeClades(marineFg,sisters_marine,toyTrees,plotTree=F)
pathvec = tree2PathsClades(marineFgTree, toyTrees)
res = correlateWithBinaryPhenotype(mamRERw, pathvec, min.sp=10, min.pos=2,
                                   weighted="auto")

#Non-Clades-Based correlation results:
marineNoCladesFgTree = foreground2Tree(marineFg,toyTrees,plotTree=F)
nonCladePath = tree2Paths(marineNoCladesFgTree, toyTrees)
noCladeRes = correlateWithBinaryPhenotype(mamRERw, nonCladePath, min.sp=10, min.pos=2,
                                          weighted="auto")
head(noCladeRes)
head(res)


correl[order(correl$p.adj),]

# plot trees for sisterlist 

CVHCorrel = readRDS("Output/CVHRemake/CVHRemakeCorrelationFile.rds")
CVHCorrel[order(CVHCorrel$p.adj),]
dev.new()
readTest = readRDS("Output/CVHRemake/CVHRemakeBinaryForegroundTree.rds")
testTreeDisplayable = readTest
replace(testTreeDisplayable$edge.length, testTreeDisplayable$edge.length==0, 0.5)
replace(testTreeDisplayable$edge.length, testTreeDisplayable$edge.length==1, 4)

binaryTreePdfname = paste(outputFolderName, filePrefix, "BinaryForegroundTree.pdf", sep="")
pdf(binaryTreePdfname, width=8, height = 14)
plotTreeHighlightBranches(testTreeDisplayable, hlspecies=which(readTest$edge.length==1), hlcols="blue",)
dev.off()

binaryTree = readRDS(binaryTreeFilename)
plotTree(testTreeDisplayable)

?pruneTree()


domesticationTreeOld = readRDS("Output/Domestication/DomesticationBinaryForegroundTree.rds")
#write_tree(domesticationTreeOld, "Output/Domestication/DomesticationBinaryTreeNewick.txt")
plotTree(domesticationTreeOld)
domesticationTreeOldNew = read.tree("Output/Domestication/DomesticationBinaryTreeNewick.txt")


domesticationTreeOldNew2 = read.tree("Output/Domestication/DomesticationBinaryTreeNewick.txt")

readTest = domesticationTreeOldNew

testTreeDisplayable = readTest
testTreeDisplayable$edge.length = replace(testTreeDisplayable$edge.length, testTreeDisplayable$edge.length==0, 0.5)
testTreeDisplayable$edge.length = replace(testTreeDisplayable$edge.length, testTreeDisplayable$edge.length==1, 4)
testTreeDisplayableIntermediate = testTreeDisplayable
plotTree(testTreeDisplayable)

source("Src/Reu/ZoonomTreeNameToCommon.R")
readTest = domesticationTreeOldNew2
testTreeDisplayable = readTest
testTreeDisplayable$edge.length = replace(testTreeDisplayable$edge.length, testTreeDisplayable$edge.length==0, 0.5)
testTreeDisplayable$edge.length = replace(testTreeDisplayable$edge.length, testTreeDisplayable$edge.length==1, 4)

plotTreeHighlightBranches(testTreeDisplayable, hlspecies=which(readTest$edge.length==1), hlcols="blue",)
plotTree(testTreeDisplayable)
png("Output/Domestication/commonNameTree.png", width = 1980, height = 1040)
ZoonomTreeNameToCommon(domesticationTreeOldNew2)
dev.off()
source("Src/Reu/categorizePaths.R")
zoonomiaMasterObject = readRDS("Data/RemadeTreesAllZoonomiaSpecies.rds")
domesticationTreeFixed = categorizePaths(domesticationTreeOldNew2, zoonomiaMasterObject, "Domestication")
dev.off()
manualTree = readRDS("Results/DomesticationManualFGTree.rds")
ZoonomTreeNameToCommon(manualTree)


domesticationTree= readRDS("Output/Domestication/DomesticationBinaryForegroundTree.rds")
source("Src/Reu/ZoonomTreeNameToCommon.R")
ZoonomTreeNameToCommon(domesticationTree)





function (realcor, permvals) 
{
  permcor = permvals$corRho
  realstat = realcor$Rho
  names(realstat) = rownames(realcor)
  permcor = permcor[match(names(realstat), rownames(permcor)), 
  ]
  permpvals = vector(length = length(realstat))
  names(permpvals) = names(realstat)
  permstats = vector(length = length(realstat))
  names(permstats) = names(realstat)
  count = 1
  while (count <= length(realstat)) {
    if (is.na(realstat[count])) {
      permpvals[count] = NA
    }
    else {
      permcor_i = permcor[count, ]
      permcor_i = permcor_i[!is.na(permcor_i)]
      if (length(permcor_i) == 0) {
        permpvals[count] = NA
        permstats[count] = NA
      }
      else {
        median_permcor = median(permcor_i)
        if (realstat[count] >= median_permcor) {
          num = length(which(permcor_i >= realstat[count]))
          denom = length(which(permcor_i >= median_permcor))
        }
        else {
          num = length(which(permcor_i <= realstat[count]))
          denom = length(which(permcor_i <= median_permcor))
        }
        permpvals[count] = (num + 1)/(denom + 1)
        permstats[count] = -log10(permpvals[count]) * 
          sign(realstat[count] - median_permcor)
      }
    }
    count = count + 1
  }
  permstats[which(is.na(permpvals))] = NA
  out = data.frame(permpval = permpvals, permstats = permstats)
  out
}







