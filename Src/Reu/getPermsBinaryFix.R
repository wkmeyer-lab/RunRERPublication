function (numperms, fg_vec, sisters_list, root_sp, RERmat, trees, 
          mastertree, permmode = "cc", method = "k", min.pos = 2, trees_list = NULL, 
          calculateenrich = F, annotlist = NULL) 
{
  if(!is.numeric(numperms)){
    stop("numperms is not a numeric value")
  }
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