getPermPvalsCategorical2= function (realCors, nullPhens, phenvals, treesObj, RERmat, method = "kw", 
          min.sp = 10, min.pos = 2, winsorizeRER = NULL, winsorizetrait = NULL, 
          weighted = F, extantOnly = FALSE) 
{
  binary = FALSE
  if (method != "kw" & method != "aov") {
    message("Binary method provided. Setting binary to TRUE. Note: binary phenotype values should be TRUE and FALSE for correct results.")
    binary = TRUE
  }
  tree = treesObj$masterTree
  keep = intersect(names(phenvals), tree$tip.label)
  tree = pruneTree(tree, keep)
  if (is.rooted(tree)) {
    tree = unroot(tree)
  }
  if (!extantOnly) {
    if (!binary) {
      message("Generating null paths")
      nullPaths = lapply(nullPhens, function(x) {
        tr = tree
        tr$edge.length = c(x$tips, x$nodes)[tr$edge[, 
                                                    2]]
        tree2Paths(tr, treesObj, categorical = TRUE, 
                   useSpecies = names(phenvals))
        message("one")
      })
    }
    else {
      message("Generating null paths")
      nullPaths = lapply(nullPhens, function(x) {
        tr = tree
        tr$edge.length = c(x$tips, x$nodes)[tr$edge[, 
                                                    2]]
        tree2Paths(tr, treesObj, categorical = TRUE, 
                   useSpecies = names(phenvals)) - 1
      })
    }
  }
  message("Calculating correlation statistics")
  if (!extantOnly) {
    corsMatPvals = matrix(nrow = nrow(RERmat), ncol = length(nullPhens), 
                          dimnames = list(rownames(RERmat), NULL))
    corsMatEffSize = matrix(nrow = nrow(RERmat), ncol = length(nullPhens), 
                            dimnames = list(rownames(RERmat), NULL))
    message("Matrixes")
  }
  else {
    corsMatPvals = matrix(nrow = nrow(RERmat), ncol = nrow(nullPhens), 
                          dimnames = list(rownames(RERmat), NULL))
    corsMatEffSize = matrix(nrow = nrow(RERmat), ncol = nrow(nullPhens), 
                            dimnames = list(rownames(RERmat), NULL))
  }
  if (!binary) {
    Ppvals = lapply(1:length(realCors[[2]]), matrix, data = NA, 
                    nrow = nrow(RERmat), ncol = length(nullPhens), dimnames = list(rownames(RERmat), 
                                                                                   NULL))
    names(Ppvals) = names(realCors[[2]])
    Peffsize = lapply(1:length(realCors[[2]]), matrix, data = NA, 
                      nrow = nrow(RERmat), ncol = length(nullPhens), dimnames = list(rownames(RERmat), 
                                                                                     NULL))
    names(Peffsize) = names(realCors[[2]])
    message("pVals")
  }
  if (extantOnly) {
    if (!binary) {
      for (i in 1:nrow(nullPhens)) {
        cors = getAllCorExtantOnly(RERmat, nullPhens[i, 
        ], method = method, min.sp = min.sp, min.pos = min.pos, 
        winsorizeRER = winsorizeRER, winsorizetrait = winsorizetrait)
        corsMatPvals[, i] = cors[[1]]$P
        corsMatEffSize[, i] = cors[[1]]$Rho
        for (j in 1:length(cors[[2]])) {
          Ppvals[[names(cors[[2]])[j]]][, i] = cors[[2]][[j]]$P
          Peffsize[[names(cors[[2]])[j]]][, i] = cors[[2]][[j]]$Rho
        }
      }
    }
    else {
      for (i in 1:nrow(nullPhens)) {
        print(i)
        cors = getAllCorExtantOnly(RERmat, nullPhens[i, 
        ], method = method, min.sp = min.sp, min.pos = min.pos, 
        winsorizeRER = winsorizeRER, winsorizetrait = winsorizetrait)
        corsMatPvals[, i] = cors$P
        corsMatEffSize[, i] = cors$Rho
      }
    }
  }
  else {
    if (!binary) {
      for (i in 1:length(nullPaths)) {
        cors = getAllCor(RERmat, nullPaths[[i]], method = method, 
                         min.sp = min.sp, min.pos = min.pos, winsorizeRER = winsorizeRER, 
                         winsorizetrait = winsorizetrait, weighted = weighted)
        corsMatPvals[, i] = cors[[1]]$P
        corsMatEffSize[, i] = cors[[1]]$Rho
        for (j in 1:length(cors[[2]])) {
          Ppvals[[names(cors[[2]])[j]]][, i] = cors[[2]][[j]]$P
          Peffsize[[names(cors[[2]])[j]]][, i] = cors[[2]][[j]]$Rho
        }
      }
    }
    else {
      for (i in 1:length(nullPaths)) {
        cors = getAllCor(RERmat, nullPaths[[i]], method = method, 
                         min.sp = min.sp, min.pos = min.pos, winsorizeRER = winsorizeRER, 
                         winsorizetrait = winsorizetrait, weighted = weighted)
        corsMatPvals[, i] = cors$P
        corsMatEffSize[, i] = cors$Rho
      }
    }
  }
  message("Obtaining permulations p-values")
  if (!binary) {
    N = nrow(realCors[[1]])
    realCors[[1]]$permP = rep(NA, N)
    for (j in 1:length(realCors[[2]])) {
      realCors[[2]][[j]]$permP = rep(NA, N)
    }
    for (gene in 1:N) {
      if (is.na(realCors[[1]]$Rho[gene])) {
        p = NA
      }
      else {
        p = sum(corsMatEffSize[gene, ] > realCors[[1]]$Rho[gene], 
                na.rm = TRUE)/sum(!is.na(corsMatEffSize[gene, 
                ]))
      }
      realCors[[1]]$permP[gene] = p
      for (j in 1:length(realCors[[2]])) {
        if (is.na(realCors[[2]][[j]]$Rho[gene])) {
          p = NA
        }
        else {
          p = sum(abs(Peffsize[[names(realCors[[2]][j])]][gene, 
          ]) > abs(realCors[[2]][[j]]$Rho[gene]), na.rm = TRUE)/sum(!is.na(Peffsize[[names(realCors[[2]][j])]][gene, 
          ]))
        }
        realCors[[2]][[j]]$permP[gene] = p
      }
    }
  }
  else {
    N = nrow(realCors)
    realCors$permP = rep(NA, N)
    for (gene in 1:N) {
      if (is.na(realCors$Rho[gene])) {
        p = NA
      }
      else {
        p = sum(abs(corsMatEffSize[gene, ]) > abs(realCors$Rho[gene]), 
                na.rm = TRUE)/sum(!is.na(corsMatEffSize[gene, 
                ]))
      }
      realCors$permP[gene] = p
    }
  }
  message("Done")
  if (!binary) {
    return(list(res = realCors, pvals = list(corsMatPvals, 
                                             Ppvals), effsize = list(corsMatEffSize, Peffsize)))
  }
  else {
    return(list(res = realCors, pvals = corsMatPvals, effsize = corsMatEffSize))
  }
}

#--------------------------------------------------
#-----------------------------------------------------
#----------------------------------------------------

getPermPvalsCategorical3 = function (realCors, nullPhens, phenvals, treesObj, RERmat, method = "kw", 
          min.sp = 10, min.pos = 2, winsorizeRER = NULL, winsorizetrait = NULL, 
          weighted = F, extantOnly = FALSE) 
{
  binary = FALSE
  if (method != "kw" & method != "aov") {
    message("Binary method provided. Setting binary to TRUE. Note: binary phenotype values should be TRUE and FALSE for correct results.")
    binary = TRUE
  }
  tree = treesObj$masterTree
  keep = intersect(names(phenvals), tree$tip.label)
  tree = pruneTree(tree, keep)
  if (is.rooted(tree)) {
    tree = unroot(tree)
  }
  if (!extantOnly) {
    if (!binary) {
      message("Generating null paths")
      nullPaths = lapply(nullPhens, function(x) {
        tr = tree
        tr$edge.length = c(x$tips, x$nodes)[tr$edge[, 
                                                    2]]
        tree2Paths(tr, treesObj, categorical = TRUE, 
                   useSpecies = names(phenvals))
      })
      message("one")
    }
    else {
      message("Generating null paths")
      nullPaths = lapply(nullPhens, function(x) {
        tr = tree
        tr$edge.length = c(x$tips, x$nodes)[tr$edge[, 
                                                    2]]
        tree2Paths(tr, treesObj, categorical = TRUE, 
                   useSpecies = names(phenvals)) - 1
      })
    }
  }
  message("Calculating correlation statistics")
  if (!extantOnly) {
    corsMatPvals = matrix(nrow = nrow(RERmat), ncol = length(nullPhens), 
                          dimnames = list(rownames(RERmat), NULL))
    corsMatEffSize = matrix(nrow = nrow(RERmat), ncol = length(nullPhens), 
                            dimnames = list(rownames(RERmat), NULL))
    message("Matrixes")
  }
  else {
    corsMatPvals = matrix(nrow = nrow(RERmat), ncol = nrow(nullPhens), 
                          dimnames = list(rownames(RERmat), NULL))
    corsMatEffSize = matrix(nrow = nrow(RERmat), ncol = nrow(nullPhens), 
                            dimnames = list(rownames(RERmat), NULL))
  }
  if (!binary) {
    Ppvals = lapply(1:length(realCors[[2]]), matrix, data = NA, 
                    nrow = nrow(RERmat), ncol = length(nullPhens), dimnames = list(rownames(RERmat), 
                                                                                   NULL))
    names(Ppvals) = names(realCors[[2]])
    Peffsize = lapply(1:length(realCors[[2]]), matrix, data = NA, 
                      nrow = nrow(RERmat), ncol = length(nullPhens), dimnames = list(rownames(RERmat), 
                                                                                     NULL))
    names(Peffsize) = names(realCors[[2]])
  }
  if (extantOnly) {
    if (!binary) {
      for (i in 1:nrow(nullPhens)) {
        cors = getAllCorExtantOnly(RERmat, nullPhens[i, 
        ], method = method, min.sp = min.sp, min.pos = min.pos, 
        winsorizeRER = winsorizeRER, winsorizetrait = winsorizetrait)
        corsMatPvals[, i] = cors[[1]]$P
        corsMatEffSize[, i] = cors[[1]]$Rho
        message("One top")
        for (j in 1:length(cors[[2]])) {
          Ppvals[[names(cors[[2]])[j]]][, i] = cors[[2]][[j]]$P
          Peffsize[[names(cors[[2]])[j]]][, i] = cors[[2]][[j]]$Rho
          message("one bottom")
        }
      }
    }
    else {
      for (i in 1:nrow(nullPhens)) {
        print(i)
        cors = getAllCorExtantOnly(RERmat, nullPhens[i, 
        ], method = method, min.sp = min.sp, min.pos = min.pos, 
        winsorizeRER = winsorizeRER, winsorizetrait = winsorizetrait)
        corsMatPvals[, i] = cors$P
        corsMatEffSize[, i] = cors$Rho
      }
    }
  }
  else {
    if (!binary) {
      for (i in 1:length(nullPaths)) {
        cors = getAllCor(RERmat, nullPaths[[i]], method = method, 
                         min.sp = min.sp, min.pos = min.pos, winsorizeRER = winsorizeRER, 
                         winsorizetrait = winsorizetrait, weighted = weighted)
        corsMatPvals[, i] = cors[[1]]$P
        corsMatEffSize[, i] = cors[[1]]$Rho
        for (j in 1:length(cors[[2]])) {
          Ppvals[[names(cors[[2]])[j]]][, i] = cors[[2]][[j]]$P
          Peffsize[[names(cors[[2]])[j]]][, i] = cors[[2]][[j]]$Rho
        }
      }
    }
    else {
      for (i in 1:length(nullPaths)) {
        cors = getAllCor(RERmat, nullPaths[[i]], method = method, 
                         min.sp = min.sp, min.pos = min.pos, winsorizeRER = winsorizeRER, 
                         winsorizetrait = winsorizetrait, weighted = weighted)
        corsMatPvals[, i] = cors$P
        corsMatEffSize[, i] = cors$Rho
      }
    }
  }
  message("Obtaining permulations p-values")
  if (!binary) {
    N = nrow(realCors[[1]])
    realCors[[1]]$permP = rep(NA, N)
    for (j in 1:length(realCors[[2]])) {
      realCors[[2]][[j]]$permP = rep(NA, N)
    }
    for (gene in 1:N) {
      if (is.na(realCors[[1]]$Rho[gene])) {
        p = NA
      }
      else {
        p = sum(corsMatEffSize[gene, ] > realCors[[1]]$Rho[gene], 
                na.rm = TRUE)/sum(!is.na(corsMatEffSize[gene, 
                ]))
      }
      realCors[[1]]$permP[gene] = p
      for (j in 1:length(realCors[[2]])) {
        if (is.na(realCors[[2]][[j]]$Rho[gene])) {
          p = NA
        }
        else {
          p = sum(abs(Peffsize[[names(realCors[[2]][j])]][gene, 
          ]) > abs(realCors[[2]][[j]]$Rho[gene]), na.rm = TRUE)/sum(!is.na(Peffsize[[names(realCors[[2]][j])]][gene, 
          ]))
        }
        realCors[[2]][[j]]$permP[gene] = p
      }
    }
  }
  else {
    N = nrow(realCors)
    realCors$permP = rep(NA, N)
    for (gene in 1:N) {
      if (is.na(realCors$Rho[gene])) {
        p = NA
      }
      else {
        p = sum(abs(corsMatEffSize[gene, ]) > abs(realCors$Rho[gene]), 
                na.rm = TRUE)/sum(!is.na(corsMatEffSize[gene, 
                ]))
      }
      realCors$permP[gene] = p
    }
  }
  message("Done")
  if (!binary) {
    return(list(res = realCors, pvals = list(corsMatPvals, 
                                             Ppvals), effsize = list(corsMatEffSize, Peffsize)))
  }
  else {
    return(list(res = realCors, pvals = corsMatPvals, effsize = corsMatEffSize))
  }
}
