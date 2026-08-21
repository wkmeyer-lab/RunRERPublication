library(dplyr)

CategoricalPermulationGetCor =  function (realCors, nullPhens, phenvals, treesObj, RERmat, method = "kw", 
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
    tryCatch({
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
    }, error = function(i, corsMatPvals, corsMatEffSize, Ppvals, Peffsize){
      message("This error occurred on this Correlation:") 
      message(i)
      message("Skipping.")
      corsMatPvals[, i] = NA
      corsMatEffSize[, i] = NA
      for (j in 1:length(cors[[2]])) {
        Ppvals[[names(cors[[2]])[j]]][, i] = NA
        Peffsize[[names(cors[[2]])[j]]][, i] = NA
      }
    } )

  }
  output = list(corsMatEffSize, Peffsize, corsMatPvals, Ppvals)
  names(output) = c("corsMatEffSize", "Peffsize", "corsMatPvals", "Ppvals")
  return(output)
}


CategoricalCalculatePermulationPValues = function(realCors, intermediateList, start=1, end=NULL, report=F){
  {totalStart = Sys.time()}
  corsMatEffSize = intermediateList[[1]]
  Peffsize = intermediateList[[2]]
  corsMatPvals = intermediateList[[3]]
  Ppvals = intermediateList[[4]]
  message("Obtaining permulations p-values")
  N = nrow(realCors[[1]]) #
  if(start == 1){ #Only do this if start = 1, because otherwise it's already made and you'll overwrite the old script's results 
    realCors[[1]]$permP = rep(NA, N) #Make a column for permP values in all of the dataframes 
    for (j in 1:length(realCors[[2]])) {
      realCors[[2]][[j]]$permP = rep(NA, N) #Make a column for permP values in all of the dataframes 
    }
  }
  
  #Start updating the correlations
  if(is.null(end)){ #if no end specified
    stop = N
  }else{
    stop = end
  }
  for (gene in start:stop) {
    if(report){geneStart = Sys.time()}
    if (is.na(realCors[[1]]$Rho[gene])) {
      p = NA
    }
    else {
      signVal = sign(realCors[[1]]$Rho[gene])
      MatEffSizes = corsMatEffSize[gene, ]
      signedMatEffSizes = MatEffSizes[which(sign(MatEffSizes) == signVal)]
      p = sum(abs(signedMatEffSizes) > abs(realCors[[1]]$Rho[gene]), na.rm = TRUE)/(sum(!is.na(signedMatEffSizes))+1)
    }
    realCors[[1]]$permP[gene] = p
    for (j in 1:length(realCors[[2]])) {
      if (is.na(realCors[[2]][[j]]$Rho[gene])) {
        p = NA
      }
      else {
        realValue = realCors[[2]][[j]]$Rho[gene]
        signValue = sign(realValue)
        peffValues = Peffsize[[names(realCors[[2]][j])]][gene, ]
        signedPeffValues = peffValues[which( sign(peffValues) == signValue)]
        p = (sum(abs(signedPeffValues) > abs(realValue), na.rm = TRUE)+1)/ (sum(!is.na(signedPeffValues))+1)
      }
      realCors[[2]][[j]]$permP[gene] = p
    }
    if(report){geneEnd = Sys.time(); geneDuration = geneEnd - geneStart;message(paste("Completed Gene", gene, "Duration", geneDuration, attr(geneDuration, "units")))}
  }
  message("Done")
  {totalEnd = Sys.time(); totalDuration = totalEnd - totalStart;message(paste("Completed p-Values; Duration", totalDuration, attr(totalDuration, "units")))}
  return(list(res = realCors, pvals = list(corsMatPvals, Ppvals), effsize = list(corsMatEffSize, Peffsize)))
  
}

CategoricalCollectIntermediateResults = function(realCors, intermediateList, initial = F, start=1, end=NULL, report=F){
  {totalStart = Sys.time()}
  corsMatEffSize = intermediateList[[1]]
  Peffsize = intermediateList[[2]]
  corsMatPvals = intermediateList[[3]]
  Ppvals = intermediateList[[4]]
  message("Collecting permulation intermediates")
  N = nrow(realCors[[1]]) #
  if(initial){ #Only do this if start = 1, because otherwise it's already made and you'll overwrite the old script's results 
    realCors[[1]]$numMoreExtremePerms = rep(0, N) #Make a column for permP values in all of the dataframes 
    realCors[[1]]$numTotalPerms = rep(0, N) #Make a column for permP values in all of the dataframes 
    realCors[[1]]$permP = rep(NA, N) #Make a column for permP values in all of the dataframes 
    for (j in 1:length(realCors[[2]])) {
      realCors[[2]][[j]]$numMoreExtremePerms = rep(0, N) #Make a column for permP values in all of the dataframes 
      realCors[[2]][[j]]$numTotalPerms = rep(0, N) #Make a column for permP values in all of the dataframes 
      realCors[[2]][[j]]$permP = rep(NA, N) #Make a column for permP values in all of the dataframes 
    }
  }
  
  #Start updating the correlations
  if(is.null(end)){ #if no end specified
    stop = N
  }else{
    stop = end
  }
  for (gene in start:stop) {
    if(report){geneStart = Sys.time()}
    if (is.na(realCors[[1]]$Rho[gene])) {
      p = NA
    }
    else {
      signVal = sign(realCors[[1]]$Rho[gene])
      MatEffSizes = corsMatEffSize[gene, ]
      signedMatEffSizes = MatEffSizes[which(sign(MatEffSizes) == signVal)]
      newMoreExtreme = sum(abs(signedMatEffSizes) > abs(realCors[[1]]$Rho[gene]), na.rm = TRUE)
      newTotal = (sum(!is.na(signedMatEffSizes)))
      realCors[[1]]$numMoreExtremePerms[gene] = realCors[[1]]$numMoreExtremePerms[gene] + newMoreExtreme 
      realCors[[1]]$numTotalPerms[gene] =  realCors[[1]]$numTotalPerms[gene] + newTotal
    }
    for (j in 1:length(realCors[[2]])) {
      if (is.na(realCors[[2]][[j]]$Rho[gene])) {
        p = NA
      }
      else {
        realValue = realCors[[2]][[j]]$Rho[gene]
        signValue = sign(realValue)
        peffValues = Peffsize[[names(realCors[[2]][j])]][gene, ]
        signedPeffValues = peffValues[which( sign(peffValues) == signValue)]
        newMoreExtreme = (sum(abs(signedPeffValues) > abs(realValue), na.rm = TRUE))
        newTotal = (sum(!is.na(signedPeffValues))+1)
        realCors[[2]][[j]]$numMoreExtremePerms[gene] = realCors[[2]][[j]]$numMoreExtremePerms[gene] + newMoreExtreme
        realCors[[2]][[j]]$numTotalPerms[gene] = realCors[[2]][[j]]$numTotalPerms[gene] + newTotal
      }
      
    }
    if(report){geneEnd = Sys.time(); geneDuration = geneEnd - geneStart;message(paste("Completed Gene", gene, "Duration", geneDuration, attr(geneDuration, "units")))}
  }
  message("Done")
  {totalEnd = Sys.time(); totalDuration = totalEnd - totalStart;message(paste("Completed p-Values; Duration", totalDuration, attr(totalDuration, "units")))}
  return(realCors)
  
}





CategoricalCalculatePValueFromCollectedIntermediates = function(CollectedIntermediates, start=1, end=NULL, report=F){
  realCors = CollectedIntermediates
  N = nrow(realCors[[1]]) #
  totalStart = Sys.time()
  
  realCors[[1]]$permP.adj = rep(NA, N) #Make a column for permP values in all of the dataframes 
  realCors[[1]]$permRank = rep(0, N) #Make a column for permP values in all of the dataframes 
  realCors[[1]]$unpermRank = rep(0, N) #Make a column for permP values in all of the dataframes 
  realCors[[1]]$rhoRank = rep(0, N) #Make a column for permP values in all of the dataframes 
  realCors[[1]]$absRhoRank = rep(0, N) #Make a column for permP values in all of the dataframes 
  realCors[[1]]$index = seq_len(nrow(realCors[[1]])) #Make a column for permP values in all of the dataframes 
  for (j in 1:length(realCors[[2]])) {
    realCors[[2]][[j]]$permP.adj = rep(NA, N) #Make a column for permP values in all of the dataframes 
    realCors[[2]][[j]]$permRank = rep(0, N) #Make a column for permP values in all of the dataframes 
    realCors[[2]][[j]]$unpermRank = rep(0, N) #Make a column for permP values in all of the dataframes 
    realCors[[2]][[j]]$rhoRank = rep(0, N) #Make a column for permP values in all of the dataframes 
    realCors[[2]][[j]]$absRhoRank = rep(0, N) #Make a column for permP values in all of the dataframes 
    realCors[[2]][[j]]$index = seq_len(nrow(realCors[[2]][[j]])) #Make a column for permP values in all of the dataframes 
  }
  
  
  #Start updating the correlations
  if(is.null(end)){ #if no end specified
    stop = N
  }else{
    stop = end
  }
  for (gene in start:stop) {
    if(report){geneStart = Sys.time()}
    if (is.na(realCors[[1]]$Rho[gene])) {
      p = NA
    }
    else {
      p = realCors[[1]]$numMoreExtremePerms[gene]/(realCors[[1]]$numTotalPerms[gene]+1)
    }
    realCors[[1]]$permP[gene] = p
    
    for (j in 1:length(realCors[[2]])) {
      if (is.na(realCors[[2]][[j]]$Rho[gene])) {
        p = NA
      }
      else {
        p = (realCors[[2]][[j]]$numMoreExtremePerms[gene])/ (realCors[[2]][[j]]$numTotalPerms[gene]+1)
      }
      realCors[[2]][[j]]$permP[gene] = p
    }
    if(report){geneEnd = Sys.time(); geneDuration = geneEnd - geneStart;message(paste("Completed Gene", gene, "Duration", geneDuration, attr(geneDuration, "units")))}
  }
  message("Done")
  message("Perforing MHTC and ranking")
  
  realCors[[1]] = processPermulatedPValue(realCors[[1]])
  for (j in 1:length(realCors[[2]])) {
    realCors[[2]][[j]] = processPermulatedPValue(realCors[[2]][[j]])
  }
  
  
  
  {totalEnd = Sys.time(); totalDuration = totalEnd - totalStart;message(paste("Completed p-Values; Duration", totalDuration, attr(totalDuration, "units")))}
  CollectedIntermediates = realCors
  
  return(CollectedIntermediates)
}


processPermulatedPValue = function(inputData, adjustMetod = "BH"){
  inputData$permP.adj = p.adjust(inputData$permP, adjustMetod) 
  
  inputData = inputData[order(inputData$permP.adj),]
  inputData$permRank = seq_len(nrow(inputData))
  
  inputData = inputData[order(inputData$p.adj),]
  inputData$unpermRank = seq_len(nrow(inputData))
  
  inputData = inputData[order(inputData$Rho),]
  inputData$rhoRank = seq_len(nrow(inputData))
  
  inputData = inputData[order(abs(inputData$Rho)),]
  inputData$absRhoRank = seq_len(nrow(inputData))
  
  inputData = inputData[order(abs(inputData$index)),]
  
  inputData = inputData %>% relocate(numMoreExtremePerms, .after = last_col())
  inputData = inputData %>% relocate(numTotalPerms, .after = last_col())
  
  
  inputData
}

#