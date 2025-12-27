library(RERconverge)
#Usage:
#Shows the RER plots of the genes with the top correlations, ticking through one every four seconds. 

showRERPlots = function(RERFile, PathsFile, CorrelationFile, usePerm = F, start = 1, sort = F, end = 15, sleep = 4){
  if(usePerm){
    message('Expects the collumn containing the permulation P values to be named "permPVal"')
    genesRankedPermP = rownames(CorrelationFile[order(CorrelationFile$permPVal),])
  }else{
    genesRankedPermP = rownames(CorrelationFile[order(CorrelationFile$p.adj),])
  }
  
  for(i in start:end){
    plotRers(RERFile, genesRankedPermP[i], PathsFile, sortrers = sort)
    message(i)
    Sys.sleep(sleep)
  }
  
}

