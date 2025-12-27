library(RERconverge)
#Usage:
#Shows the RER plots of the genes with the top correlations, ticking through one every four seconds. 

showRERTreePlots = function(mainTrees, RERFile, CorrelationFile, usePerm = F, start = 1, style = "c"){
  if(usePerm){
    message('Expects the collumn containing the permulation P values to be named "permPVal"')
    genesRankedPermP = rownames(CorrelationFile[order(CorrelationFile$permPVal),])
  }else{
    genesRankedPermP = rownames(CorrelationFile[order(CorrelationFile$p.adj),])
  }
  
  for(i in start:100){
    treePlotRers(mainTrees, RERFile, genesRankedPermP[i], type = style)
    message(i)
    Sys.sleep(4)
  }
  
}