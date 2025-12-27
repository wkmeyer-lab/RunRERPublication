#This function will combine the list of the intermediates produced by CategoricalPermulationGetCor; to be fed into CategoricalCalculatePermulationPValues

combineCategoricalPermulationIntermediates = function(intermediate1, intermediate2){
  intermediate1[[1]] = cbind(intermediate1[[1]], intermediate2[[1]])
  for(i in 1:length(intermediate1[[2]])){
    intermediate1[[2]][[i]] = cbind(intermediate1[[2]][[i]], intermediate2[[2]][[i]])
  }
  intermediate1[[3]] = cbind(intermediate1[[3]], intermediate2[[3]])
  for(i in 1:length(intermediate1[[4]])){
    intermediate1[[4]][[i]] = cbind(intermediate1[[4]][[i]], intermediate2[[4]][[i]])
  }
  return(intermediate1)
}