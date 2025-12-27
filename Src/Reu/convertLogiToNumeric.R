#Usage: Convert any vectors in a dataframe or list of dataframes into numeric vectors.
#Designed for use when dataframe contains vector of entirely NA, and thus is interpreted as logical when it should be numeric.
#Has both a Dataframe version and a List version. the List version is dependent on the dataframe version. 



convertLogiToNumericDataframe = function(inputDataframe){
  collumnsLogical = sapply(inputDataframe, is.logical)
  inputDataframe[,collumnsLogical] = lapply(inputDataframe[,collumnsLogical], as.numeric)
  return(inputDataframe)
}


convertLogiToNumericList = function(inputList){
  framesInList = sapply(inputList, is.data.frame)
  inputList[framesInList] = lapply(inputList[framesInList], convertLogiToNumericDataframe)
  return(inputList)
}