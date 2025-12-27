#Usage: 
#Provide the function with he binary tree of interest, the pre-loaded object of the full zoonomia dataset, and a prefix to use when saving the tree. 
categorizePaths = function(binaryPhenotypeTree, zonomMasterObject, prefix, overwrite = F){
  filename = paste("Results/", prefix, "ManualFGTree.rds", sep='') #use prefix argument to generate filename
  message(filename)                                                #Send the filename 
  if(file.exists(filename) && overwrite == F){                     #If the file of that filename exists, and not told to overwrite
    manualFGTree = readRDS(filename)                               #Use that file, go to return step 
  }else{                                                           #else
    source("Src/Reu/ZonomNameConvertVector.R")                     
    binaryTree = binaryPhenotypeTree                              #use binary tree provided 
    renamedTree = binaryTree                                      #make a copy of the tree to rename the tips on 
    renamedTree$tip.label = ZonomNameConvertVectorCommon(binaryTree$tip.label) #rename tips for readability
    manualFGTree = click_select_foreground_branches(renamedTree)  #use click_select to have the user select a manual path
    manualFGTree$tip.label = binaryTree$tip.label                 #reset the tiplabel to the zoonomia name
    saveRDS(manualFGTree, filename)                               #save the tree to the filename
  }
  manualPaths = tree2Paths(manualFGTree, zonomMasterObject)       #Convert the tree to a path using zoonomia dataset
  return(manualPaths)
}