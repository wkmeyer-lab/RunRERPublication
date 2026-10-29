library(RERconverge)
source("Src/Reu/ZonomNameConvertMatrixCommon.R")
source("Src/Reu/ZonomNameConvertVector.R")
source("Src/Reu/categorizePaths.R")

#Read in the correlation file, RERs, and Paths produced by the original run RERs step 
allInsectivoryData = read.csv("Results/allInsectivoryCorrelationFile.csv")
RERs = readRDS("Results/allInsectivoryRERFile.rds")
#Note that this path is uncategorized. 
Paths = readRDS("Results/allInsectivoryPathsFile.rds")

#plot the uncategorized paths for 
plotRers(RERs,"KIAA0825", Paths )

#convert the zoonomia names to common names
CNRers = ZonomNameConvertMatrixCommon(RERs)

#Find indexes with "bat" in the species name 
CNNames = CNRers[1,]
names(CNNames)
bats = grep("bat", names(CNNames))

#Read in the binary tree (no paths)
binaryTree = readRDS("Results/allInsectivoryBinaryForegroundTree.rds")
plotTree(binaryTree)

#Read in the zoonomia master file. Ensure that this file is stored in the "Data" directory
#It is very large, and thus is not included in the Github! 
zonomMaster = readRDS("Data/RemadeTreesAllZoonomiaSpecies.rds")

#This line runs the manual selection 
#The "functionPaths" name can be changed, if a different file wants to be saved. 
#A tree with the manual annotations used is provided below, this code is present to show how it is generated.
#This function saves the *tree* to a file, and returns the path when called
functionPaths = categorizePaths(binaryTree, zonomMaster, "functionPath", overwrite = T) 

#This code loads in the premade tree, with the manual annotations I used. 
premadefucntionTree = readRDS("Results/premadefunctionPathManualFGTree.rds")
plotTree(premadefucntionTree)

#Convert the tree to a path; as the tree (not path) is what is saved 
premadefunctionPaths = tree2Paths(readRDS("Results/premadefunctionPathManualFGTree.rds"), zonomMaster)

#Plot the RERs 
plotRers(CNRers,"KIAA0825", Paths, sortrers = T)
plotRers(CNRers,"KIAA0825", functionPaths, sortrers = T)
plotRers(CNRers,"KIAA0825", premadefunctionPaths, sortrers = T)


plotRers(CNRers,"ZNF292", Paths, sortrers = T)
plotRers(CNRers,"ZNF292", functionPaths, sortrers = T)
plotRers(CNRers,"ZNF292", premadefunctionPaths, sortrers = T)


