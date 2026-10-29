#This is a home for a variety of once-off or manually operated code
#that interacts with the other scripts in this repo. 
#Each titled section is written with a different goal, and many 
#are not actually able to be run on their own, as they interact 
#or depend on other scripts in the repo, sometimes partially-run scripts 
#(for example, code used to debug, or tweak a script while it is running).

#Many of these blocks will not function on their own, expecting either 
#partially run other scripts, other scratch blocks, or even older now-removed 
#versions of other parts of their own block. As this is a 'scratch space' for
#one-off code, it is not nearly as polished as anything else. 
#While some of the blocks can be run as a self-contained mini-script, most are
#better thought of as effectively a premenant record of lines of code that
#were manually run through the console. That was the original purpose of these 
#documents -- instead of running a line of code in the ephemeral console, it was
#written in one of these documents, and then executed. 

#Because many of these scripts are therefore disjointed, unfinished, or 
#otherwise break good practices because they were meant to be run with 
#constant manual checking, running them unsupervised could cause Many 
#problems, notably including OVERWRITTING STANDARD OUTPUT FILES (as this
#is an anything-goes record of code run during development). While these
#most dangerous lines are often commented out, they aren't always. 

#As such, there is a line that prevents automatic execution of this file at 
#start of the code, because this was never meant to be run that way, but it 
#is easy to hit a "run-all" shortcut. 


a = b #prevent full runs, as each of the bracketed sections is mea
library(RERconverge)
library(tools)
library(scales)
library(data.table)
library(geiger)

source("Src/Reu/ZoonomTreeNameToCommon.R")
source("Src/Reu/RERCOnvergeFunctions.R")

par(mfrow=c(1,1))
palette(c("#1B9E77","#7570B3", "#000000", "#E7298A"))



palette(c("#1B9E77", "#000000", "#7570B3", "#E7298A"))


packageVersion("RERconverge")




#---------------------------------------------------------------------
# --- Longevity results  --- 
# --------------------------------------------------------------------

source("Src/Reu/CategoricalDroptip.R")
source("Src/Reu/ZoonomTreeNameToCommon.R")
source("Src/Reu/convertColNameToIndex.R")
source('Src/Reu/RERCOnvergefunctions.R')
palette(c( "darkgreen", "darkblue","black", "red", "yellow"))
nameColumn = "ZoonomiaTip"

mergedData = read.csv("data/mergedData.csv")
mainAnalysisSpecies = readRDS("Output/ComplexDietCentralAnalysisSimplify2/ComplexDietCentralAnalysisSimplify2SpeciesFilter.rds")
mainAllSpecies = readRDS("Output/ComplexDietCentralAnalysisSimplifyAllSpecies/ComplexDietCentralAnalysisSimplifyAllSpeciesSpeciesFilter.rds")
ageDataSheet = read.csv("Data/MaturityLifespanData.csv")
mainTrees = readRDS("Data/zoonomiaAllMammalsTrees.rds")

mainCategoricalTree = readRDS("Output/ComplexDietCentralAnalysisSimplify2/ComplexDietCentralAnalysisSimplify2CategoricalTree.rds")
plotTreeCategorical(mainCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = masterTree)

TEDatasheet = read.csv("Data/TEContentOsmanski2023TS8.csv")
TEDatasheet$Zoonomia = mergedData$ZoonomiaTip[match(tolower(TEDatasheet$Species), mergedData$Scientific_Binomial)]

cancerDatasheet = read.csv("Data/species-cancer-prevalence-data.csv")
cancerDatasheet$Zoonomia = mergedData$ZoonomiaTip[match(gsub(" ", "_", tolower(cancerDatasheet$Species)), mergedData$Scientific_Binomial)]




masterTree = mainTrees$masterTree
categoricalTree = mainCategoricalTree
TreeTipCol = "Zoonomia"
DataColumn = "MaximumAge"
DataColumn = "combinedBodysize"




makeSubsetCategoricalTreeToContiousData = function(DataColumn = get("DataColumn", envir = .GlobalEnv), metric = "diff",
                                                   categoricalTree = get("categoricalTree", envir = .GlobalEnv), masterTree = get("masterTree", envir = .GlobalEnv), 
                                                   TreeTipCol = get("TreeTipCol", envir = .GlobalEnv), continousDataSheet = get("continousDataSheet", envir = .GlobalEnv)  
                                                   ){
  categoricalSpecies = categoricalTree$tip.label
  continousSpecies = continousDataSheet[which(!is.na(continousDataSheet[,convertColnameToIndex(DataColumn, continousDataSheet)])), convertColnameToIndex(TreeTipCol, continousDataSheet)]
  continousSpecies = continousSpecies[-which(is.na(continousSpecies))]
  overlapSpecies = continousSpecies[which(continousSpecies %in% categoricalSpecies)]

  message(paste("Overlap size:", length(overlapSpecies)))
  
  tipsToRemove = categoricalTree$tip.label[!categoricalTree$tip.label %in% overlapSpecies]
  overlapTree = categoricalDropTip(categoricalTree, tipsToRemove)
  
  masterTreeTipsRemove = masterTree$tip.label[!masterTree$tip.label %in% overlapSpecies]
  overlapMasterTree = drop.tip(masterTree, masterTreeTipsRemove)
  
  
  
  
  #re-order the overlap (phenotype) tree to the same edge order as the master tree 
  overlapTipToMasterTip = match(overlapTree$tip.label, overlapMasterTree$tip.label)
  names(overlapTipToMasterTip) = seq_along(overlapTipToMasterTip)
  overlapTipToMasterTip = overlapTipToMasterTip[order(overlapTipToMasterTip)]
  
  reconstructedOverlapTree = overlapTree
  reconstructedOverlapTree$tip.label = reconstructedOverlapTree$tip.label[as.integer(names(overlapTipToMasterTip))]
  
  
  nodeMatch = matchAllNodes(overlapTree, overlapMasterTree)
  
  #We don't need to rearrange the edge length values, because we already changed the tip label positions. 
  edgeReconstruction = data.frame(seq_along(overlapTree$edge.length))
  names(edgeReconstruction) = "EdgeIndex"
  edgeReconstruction$originalStart = overlapTree$edge[,1]
  edgeReconstruction$originalEnd = overlapTree$edge[,2]
  edgeReconstruction$newStart = nodeMatch[match(overlapTree$edge[,1], nodeMatch[, 1]), 2]
  edgeReconstruction$newEnd = nodeMatch[match(overlapTree$edge[,2], nodeMatch[, 1]), 2]
  
  
  
  reconstructedOverlapTree$edge[,1] = edgeReconstruction$newStart
  reconstructedOverlapTree$edge[,2] = edgeReconstruction$newEnd

  commonOverlapTree = ZoonomTreeNameToCommon(overlapTree, tipCol = nameColumn, plot = T)
  commonOverlapMasterTree = ZoonomTreeNameToCommon(overlapMasterTree, tipCol = nameColumn, plot = T)
  commonReconstructedTree =  ZoonomTreeNameToCommon(reconstructedOverlapTree, tipCol = nameColumn, plot = T)
  

  # Plot the trees to confirm the conversion was successful
  par(mfrow = c(1, 2))
  
  tryCatch({
    plotTreeCategorical(commonOverlapTree,c("Herbivore", "Insectivore", "Omnivore", "Vertivore"),master = commonOverlapMasterTree)
    plotTreeCategorical(commonReconstructedTree,c("Herbivore", "Insectivore", "Omnivore", "Vertivore"),master = commonOverlapMasterTree)
  }, error = function(e) {
    message("Non-Critical Error plotting trees: ", conditionMessage(e))
  }, finally = {
    par(mfrow = c(1, 1))
  })
  
  
  
  #Make the tree with continuous information
  overlapData = continousDataSheet[continousDataSheet[,convertColnameToIndex(TreeTipCol, continousDataSheet)] %in% overlapSpecies, c(convertColnameToIndex(TreeTipCol, continousDataSheet), convertColnameToIndex(DataColumn, continousDataSheet))]
  overlapData = overlapData[order(match(overlapData[,1], overlapMasterTree$tip.label)),]
  overlapDataVec = overlapData[,2]
  names(overlapDataVec) = overlapData[,1]
  
  continousTree = edgeVars(overlapMasterTree, overlapDataVec, metric = metric)
  
  
  
  if(!all.equal(overlapMasterTree$tip.label, reconstructedOverlapTree$tip.label)&
  all.equal(overlapMasterTree$edge, continousTree$edge)&
  all.equal(overlapMasterTree$edge, continousTree$edge)&
  all.equal(overlapMasterTree$edge, reconstructedOverlapTree$edge)){
    stop("Error, trees do not match")
  }
  
  
  edgeTable = data.frame(overlapMasterTree$edge.length, reconstructedOverlapTree$edge.length, continousTree$edge.length)
  names(edgeTable) = c("Master Tree length", "Phenotype", "Continuous Value")
  
  TreeOutputs = list(edgeTable, overlapMasterTree, reconstructedOverlapTree, continousTree)
  names(TreeOutputs) = c("edgeTable", "overlapMasterTree", "phenotypeTree", "continousTree")
  TreeOutputs
}



AgeDiffTrees = makeSubsetCategoricalTreeToContiousData("MaximumAge", "diff", continousDataSheet = ageDataSheet)
AgeMeanTrees = makeSubsetCategoricalTreeToContiousData("MaximumAge", "mean", continousDataSheet = ageDataSheet)

TeDiffTrees = makeSubsetCategoricalTreeToContiousData("Young_DNA_proportion", "diff", continousDataSheet = TEDatasheet)
TeMeanTrees = makeSubsetCategoricalTreeToContiousData("Young_DNA_proportion", "mean", continousDataSheet = TEDatasheet)

CancerDiffTrees = makeSubsetCategoricalTreeToContiousData("NeoplasiaPrevalence", "diff", continousDataSheet = cancerDatasheet)
CancerMeanTrees = makeSubsetCategoricalTreeToContiousData("NeoplasiaPrevalence", "mean", continousDataSheet = cancerDatasheet)



comboOutput = list(AgeDiffTrees, AgeMeanTrees, TeDiffTrees, TeMeanTrees, CancerDiffTrees, CancerMeanTrees)
names(comboOutput) = c("AgeDiffTrees", "AgeMeanTrees", "TeDiffTrees", "TeMeanTrees", "CancerDiffTrees", "CancerMeanTrees")

saveRDS(comboOutput, "Results/DnaRepairResultsForWynn.rds")


diffBodysizeTrees = makeSubsetCategoricalTreeToContiousData("combinedBodysize", "diff")
meanBodysizeTrees = makeSubsetCategoricalTreeToContiousData("combinedBodysize", "mean")



#
























#---------------------------------------------------------------------
# --- Misc paper checks--- 
# --------------------------------------------------------------------

mergedData = read.csv("Data/mergedData.csv")

mergedData$CommonName[which(mergedData$SimplifiedDietConvertEqualDrop =="zMixedPredator")]

supData1 = read.csv("Results/SupplementaryTable1 -- SpeciesNamesTipsAndDiet.csv")



supData1$HerbivoreType = rep("TBD")

supData1$HerbivoreType[!supData1$ThreeDietAnalysis %in% "Herbivore"] = NA



#supData1$HerbivoreType[which( (supData1$Diet.Fruit + supData1$Diet.Nect + supData1$Diet.Seed) >= 90)] = "Potential Specialist Herbivore"

supData1$HerbivoreType[which(supData1$Diet.Seed > 50)] = "Granivore"
supData1$HerbivoreType[which(supData1$Diet.Nect > 50)] = "Nectivore"
supData1$HerbivoreType[which(supData1$Diet.Fruit > 50)] = "Frugivore"
supData1$HerbivoreType[which(supData1$Diet.PlantO > 50)] = "Loose Generalist Herbivore"

supData1$HerbivoreType[which(supData1$Diet.PlantO >= 90)] = "Generalist Herbivore"
supData1$HerbivoreType[which(supData1$Diet.Seed >= 90)] = "Strict Granivore"
supData1$HerbivoreType[which(supData1$Diet.Nect >= 90)] = "Strict Nectivore"
supData1$HerbivoreType[which(supData1$Diet.Fruit >= 90)] = "Strict Frugivore"

supData1$HerbivoreType[!supData1$ThreeDietAnalysis %in% "Herbivore"] = NA


supData1$HerbivoreType[supData1$TreeTip == "vs_NA"]


supData1$HerbivoreType[which(supData1$HerbivoreType == "TBD")] = "Mixed Herbivore"

table(supData1$HerbivoreType)


supData1$HerbivoreTypeInMain = supData1$HerbivoreType

supData1$HerbivoreTypeInMain[!supData1$InMainAnalysis] = NA

table(supData1$HerbivoreTypeInMain)





supData1 = supData1[,c(2:11, 24, 12:23)]

#write.csv(supData1, "Results/SupplementaryTable1 -- SpeciesNamesTipsAndDiet.csv")


supData1 = mergedData

supData1$HerbivoreType = rep("TBD")

supData1$HerbivoreType[!supData1$SimplifiedDietConvert %in% "Herbivore"] = NA

supData1$HerbivoreType[which(supData1$Diet.PlantO >= 90)] = "Generalist Herbivore"

#supData1$HerbivoreType[which( (supData1$Diet.Fruit + supData1$Diet.Nect + supData1$Diet.Seed) >= 90)] = "Potential Specialist Herbivore"

supData1$HerbivoreType[which(supData1$Diet.Seed > 50)] = "Granivore"
supData1$HerbivoreType[which(supData1$Diet.Nect > 50)] = "Nectivore"
supData1$HerbivoreType[which(supData1$Diet.Fruit > 50)] = "Frugivore"

supData1$HerbivoreType[which(supData1$Diet.Seed >= 90)] = "Strict Granivore"
supData1$HerbivoreType[which(supData1$Diet.Nect >= 90)] = "Strict Nectivore"
supData1$HerbivoreType[which(supData1$Diet.Fruit >= 90)] = "Strict Frugivore"

supData1$HerbivoreType[!supData1$SimplifiedDietConvert %in% "Herbivore"] = NA




supData1$HerbivoreType[which(supData1$HerbivoreType == "TBD")] = "Mixed Herbivore"

mergedData = supData1

mergedData$ComplexHerbDiet = mergedData$SimplifiedDietConvert

mergedData$ComplexHerbDiet[which(mergedData$ComplexHerbDiet == "Herbivore")] = mergedData$HerbivoreType[which(mergedData$ComplexHerbDiet == "Herbivore")]

table(mergedData$ComplexHerbDiet)

write.csv(mergedData, "Data/mergedData.csv")


palette(c("hotpink", "darkgreen", "orange", "blue", "yellow", "black", "red"))
plotTreeCategorical(commonCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore", "zMixedPredator"), master = commonMainTrees$masterTree)
plotTreeCategorical(categoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore", "zMixedPredator"), master = mainTrees$masterTree)



#---------------------------------------------------------------------
# --- Make new alt reconstrcution --- 
# --------------------------------------------------------------------

mainTreesBackup = mainTrees
commonMainTreesBackup = commonMainTrees

mainTrees$masterTree$edge.length = rep(1, length(mainTrees$masterTree$edge.length))
commonMainTrees$masterTree$edge.length = rep(1, length(commonMainTrees$masterTree$edge.length))

plotTreeCategorical(commonCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = commonMainTreesBackup$masterTree)
plotTreeCategorical(categoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = mainTreesBackup$masterTree)


equalCategoricalTree = readRDS("output/ComplexDietCentralAnalysisSimplifyAltReconstruction/ComplexDietCentralAnalysisSimplifyAltReconstructionCategoricalTree.rds")
mainCategoricalTree = readRDS("output/ComplexDietCentralAnalysisSimplify2/ComplexDietCentralAnalysisSimplify2CategoricalTree.rds")

equalCategoricalTree$edge.length == mainCategoricalTree$edge.length

mismatchBranches = which(!equalCategoricalTree$edge.length == mainCategoricalTree$edge.length)
length(mismatchBranches)

edgelabels(edge = mismatchBranches, cex=0.4, col= "purple", frame="none")

clades$phylopic


mainResults = readRDS("Output/ComplexDietCentralAnalysisSimplify2/ComplexDietCentralAnalysisSimplify2combinedGeneResultsWithAlternates.rds")

length(which(mainResults$`IV-Rho`[which(mainResults$`IV-significantRobust`)] <0 ))


#---------------------------------------------------------------------
# --- New PhenotypeBInary Analyses --- 
# --------------------------------------------------------------------

#Load central files 
centralRER = readRDS("Output/ComplexDietBinaryCentral/ComplexDietBinaryCentralRERFile.rds")
centralSpeciesFilter = readRDS("Output/ComplexDietBinaryCentral/ComplexDietBinaryCentralSpeciesFilter.rds")
centralCategoricalTree = readRDS("Output/ComplexDietBinaryCentral/ComplexDietBinaryCentralCategoricalTree.rds")
centralPhenotypeVector = readRDS("Output/ComplexDietBinaryCentral/ComplexDietBinaryCentralCategoricalPhenotypeVector.rds")

mainTrees = readRDS("data/zoonomiaAllMammalsTrees.rds")


#Make binary specific outputs

infileprefix = "ComplexDietBinaryCarnivore"
phenName = "Invertivore"
newPhen = NULL
plot = T

phenName = c("Invertivore", "Vertivore")
newPhen = "Carnivore"

makeBinaryPhen = function(infileprefix, phenName, newPhen = NULL, plot = T){

  filePrefix = infileprefix
  if(is.null(newPhen)){newPhen = phenName}
  
  #Make output directory
  #  Output Directory 
  if(!dir.exists("Output")){                                      #Make output directory if it does not exist
    dir.create("Output")
  }
  outputFolderNameNoSlash = paste("Output/",filePrefix, sep = "") #Set the prefix sub directory
  if(!dir.exists(outputFolderNameNoSlash)){                       #create that directory if it does not exist
    dir.create(outputFolderNameNoSlash)
  }
  outputFolderName = paste("Output/",filePrefix,"/", sep = "")
  
  newPhenotypeVector = centralPhenotypeVector
  newPhenotypeVector[which(!newPhenotypeVector %in% phenName)] = "Background"
  newPhenotypeVector[which(newPhenotypeVector %in% phenName)] = newPhen
  
  
  #get the selected phen index
  phenlist = names(table((centralPhenotypeVector)))
  
  phenIndex = which(phenlist %in% phenName)
  
  newCategoricalTree = centralCategoricalTree
  
  #temporarily store the background 
  newCategoricalTree$edge.length[which(!newCategoricalTree$edge.length %in% phenIndex)] = 1000
  
  #temporarily store the foreground
  newCategoricalTree$edge.length[which(newCategoricalTree$edge.length %in% phenIndex)] = 2000
  
  newCategoricalTree$edge.length[which(newCategoricalTree$edge.length == 1000)] = 1
  newCategoricalTree$edge.length[which(newCategoricalTree$edge.length == 2000)] = 2

  
  if(plot){
    plotTreeCategorical(newCategoricalTree, c("Background", newPhen), master = mainTrees$masterTree)
    
    treeImageFilename = paste(outputFolderName, filePrefix, "CategoricalTree.pdf", sep="")
    pdf(treeImageFilename, height = length(newPhenotypeVector)/18, width = 10)  
    plotTreeCategorical(newCategoricalTree, c("Background", newPhen), master = mainTrees$masterTree)
    dev.off()
    
  }
  
  #Save the files
  speciesFilterFilename = paste(outputFolderName, filePrefix, "SpeciesFilter.rds",sep="")
  saveRDS(centralSpeciesFilter, speciesFilterFilename)
  
  RERFileName = paste(outputFolderName, filePrefix, "RERFile.rds", sep= "") 
  saveRDS(centralRER, RERFileName)
  
  phenotypeVectorFilename = paste(outputFolderName, filePrefix, "CategoricalPhenotypeVector.rds",sep="") #make a filename based on the prefix
  saveRDS(newPhenotypeVector, file = phenotypeVectorFilename)  
  
  categoricalTreeFilename = paste(outputFolderName, filePrefix, "CategoricalTree.rds", sep="") #make a filename based on the prefix
  saveRDS(newCategoricalTree, categoricalTreeFilename) 
}

  
makeBinaryPhen("ComplexDietBinaryInvertivore", "Invertivore")
makeBinaryPhen("ComplexDietBinaryVertivore", "Vertivore")
makeBinaryPhen("ComplexDietBinaryHerbivore", "Herbivore")  
makeBinaryPhen("ComplexDietBinaryOmnivore", "Omnivore")
makeBinaryPhen("ComplexDietBinaryCarnivore", c("Invertivore", "Vertivore"), "Carnivore")

# ------ 

testCarVec = readRDS(phenotypeVectorFilename)

filenamescafold = function(prefix){
  name = paste0("Output/ComplexDietBinary", prefix, "/Background-", prefix,"/ComplexDietBinary", prefix, "Background-", prefix, "CorrelationFile.rds")
  name
}

carnCor = readRDS(filenamescafold("Carnivore"))
invCor = readRDS(filenamescafold("Invertivore"))
vertCor = readRDS(filenamescafold("Vertivore"))
herbCor = readRDS(filenamescafold("Herbivore"))
omniCor = readRDS(filenamescafold("Omnivore"))

names(carnCor) = paste0("C-", names(carnCor))
names(invCor) = paste0("I-", names(invCor))
names(vertCor) = paste0("V-", names(vertCor))
names(herbCor) = paste0("H-", names(herbCor))
names(omniCor) = paste0("O-", names(omniCor))

binaryCors = cbind(carnCor, invCor, vertCor, herbCor, omniCor)

combinedBinaryCorsFilename = "Output/ComplexDietBinaryCentral/ComplexBinaryCentralCombinedGeneCorrelations"
saveRDS(binaryCors, paste0(combinedBinaryCorsFilename, ".rds"))
write.csv(binaryCors, paste0(combinedBinaryCorsFilename, ".csv"))

#---------------------------------------------------------------------
# --- Make supplementary Files --- 
# --------------------------------------------------------------------

convertColnameToIndex = function(name, df){
  output = which(colnames(df)==name)
  output
}

#Supplementary Table 1 

mergedData = read.csv("Data/MergedData.csv")

colsToKeep = c("CommonName", "ScientificNameFull", "ZoonomiaTip", "MSWC_Order", "MSWC_Family", "Divider", "SimplifiedDietConvert", "SimplifiedDietConvert", "SimplifiedDietConvertStrictPred", "Divider2", "Diet.Inv", "Diet.Vend", "Diet.Vect", "Diet.Vfish", "Diet.Vunk", "Diet.Scav", "Diet.Fruit", "Diet.Nect", "Diet.Seed", "Diet.PlantO")
colsToKeepIndex <- match(colsToKeep, names(mergedData))

SupData1 = mergedData[,colsToKeepIndex]

names(SupData1) = c("CommonName", "ScientificNameFull", "TreeTip", "MSWC_Order", "MSWC_Family", "Divider", "ThreeDietAnalysis", "FourDietAnalysis", "StrictFourDietAnalysis", "Divider2", "Diet.Inv", "Diet.Vend", "Diet.Vect", "Diet.Vfish", "Diet.Vunk", "Diet.Scav", "Diet.Fruit", "Diet.Nect", "Diet.Seed", "Diet.PlantO")

SupData1$ThreeDietAnalysis[SupData1$ThreeDietAnalysis %in% c("Vertivore", "Invertivore")] = "Carnivore"

mainAnalysisSpecies = readRDS("Output/ComplexDietCentralAnalysisSimplify2/ComplexDietCentralAnalysisSimplify2SpeciesFilter.rds")

SupData1$InMainAnalysis = (SupData1$TreeTip %in% mainAnalysisSpecies)


files <- list.files("Output/ComplexDietCentralAnalysisSimplify2/Alternates", full.names = TRUE)

alternateSpeciesFiles = files[grep("SpeciesFilter.rds", files)]

alternateSpecies = list()

for(i in 1:length(alternateSpeciesFiles)){
  currentSet = readRDS(alternateSpeciesFiles[i])
  alternateSpecies[[i]] = currentSet
}

alternateSpecies = unlist(alternateSpecies)

nameCounts <- as.data.frame(table(alternateSpecies))

names(nameCounts) <- c("Name", "Count")

#Add any missing species to the count of alternates
speciesInTableButNoAlternates = SupData1$TreeTip[which(!SupData1$TreeTip %in% nameCounts$Name)]
emptyCOunts = rep(0, length(speciesInTableButNoAlternates))
addTable = data.frame(toAdd, emptyCOunts)
names(addTable) <- c("Name", "Count")
nameCounts = rbind (nameCounts, addTable)


nameCounts = nameCounts[order(match(nameCounts$Name, SupData1$TreeTip)),]

SupData1$NumberOfIncludingAlternAlternates = nameCounts$Count

SupData1 = SupData1[,c(1:9, 21,22, 10:20)]

write.csv(SupData1, "Results/SupplementaryTable1 -- SpeciesNamesTipsAndDiet.csv")


# Supplementary Table 3


# Use DietThresholding.R script
  #Use the version that creates a compare script with no subcategories, and then trims to the 3 diet

write.csv(CompareDiets, "Results/SupplementaryTable3 -- SpeciesThresholdComparisions2.csv")



#Supplementary Data 5

outputFile = "Output/ComplexDietCentralAnalysisSimplify2/SupplementaryData5--Trees.txt" 

mainTrees = readRDS("data/zoonomiaAllMammalsTrees.rds")
masterTree = mainTrees$masterTree

fourCatTree = readRDS("Output/ComplexDietCentralAnalysisSimplify2/ComplexDietCentralAnalysisSimplify2CategoricalTree.rds")
threeCatTree = readRDS("Output/ComplexDietCentralAnalysisSimplify2/ComplexDietCentralAnalysisSimplify2MergedCategoricalTree.rds")
equalTree = readRDS("Output/ComplexDietCentralAnalysisSimplifyAltReconstruction/ComplexDietCentralAnalysisSimplifyAltReconstructionCategoricalTree.rds")
fourDietStrictTree = readRDS("Output/ComplexDietCentralAnalysisSimplifyExclusivePred/ComplexDietCentralAnalysisSimplifyExclusivePredCategoricalTree.rds")

BovTree = readRDS("Output/CladeBinaryBovidae/CladeBinaryBovidaeCategoricalTree.rds")
CerTree = readRDS("Output/CladeBinaryCervidae/CladeBinaryCervidaeCategoricalTree.rds")
CriTree = readRDS("Output/CladeBinaryCricetidae/CladeBinaryCricetidaeCategoricalTree.rds")
HysTree = readRDS("Output/CladeBinaryHystricognathi/CladeBinaryHystricognathiCategoricalTree.rds")
PerTree = readRDS("Output/CladeBinaryPeropdidae/CladeBinaryPeropdidaeCategoricalTree.rds")
VesTree = readRDS("Output/CladeBinaryVespertilionidae/CladeBinaryVespertilionidaeCategoricalTree.rds")

writeLines("----Master Tree----", outputFile)
write.tree(masterTree, file = outputFile, append = T)
cat("\n----Four Category Tree----\n", file = outputFile, append = T)
write.tree(fourCatTree, file = outputFile, append = T)
cat("\n----Three Category Tree----\n",  file =outputFile, append = T)
write.tree(threeCatTree, file = outputFile, append = T)


cat("\n----Alternate Reconstruction Tree----\n",  file =outputFile, append = T)
write.tree(equalTree, file = outputFile, append = T)
cat("\n----Strict Vertivore Invertivore Tree----\n",  file =outputFile, append = T)
write.tree(fourDietStrictTree, file = outputFile, append = T)

cat("\n\n----Clade Binary Trees----\n",  file =outputFile, append = T)
cat("\n----Bovidae Tree----\n",  file =outputFile, append = T)
write.tree(BovTree, file = outputFile, append = T)
cat("\n----Cervidae Tree----\n",  file =outputFile, append = T)
write.tree(CerTree, file = outputFile, append = T)
cat("\n----Cricetidae Tree----\n",  file =outputFile, append = T)
write.tree(CriTree, file = outputFile, append = T)
cat("\n----Hystricognathi Tree----\n",  file =outputFile, append = T)
write.tree(HysTree, file = outputFile, append = T)
cat("\n----Peropdidae Tree----\n",  file =outputFile, append = T)
write.tree(PerTree, file = outputFile, append = T)
cat("\n----Vespertilionidae Tree----\n",  file =outputFile, append = T)
write.tree(VesTree, file = outputFile, append = T)



cat("\n\n----Alternate Trees----\n",  file =outputFile, append = T)

for(i in 1:100){
  alternateTree = readRDS(paste0("Output/ComplexDietCentralAnalysisSimplify2/Alternates/Alternate", i, "ComplexDietCentralAnalysisSimplify2CategoricalTree.rds"))
  
  cat(paste("\n----Alternate", i, "Trees----\n"),  file =outputFile, append = T)
  write.tree(alternateTree, file = outputFile, append = T)
  
  
}


#Supplementary data 7 
paperAnalysis = readRDS("Output/ComplexDietCentralAnalysisSimplify2/ComplexDietCentralAnalysisSimplify2PairwiseCorrelationFile.rds")
ELAnalysis = readRDS("Output/ComplexDietCentralAnalysisSimplifyAltReconstruction/ComplexDietCentralAnalysisSimplifyAltReconstructionPairwiseCorrelationFile.rds")

paperAnalysis = paperAnalysis[-c(2,3,6,8)]
ELAnalysis = ELAnalysis[-c(2,3,6,8)]


#this function doesn't work 
compareAnalyses = function(ogAnalysis, newAnalysis, prefix = NULL, pCuttof = 0.05){
  
  if(!is.null(prefix)){
    outputOGAnalysis = list()
    outputNewAnalysis = list()
    for(i in 1:length(prefix)){
      currentCols = grep(prefix[i], names(ogAnalysis))
      currentCols = currentCols[1:3]
      
      tempOGAnalysis = ogAnalysis[,currentCols]
      tempNewAnalysis = newAnalysis[,currentCols]
      
      tempOGAnalysis = list(tempOGAnalysis)
      names(tempOGAnalysis) = prefix[i]
      
      tempNewAnalysis = list(tempNewAnalysis)
      names(tempNewAnalysis) = prefix[i]
      
      
      outputOGAnalysis = append(outputOGAnalysis, tempOGAnalysis)
      outputNewAnalysis = append(outputNewAnalysis, tempNewAnalysis)
    }
    
    newAnalysis = outputNewAnalysis
    ogAnalysis = outputOGAnalysis
  }
  
  
  
  results <- tibble(
    index = integer(),
    correlation = numeric(),
    rhoRankcorrelation = numeric(),
    numMismatchedNAs = integer(),
    pCorrelation = numeric(),
    numMismatchedPs = integer(),
    padjCorrelation = numeric(),
    numMismatchedPadjs = integer(),
    totalOgSigPadjs = integer(), 
    totalNewSigPadjs = integer(),
    missingFraction = integer(),
    mismatchedNAs = I(list()),
    mismatchedPs = I(list()),
    mismatchedPadjs = I(list())
  )
  
  for (i in seq_along(ogAnalysis)) {
    
    correlation <- cor(ogAnalysis[[i]][[1]], newAnalysis[[i]][[1]], use = "complete.obs")
    
    mismiatchedNAs <- which(!is.na(ogAnalysis[[i]][[1]]) %in% is.na(newAnalysis[[i]][[1]]))
    
    ogAnalysis[[i]] = ogAnalysis[[i]] %>% mutate(rhoRank = rank(ogAnalysis[[i]][[1]])) 
    newAnalysis[[i]] = newAnalysis[[i]] %>% mutate(rhoRank = rank(newAnalysis[[i]][[1]]))
    
    rhoRankcorrelation <- cor(ogAnalysis[[i]][[4]], newAnalysis[[i]][[4]], use = "complete.obs")
    
    pCorrelation <- cor(ogAnalysis[[i]][[2]], newAnalysis[[i]][[2]], use = "complete.obs")
    
    mismiatchedPs <- which(ogAnalysis[[i]][[2]] < pCuttof)[which(!which(ogAnalysis[[i]][[2]] < pCuttof) %in% 
                                                                   which(newAnalysis[[i]][[2]] < pCuttof))]
    
    padjCorrelation <- cor(ogAnalysis[[i]][[3]], newAnalysis[[i]][[3]], use = "complete.obs")
    
    mismiatchedPadjs <- which(ogAnalysis[[i]][[3]] < pCuttof)[which(!which(ogAnalysis[[i]][[3]] < pCuttof) %in% 
                                                                      which(newAnalysis[[i]][[3]] < pCuttof))]
    
    totalOgSigPadjs = length(which(ogAnalysis[[i]][[3]] < pCuttof))
    totalNewSigPadjs = length(which(newAnalysis[[i]][[3]] < pCuttof))
    
    missingFraction = length(mismiatchedPadjs) / length(which(ogAnalysis[[i]][[3]] < pCuttof))
    
    results <- rbind(results, tibble(
      index = i,
      correlation = correlation,
      rhoRankcorrelation = rhoRankcorrelation,
      numMismatchedNAs = length(mismiatchedNAs),
      pCorrelation = pCorrelation,
      numMissingPs = length(mismiatchedPs),
      padjCorrelation = padjCorrelation,
      numMissingPadjs = length(mismiatchedPadjs),
      totalOgSigPadjs = totalOgSigPadjs,
      totalNewSigPadjs = totalNewSigPadjs,
      missingFraction = missingFraction,
      mismatchedNAs = list(mismiatchedNAs),
      missingPs = list(mismiatchedPs),
      missingPadjs = list(mismiatchedPadjs),
    ))
  }
  row.names(results) = names(ogAnalysis)
  return(results)
}


geneMissingInNewAnalysis = compareAnalyses(paperAnalysis, newAnalysis, c("HI", "HV", "IV"))

i=1

cor(paperAnalysis[[i]][[1]], ELAnalysis[[i]][[1]], use = "complete.obs") *  cor(paperAnalysis[[i]][[1]], ELAnalysis[[i]][[1]], use = "complete.obs")

i=4
correlation = cor(paperAnalysis[[i]][[3]], ELAnalysis[[i]][[3]], use = "complete.obs") *  cor(paperAnalysis[[i]][[3]], ELAnalysis[[i]][[3]], use = "complete.obs")
correlation



#---------------------------------------------------------------------
# --- Update numbers in text --- 
# --------------------------------------------------------------------



#Getting numbers 
mainCategoricalTree = readRDS("Output/ComplexDietCentralAnalysisSimplify2/ComplexDietCentralAnalysisSimplify2CategoricalTree.rds")
mainPhenv = readRDS("Output/ComplexDietCentralAnalysisSimplify2/ComplexDietCentralAnalysisSimplify2CategoricalPhenotypeVector.rds")

table(mainPhenv)

allPhenv = readRDS("Output/ComplexDietCentralAnalysisSimplifyAllSpecies/ComplexDietCentralAnalysisSimplifyAllSpeciesCategoricalPhenotypeVector.rds")
table(allPhenv)

#mainGeneOUt 
mainGeneOut = readRDS("Output/ComplexDietCentralAnalysisSimplify2/ComplexDietCentralAnalysisSimplify2CombinedGeneResultsWithAlternates.rds")
mainGoOut = readRDS("Output/ComplexDietCentralAnalysisSimplify2/ComplexDietCentralAnalysisSimplify2CombinedGOResultsWithAlternates-KeggReactome.rds")

library(stringr)
library(dplyr)
library(tidyr)
library(ggrepel)
makeBarPlot = function(inDf, colNameSuffix = "-significantRobust") {
  
  prefixes <- str_remove(
    names(inDf)[str_ends(names(inDf), "-significantRobust")],
    "-significantRobust$"
  )
  
  plotData <- purrr::map_dfr(prefixes, function(prefix) {
    significantCol <- paste0(prefix, "-significant")
    robustCol <- paste0(prefix, "-significantRobust")
    
    tibble(
      Prefix = prefix,
      Significant = inDf[[significantCol]],
      Robust = inDf[[robustCol]]
    ) %>%
      pivot_longer(
        cols = c(Significant, Robust),
        names_to = "Type",
        values_to = "Value"
      ) %>%
      filter(!is.na(Value)) %>%
      mutate(
        Type = recode(
          Type,
          Significant = "Significant",
          Robust = "Robust significant"
        )
      )
  })
  
  # Count TRUE and FALSE results for each prefix and type
  countData <- plotData %>%
    group_by(Prefix, Type) %>%
    summarise(
      TrueCount = sum(Value == TRUE),
      FalseCount = sum(Value == FALSE),
      Total = n(),
      TruePercent = mean(Value) * 100,
      .groups = "drop"
    )
  
  # Labels showing the actual TRUE counts
  labelData <- countData %>%
    mutate(
      Label = paste0(
        "TRUE=", TrueCount,
        "\nFALSE=", FalseCount
      )
    )
  
  ggplot(plotData, aes(x = Prefix, fill = Type)) +
    geom_bar(
      aes(y = as.numeric(Value)),
      stat = "summary",
      fun = sum,
      position = position_dodge(width = 0.9)
    ) +
    geom_text(
      data = labelData,
      aes(
        x = Prefix,
        y = TrueCount,
        group = Type,
        label = paste0(TrueCount, " (", sprintf("%.1f", TruePercent), "%)")
      ),
      inherit.aes = FALSE,
      position = position_dodge(width = 0.9),
      vjust = -0.5
    ) +
    labs(
      x = "Prefix",
      y = "Number of TRUE results",
      fill = "Column",
      title = colNameSuffix
    ) +
    scale_fill_manual(
      values = c(
        "Significant" = "steelblue",
        "Robust significant" = "tomato"
      )
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.15))
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}


makeBarPlot(mainGoOut)

CompareDiets = mainGoOut
makeBarPlot()


makeBarPlot = function(inDf, prefixesToUse = NULL, labelOffset = 0.12, labelThreshold = 0.15, colNameSuffix = "-significantRobust", xAxisValues = NULL) {
  
  # Identify available prefixes
  availablePrefixes <- sub(
    "-significantRobust$",
    "",
    names(inDf)[endsWith(names(inDf), "-significantRobust")]
  )
  
  # Use selected prefixes, or all available prefixes
  if (is.null(prefixesToUse)) {
    prefixes <- availablePrefixes
  } else {
    prefixes <- prefixesToUse
    
    missingPrefixes <- setdiff(prefixes, availablePrefixes)
    
    if (length(missingPrefixes) > 0) {
      stop(
        "Missing -significantRobust columns for prefixes: ",
        paste(missingPrefixes, collapse = ", ")
      )
    }
  }
  
  # Manually set x-axis order and displayed values
  if (is.null(xAxisValues)) {
    xAxisValues <- prefixes
  } else {
    if (length(xAxisValues) != length(prefixes)) {
      stop("xAxisValues must have the same length as prefixes.")
    }
    
    names(xAxisValues) <- prefixes
  }
  
  
  # Prepare data
  plotData <- purrr::map_dfr(prefixes, function(prefix) {
    
    significantCol <- paste0(prefix, "-significant")
    robustCol <- paste0(prefix, "-significantRobust")
    
    if (!significantCol %in% names(inDf)) {
      stop("Missing column: ", significantCol)
    }
    
    significantTotal <- sum(inDf[[significantCol]] == TRUE, na.rm = TRUE)
    robustTotal <- sum(inDf[[robustCol]] == TRUE, na.rm = TRUE)
    
    tibble::tibble(
      Prefix = prefix,
      `Robust significant` = robustTotal,
      `Non-robust significant` = significantTotal - robustTotal
    ) %>%
      tidyr::pivot_longer(
        cols = c(`Robust significant`, `Non-robust significant`),
        names_to = "Type",
        values_to = "TrueCount"
      )
  })
  
  countData <- plotData %>%
    dplyr::group_by(Prefix) %>%
    dplyr::mutate(
      PrefixTotal = sum(TrueCount),
      Percent = dplyr::if_else(
        PrefixTotal > 0,
        100 * TrueCount / PrefixTotal,
        0
      ),
      Label = paste0(
        TrueCount, "\n(",
        sprintf("%.1f", Percent), "%)"
      )
    ) %>%
    dplyr::ungroup()
  
  # Threshold based on the tallest total bar
  largestColumn <- max(
    countData %>%
      dplyr::group_by(Prefix) %>%
      dplyr::summarise(Total = sum(TrueCount), .groups = "drop") %>%
      dplyr::pull(Total),
    na.rm = TRUE
  )
  
  labelCutoff <- largestColumn * labelThreshold
  
  # Plot stacked bars
  ggplot2::ggplot(
    countData,
    ggplot2::aes(x = Prefix, y = TrueCount, fill = Type)
  ) +
    ggplot2::geom_col() +
    # Labels inside sufficiently large segments
    ggplot2::geom_text(
      data = dplyr::filter(
        countData,
        TrueCount > 0,
        TrueCount >= labelCutoff
      ),
      ggplot2::aes(label = Label),
      position = ggplot2::position_stack(vjust = 0.5),
      color = "white",
      fontface = "bold",
      size = 3.5
    ) +
    
    # Labels for small segments, positioned above the stacked bar
    ggrepel::geom_text_repel(
      data = dplyr::filter(
        countData,
        TrueCount > 0,
        TrueCount < labelCutoff
      ) %>%
        dplyr::group_by(Prefix) %>%
        dplyr::mutate(
          BarTotal = sum(TrueCount),
          LabelY = BarTotal + largestColumn * labelOffset
        ) %>%
        dplyr::ungroup(),
      ggplot2::aes(
        x = Prefix,
        y = LabelY,
        label = Label
      ),
      inherit.aes = FALSE,
      direction = "y",
      box.padding = 0.5,
      point.padding = 0.2,
      force = 2,
      max.overlaps = Inf,
      color = "black",
      fontface = "bold",
      size = 3.5,
      segment.color = "grey50"
    ) +
    ggplot2::labs(
      x = "Comparision",
      y = "Number of significant results",
      fill = "Robustness",
      title = colNameSuffix
    ) +
    ggplot2::scale_x_discrete(
      limits = prefixes,
      labels = xAxisValues
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "Non-robust significant" = "steelblue",
        "Robust significant" = "tomato"
      )
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.08))
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      )
    )
}



chosenPrefixes = c("CH", "HI", "HV", "IV")
xAxisLabs = c("CombinedCarnivore-Herbivore", "Invertivore-Herbivore", "Vertivore-Herbivore", "Invertivore-Vertivore")

geneRobustPlot = makeBarPlot(mainGeneOut, chosenPrefixes, 0.15, 0.1, "Genes", xAxisValues =xAxisLabs)
goRobustPlot = makeBarPlot(mainGoOut, chosenPrefixes, 0.20, 0.15, "K-R Categories", xAxisValues= xAxisLabs)

mean(mainGeneOut$`CH-PadjNumSignificant`[which(mainGeneOut$`CH-significant`)])
mean(mainGeneOut$`HI-PadjNumSignificant`[which(mainGeneOut$`HI-significant`)])
mean(mainGeneOut$`HV-PadjNumSignificant`[which(mainGeneOut$`HV-significant`)])
mean(mainGeneOut$`IV-PadjNumSignificant`[which(mainGeneOut$`IV-significant`)])

mean(mainGoOut$`CH-PadjNumSignificant`[which(mainGoOut$`CH-significant`)])
mean(mainGoOut$`HI-PadjNumSignificant`[which(mainGoOut$`HI-significant`)])
mean(mainGoOut$`HV-PadjNumSignificant`[which(mainGoOut$`HV-significant`)])
mean(mainGoOut$`IV-PadjNumSignificant`[which(mainGoOut$`IV-significant`)])

barplotFilename = ("Output/ComplexDietCentralAnalysisSimplify2/Figure S3 -- Robustness to Alternates.png")

png(barplotFilename, height = 6, width = 12, units = "in", res = 216)
grid.arrange(geneRobustPlot, goRobustPlot, ncol = 2)
grid.text("A", x = unit(0.02, "npc"), y = unit(0.99, "npc"), just = c("left", "top"), gp = gpar(fontsize = 24, fontface = "bold"))
grid.text("B", x = unit(0.5, "npc"), y = unit(0.99, "npc"), just = c("left", "top"), gp = gpar(fontsize = 24, fontface = "bold"))
dev.off()

length(which(mainGeneOut$`IV-significantRobust`))

length(which(mainGoOut$`CH-significant`))


strictGeneOut = readRDS("Output/ComplexDietCentralAnalysisSimplifyExclusivePred/ComplexDietCentralAnalysisSimplifyExclusivePredCombinedGeneResultsWithAlternates.rds")
strictGoOut = readRDS("Output/ComplexDietCentralAnalysisSimplifyExclusivePred/ComplexDietCentralAnalysisSimplifyExclusivePredCombinedGOResultsWithAlternates-KeggReactome.rds")

mainIVGenes = rownames(mainGeneOut)[which(mainGeneOut$`IV-significantRobust`)]
mainIVGo = rownames(mainGoOut)[which(mainGoOut$`IV-significantRobust`)]

strictIVGenes = rownames(strictGeneOut)[which(strictGeneOut$`IV-significantRobust`)]
strictIVGo = rownames(strictGoOut)[which(strictGoOut$`IV-significantRobust`)]

overlapIVGenes = mainIVGenes[which(mainIVGenes %in% strictIVGenes)]
overlapIVGo = mainIVGo[which(mainIVGo %in% strictIVGo)]


mainHIGenes = rownames(mainGeneOut)[which(mainGeneOut$`HI-significantRobust`)]
mainHIGo = rownames(mainGoOut)[which(mainGoOut$`HI-significantRobust`)]

strictHIGenes = rownames(strictGeneOut)[which(strictGeneOut$`HI-significantRobust`)]
strictHIGo = rownames(strictGoOut)[which(strictGoOut$`HI-significantRobust`)]

overlapHIGenes = mainHIGenes[which(mainHIGenes %in% strictHIGenes)]
overlapHIGo = mainHIGo[which(mainHIGo %in% strictHIGo)]


mainHVGenes = rownames(mainGeneOut)[which(mainGeneOut$`HV-significantRobust`)]
mainHVGo = rownames(mainGoOut)[which(mainGoOut$`HV-significantRobust`)]

strictHVGenes = rownames(strictGeneOut)[which(strictGeneOut$`HV-significantRobust`)]
strictHVGo = rownames(strictGoOut)[which(strictGoOut$`HV-significantRobust`)]

overlapHVGenes = mainHIGenes[which(mainHVGenes %in% strictHVGenes)]
overlapHVGo = mainHIGo[which(mainHVGo %in% strictHVGo)]

mainCHGenes = rownames(mainGeneOut)[which(mainGeneOut$`CH-significantRobust`)]
mainCHGo = rownames(mainGoOut)[which(mainGoOut$`CH-significantRobust`)]

strictCHGenes = rownames(strictGeneOut)[which(strictGeneOut$`CH-significantRobust`)]
strictCHGo = rownames(strictGoOut)[which(strictGoOut$`CH-significantRobust`)]

overlapCHGenes = mainCHGenes[which(mainCHGenes %in% strictCHGenes)]
overlapCHGo = mainCHGo[which(mainCHGo %in% strictCHGo)]


geneCompare = list(mainIVGenes, strictIVGenes, overlapIVGenes, mainIVGo, strictIVGo, overlapIVGo, mainHIGenes, strictHIGenes, overlapHIGenes, mainHIGo, strictHIGo, overlapHIGo, mainHVGenes, strictHVGenes, overlapHVGenes, mainHVGo, strictHVGo, overlapHVGo, mainCHGenes, strictCHGenes, overlapCHGenes, mainCHGo, strictCHGo, overlapCHGo)
names(geneCompare) = c("mainIVGenes", "strictIVGenes", "overlapIVGenes", "mainIVCategories", "strictIVCategories", "overlapIVCategories", "mainHIGenes", "strictHIGenes", "overlapHIGenes", "mainHICategories", "strictHICategories", "overlapHICategories", "mainHVGenes", "strictHVGenes", "overlapHVGenes", "mainHVCategories", "strictHVCategories", "overlapHVCategories", "mainCHGenes", "strictCHGenes", "overlapCHGenes", "mainCHCategories", "strictCHCategories", "overlapCHCategories")

makeMainStrictVenn = function(analysisName, compareList = geneCompare){
  relevantResults = grep(analysisName, names(compareList))
  names(relevantResults) = names(compareList)[grep(analysisName, names(compareList))]
  
  mainAnalysisIndex = relevantResults[grep("main", names(relevantResults))]
  strictAnalysisIndex = relevantResults[grep("strict", names(relevantResults))]
  
  venn_list <- list(
    Main = compareList[[mainAnalysisIndex]],
    Strict = compareList[[strictAnalysisIndex]]
  )
  titleName = paste(analysisName)
  vennout = ggvenn(venn_list, set_name_size = 4, text_size = 4) +  labs(title = titleName)+
    coord_cartesian(clip = "off") +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),
      plot.margin = margin(t = 40, r = 10, b = 10, l = 10)
    )
  vennout
}

CHGeneVenn = makeMainStrictVenn("CHGenes")
HVGeneVenn = makeMainStrictVenn("HIGenes")
HIGeneVenn = makeMainStrictVenn("HVGenes")
IVGeneVenn = makeMainStrictVenn("IVGenes")

GeneVenn = grid.arrange(CHGeneVenn, HVGeneVenn, HIGeneVenn, IVGeneVenn)

CHCategoriesVenn = makeMainStrictVenn("CHCategories")
HVCategoriesVenn = makeMainStrictVenn("HICategories")
HICategoriesVenn = makeMainStrictVenn("HVCategories")
IVCategoriesVenn = makeMainStrictVenn("IVCategories")

CategoriesVenn = grid.arrange(CHCategoriesVenn, HVCategoriesVenn, HICategoriesVenn, IVCategoriesVenn)

mainStictVennFilename = "Output/ComplexDietCentralAnalysisSimplify2/Supplmentary Figure 4 -- Overlap with Strict Analysis.png"
png(mainStictVennFilename, height = 6, width = 12, units = "in", res = 216)
grid.arrange(GeneVenn, CategoriesVenn, ncol = 2)
dev.off()

colnames(strictGeneOut) = paste0(colnames(strictGeneOut), "-S")
combinedGeneOut = cbind (mainGeneOut, strictGeneOut)
colnames(strictGoOut) = paste0(colnames(strictGoOut), "-S")
combinedGoOut = cbind (mainGoOut, strictGoOut)
colnames(combinedGeneOut)
plot(combinedGeneOut$`CH-Rho`, combinedGeneOut$`CH-Rho-S`)
plot(combinedGeneOut$`HI-Rho`, combinedGeneOut$`HI-Rho-S`)
plot(combinedGeneOut$`HV-Rho`, combinedGeneOut$`HV-Rho-S`)
plot(combinedGeneOut$`IV-Rho`, combinedGeneOut$`IV-Rho-S`)


makeCorPlot = function(x, y, prefix, titleVal){

plot(x, y,
     xlab = paste0(prefix," Stat Main"),
     ylab = paste0(prefix," Stat Strict"),
     main = titleVal
     )

# Calculate R-squared
model <- lm(y ~ x)
rSquared <- summary(model)$r.squared

# Add R-squared label
legend("topleft",
       legend = bquote(R^2 == .(round(rSquared, 3))),
       bty = "n")
}
par(mfrow = c(1,1))
makeCorPlot(combinedGeneOut$`CH-Rho`, combinedGeneOut$`CH-Rho-S`, "C-H", "Combined Carnivore-Herbivore Genes")
makeCorPlot(combinedGeneOut$`HI-Rho`, combinedGeneOut$`HI-Rho-S`, "H-I", "Invertivore-Herbivore Genes")
makeCorPlot(combinedGeneOut$`HV-Rho`, combinedGeneOut$`HV-Rho-S`, "H-V", "Vertivore-Herbivore Genes")
makeCorPlot(combinedGeneOut$`IV-Rho`, combinedGeneOut$`IV-Rho-S`, "I-V", "Invertivore-Vertivore Genes")

makeCorPlot(combinedGoOut$`CH-stat`, combinedGoOut$`CH-stat-S`, "C-H", "Combined Carnivore-Herbivore Categories")
makeCorPlot(combinedGoOut$`HI-stat`, combinedGoOut$`HI-stat-S`, "H-I", "Invertivore-Herbivore Categories")
makeCorPlot(combinedGoOut$`HV-stat`, combinedGoOut$`HV-stat-S`, "H-V", "Vertivore-Herbivore Categories")
makeCorPlot(combinedGoOut$`IV-stat`, combinedGoOut$`IV-stat-S`, "I-V", "Invertivore-Vertivore Categories")

mainStictCorFilename = "Output/ComplexDietCentralAnalysisSimplify2/Supplmentary Figure 4 -- Correlation with Strict Analysis.png"
png(mainStictCorFilename, height = 6, width = 12, units = "in", res = 216)
par(mfrow = c(2,4))
makeCorPlot(combinedGeneOut$`CH-Rho`, combinedGeneOut$`CH-Rho-S`, "C-H", "Combined Carnivore-Herbivore Genes")
makeCorPlot(combinedGeneOut$`HI-Rho`, combinedGeneOut$`HI-Rho-S`, "H-I", "Invertivore-Herbivore Genes")

makeCorPlot(combinedGoOut$`CH-stat`, combinedGoOut$`CH-stat-S`, "C-H", "Combined Carnivore-Herbivore Categories")
makeCorPlot(combinedGoOut$`HI-stat`, combinedGoOut$`HI-stat-S`, "H-I", "Invertivore-Herbivore Categories")

makeCorPlot(combinedGeneOut$`HV-Rho`, combinedGeneOut$`HV-Rho-S`, "H-V", "Vertivore-Herbivore Genes")
makeCorPlot(combinedGeneOut$`IV-Rho`, combinedGeneOut$`IV-Rho-S`, "I-V", "Invertivore-Vertivore Genes")

makeCorPlot(combinedGoOut$`HV-stat`, combinedGoOut$`HV-stat-S`, "H-V", "Vertivore-Herbivore Categories")
makeCorPlot(combinedGoOut$`IV-stat`, combinedGoOut$`IV-stat-S`, "I-V", "Invertivore-Vertivore Categories")

dev.off()



overlapHivGenes = mainHIGenes[which(mainHIGenes %in% mainHVGenes)]
length(overlapHivGenes)/length(mainHVGenes)
length(overlapHivGenes)/length(mainHIGenes)

overlapHivGo = mainHIGo[which(mainHIGo %in% mainHVGo)]
length(overlapHivGo)/length(mainHVGo)
length(overlapHivGo)/length(mainHIGo)





maiHVSigOnly = mainGoOut[which(mainGoOut$`HV-significantRobust`),]

maiHVSigOnly[,c(2,52,102, 105)]

which(rownames(mainGoOut) == "KEGG_TRYPTOPHAN_METABOLISM")
test = mainGoOut[1676,]

#Compare GO patterns 

aminoAcid = c("REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",
"KEGG_GLYCINE_SERINE_AND_THREONINE_METABOLISM",
"KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION",
"KEGG_TYROSINE_METABOLISM",
"KEGG_TRYPTOPHAN_METABOLISM",
"REACTOME_BRANCHED_CHAIN_AMINO_ACID_CATABOLISM",
"REACTOME_GLYOXYLATE_METABOLISM_AND_GLYCINE_DEGRADATION")

lipid = c("REACTOME_METABOLISM_OF_LIPIDS",
"REACTOME_FATTY_ACID_METABOLISM",
"KEGG_FATTY_ACID_METABOLISM",
"KEGG_STEROID_HORMONE_BIOSYNTHESIS",
"REACTOME_METABOLISM_OF_STEROIDS",
"REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION")

bile = c("REACTOME_BILE_ACID_AND_BILE_SALT_METABOLISM",
"REACTOME_RECYCLING_OF_BILE_ACIDS_AND_SALTS",
"REACTOME_SYNTHESIS_OF_BILE_ACIDS_AND_BILE_SALTS_VIA_7ALPHA_HYDROXYCHOLESTEROL")

peroxisome = c("KEGG_PEROXISOME",
"REACTOME_PEROXISOMAL_LIPID_METABOLISM",
"REACTOME_PEROXISOMAL_PROTEIN_IMPORT")

vitamin = c("REACTOME_METABOLISM_OF_WATER_SOLUBLE_VITAMINS_AND_COFACTORS",
"REACTOME_DEFECTS_IN_VITAMIN_AND_COFACTOR_METABOLISM",
"REACTOME_COBALAMIN_CBL_VITAMIN_B12_TRANSPORT_AND_METABOLISM",
"REACTOME_UPTAKE_OF_DIETARY_COBALAMINS_INTO_ENTEROCYTES",
"REACTOME_CHOLINE_CATABOLISM",
"REACTOME_METABOLISM_OF_VITAMINS_AND_COFACTORS")

etc = c("REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
"REACTOME_COMPLEX_I_BIOGENESIS")
oxidation = c("REACTOME_BIOLOGICAL_OXIDATIONS",
"REACTOME_METABOLIC_DISORDERS_OF_BIOLOGICAL_OXIDATION_ENZYMES")

drug = c("KEGG_DRUG_METABOLISM_CYTOCHROME_P450",
"REACTOME_CYTOCHROME_P450_ARRANGED_BY_SUBSTRATE_TYPE",
"REACTOME_PHASE_I_FUNCTIONALIZATION_OF_COMPOUNDS",
"REACTOME_PHASE_II_CONJUGATION_OF_COMPOUNDS")

sulfur = c("KEGG_SULFUR_METABOLISM",
"REACTOME_CYTOSOLIC_SULFONATION_OF_SMALL_MOLECULES")

transport = c("REACTOME_AMINO_ACID_TRANSPORT_ACROSS_THE_PLASMA_MEMBRANE",
"REACTOME_SLC_MEDIATED_TRANSMEMBRANE_TRANSPORT",
"REACTOME_TRANSPORT_OF_INORGANIC_CATIONS_ANIONS_AND_AMINO_ACIDS_OLIGOPEPTIDES")

allCategories <- c(
  setNames(aminoAcid, rep("aminoAcid", length(aminoAcid))),
  setNames(lipid, rep("lipid", length(lipid))),
  setNames(bile, rep("bile", length(bile))),
  setNames(peroxisome, rep("peroxisome", length(peroxisome))),
  setNames(vitamin, rep("vitamin", length(vitamin))),
  setNames(etc, rep("etc", length(etc))),
  setNames(oxidation, rep("oxidation", length(oxidation))),
  setNames(drug, rep("drug", length(drug))),
  setNames(sulfur, rep("sulfur", length(sulfur))),
  setNames(transport, rep("transport", length(transport)))
)


categoryResults = mainGoOut[which(rownames(mainGoOut) %in% allCategories), c(convertColnameToIndex("CH-stat", mainGoOut), convertColnameToIndex("HI-stat", mainGoOut), convertColnameToIndex("HV-stat", mainGoOut)) ]
categoryResults$cluster = names(allCategories)[match(rownames(categoryResults), allCategories)]

categoryResults$ratio = categoryResults$`HV-stat` /categoryResults$`HI-stat`
categoryResults$VCratio = -1*categoryResults$`HV-stat` /categoryResults$`CH-stat`

hist(categoryResults$ratio)
hist(
  categoryResults$ratio,
  breaks = c(seq(0, 1, by = 0.1),2),
  right = FALSE,
  xlab = "Ratio",
  main = "Distribution of Ratios"
)

# categoryResults = categoryResults[order(categoryResults$cluster),]

cluster_colors <- c(
  aminoAcid = "#E41A1C",
  lipid      = "#984EA3",
  bile       = "#4DAF4A",
  vitamin    = "blue",
  drug       = "#F781BF",
  
  oxidation  = "#FF7F00",
  peroxisome = "#A65628",
  sulfur     = "yellow", 
  transport  = "#00A6D6",
  
  etc        = "#555555"
)

ggplot(categoryResults, aes(x = `HI-stat`, y = `HV-stat`, color = cluster)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = cluster_colors) +
  labs(
    x = "HI Statistic",
    y = "HV Statistic",
    color = "Cluster"
  ) +
  theme_classic() +
  theme(
    legend.position = "right"
  )+
  coord_equal() +
  geom_smooth(
    aes(group = 1),
    method = "lm",
    se = FALSE,
    color = "black",
    linetype = "solid"
  ) 


clusterCorrelations <- categoryResults %>%
  group_by(cluster) %>%
  summarise(
    n = sum(complete.cases(`HI-stat`, `HV-stat`)),
    correlation = if (n >= 2 &&
                      sd(`HI-stat`, na.rm = TRUE) > 0 &&
                      sd(`HV-stat`, na.rm = TRUE) > 0) {
      cor(`HI-stat`, `HV-stat`, use = "complete.obs")
    } else {
      NA_real_
    },
    pValue = if (n >= 3 &&
                 sd(`HI-stat`, na.rm = TRUE) > 0 &&
                 sd(`HV-stat`, na.rm = TRUE) > 0) {
      cor.test(`HI-stat`, `HV-stat`)$p.value
    } else {
      NA_real_
    },
    .groups = "drop"
  )

cor(animoAcidResults$`HI-stat`, animoAcidResults$`HV-stat`)
plot(animoAcidResults$`HI-stat`, animoAcidResults$`HV-stat`)
dev.off()

#---------------------------------------------------------------------
# --- Comapre main vs exclusive --- 
# --------------------------------------------------------------------
library(tibble)
library(ggvenn)
library(gridExtra)

mainGeneData = readRDS("Output/ComplexDietCentralAnalysisSimplify2/ComplexDietCentralAnalysisSimplify2combinedGeneResultsWithAlternates.rds")
excluGeneData = readRDS("Output/ComplexDietCentralAnalysisSimplifyExclusivePred/ComplexDietCentralAnalysisSimplifyExclusivePredcombinedGeneResultsWithAlternates.rds")


all(rownames(mainGeneData) == rownames(excluGeneData))
all(names(mainGeneData) == names(excluGeneData))

#CH
currentCol = which(names(mainGeneData) == "CH-P")
  #Cor = 0.89
currentCol = which(names(mainGeneData) == "HI-P")
  #Cor = 0.88
currentCol = which(names(mainGeneData) == "HV-P")
  #Cor = 0.72
currentCol = which(names(mainGeneData) == "IV-P")
  #Cor = 0.72


mainLogBaseP = log10(mainGeneData[currentCol]) * sign(mainGeneData[(currentCol -7)]) *-1
excluLogBaseP = log10(excluGeneData[currentCol]) * sign(excluGeneData[(currentCol -7)]) *-1
cor(mainLogBaseP, excluLogBaseP, method = "spearman", use = "complete.obs")




mainGOData = readRDS("Output/ComplexDietCentralAnalysisSimplify2/ComplexDietCentralAnalysisSimplify2combinedGOResultsWithAlternates-KeggReactome.rds")
excluGOData = readRDS("Output/ComplexDietCentralAnalysisSimplifyExclusivePred/ComplexDietCentralAnalysisSimplifyExclusivePredcombinedGOResultsWithAlternates-KeggReactome.rds")


all(rownames(mainGOData) == rownames(excluGOData))
all(names(mainGOData) == names(excluGOData))


#CH
currentCol = which(names(mainGOData) == "CH-pval")
#Cor = 0.89
currentCol = which(names(mainGOData) == "HI-pval")
#Cor = 0.93
currentCol = which(names(mainGOData) == "HV-pval")
#Cor = 0.62
currentCol = which(names(mainGOData) == "IV-pval")
#Cor = 0.76


mainLogBaseP = log10(mainGOData[currentCol]) * sign(mainGOData[(currentCol -7)]) *-1
excluLogBaseP = log10(excluGOData[currentCol]) * sign(excluGOData[(currentCol -7)]) *-1
cor(mainLogBaseP, excluLogBaseP, method = "spearman", use = "complete.obs")



plot1 = plot(signedLogBaseP, signedLogPermP)
plot2 = plot(signedLogBaseP, signedLogPermP)
plot3 = plot(signedLogBaseP, signedLogPermP)

plot1




ChOverlapGO = which(mainGOData$`CH-significantRobust` & excluGOData$`CH-significantRobust` & mainGOData$`CH-stat` <0)
rownames(mainGOData)[ChOverlapGO]

HiOverlapGO = which(mainGOData$`HI-significantRobust` & excluGOData$`HI-significantRobust` & mainGOData$`HI-stat` >0)
rownames(mainGOData)[HiOverlapGO]

HvOverlapGO = which(mainGOData$`HV-significantRobust` & excluGOData$`HV-significantRobust` & mainGOData$`HV-stat` >0)
rownames(mainGOData)[HvOverlapGO]

IvOverlapGO = which(mainGOData$`IV-significantRobust` & excluGOData$`IV-significantRobust`)
rownames(mainGOData)[IvOverlapGO]


#---------------------------------------------------------------------
# --- DebugGoOrdering --- 
# --------------------------------------------------------------------



tempGeneset1 = data.frame(a=c(1,2,3), b=c(1,2,3))
rownames(tempGeneset1) = c("one", "two", "three")
tempGeneset2 = data.frame(a=c(1,4), b=c(1,4))
rownames(tempGeneset2) = c("one", "four")
rownames(tempGeneset2) = c("one", "three")

GOResults = list(tempGeneset1, tempGeneset2)
mainRowNames = rownames(EnrichmentRuns[[1]][[1]])

errorSet = NA
for(i in 1:100){
  for(j in 1:length(EnrichmentRuns[[i]])){
    newRowNames = rownames(EnrichmentRuns[[i]][[j]]) 
    isEqual = all(newRowNames ==mainRowNames)
    if(!isEqual){
      errorCode = paste0(i, "-", j)
      names(errorCode) = i
      errorSet = append(errorSet, errorCode)
    }
  }
}
length(unique(names(errorSet)))

#This confirms that at least the row names are always consisent within an alternate
errorSet = NA
for(i in 1:100){
  mainRowNames = rownames(EnrichmentRuns[[i]][[1]])
  for(j in 2:length(EnrichmentRuns[[i]])){
    newRowNames = rownames(EnrichmentRuns[[i]][[j]]) 
    isEqual = all(newRowNames ==mainRowNames)
    if(!isEqual){
      errorCode = paste0(i, "-", j)
      names(errorCode) = i
      errorSet = append(errorSet, errorCode)
    }
  }
}


length(unique(names(errorSet)))

newRowNames
mainRowNames

which(!newRowNames == mainRowNames)

combinedEnrichment[[1]][[1]][[279]]

testEnrichment = EnrichmentRuns[[i]][[j]]
testEnrichment$stat[[279]]

#---------------------------------------------------------------------
# --- compare working and broken alternates--- 
# --------------------------------------------------------------------
brokenAlternatesFilename = paste(outputFolderName, "Alternates/Broken/", filePrefix, "BrokenAlternateSets.rds", sep="")
brokenAlternates = readRDS(brokenAlternatesFilename)
alternateSets = readRDS(paste0(outputFolderName, filePrefix, "AlternatePruningSpecies.rds"))

workingAlternates = alternateSets


brokenVector = unlist(brokenAlternates)
mainVector = unlist(workingAlternates)

speciesInBroken = unique(brokenVector)
speciesInMain = unique(mainVector)

which(!speciesInBroken %in% speciesInMain)
#It's not a particular species that breaks things, then. 
#Yeah, I have no idea. COmparing the working and the broken ones was my best bet, but other than somehow finding pairs or something, or a better was to visualize the borken paths, 


table(brokenVector)
table(mainVector)


pathsList = list()
j=0
for(i in names(brokenAlternates)){
  message(i)
  j=j+1
  currentSet = brokenAlternates[which(names(brokenAlternates) ==i)]
  alternateFilePrefix = paste0("/Alternates/Alternate", i)
  currentCategoricalTreeFilename = paste(outputFolderName, alternateFilePrefix, filePrefix, "CategoricalTree.rds", sep="") #make a filename based on the prefix
  currentCommonCategoricalTreeFilename = paste(outputFolderName, alternateFilePrefix, filePrefix, "CategoricalCommonTree.rds", sep="") #make a filename based on the prefix
  currentPhenotypeVectorFilename = paste(outputFolderName, alternateFilePrefix, filePrefix, "CategoricalPhenotypeVector.rds",sep="") #make a filename based on the prefix
  currentSpeciesFilterFilename =  paste(outputFolderName, alternateFilePrefix, filePrefix, "SpeciesFilter.rds",sep="") #set a filename for the species filter based on the prefix 
  
  currentPhenotypeTree = readRDS(currentCategoricalTreeFilename)
  currentSpeciesFilter = readRDS(currentSpeciesFilterFilename)
  
  currentPathsObject = tree2Paths(currentPhenotypeTree, mainTrees, useSpecies = currentSpeciesFilter, categorical = TRUE)
  table(currentPathsObject)
  pathsList[[j]] = currentPathsObject
}
names(pathsList) = names(brokenAlternates)
  

table(pathsList[[1]])

#And now they don't seem to be broken?


#---------------------------------------------------------------------
# --- Make MergedTree plots--- 
# --------------------------------------------------------------------

mergedTree = readRDS("Output/ComplexDietCentralAnalysisSimplify2/ComplexDietCentralAnalysisSimplify2MergedCategoricalTree.rds")
mergedAllTree = readRDS("Output/ComplexDietCentralAnalysisSimplifyAllSpecies/ComplexDietCentralAnalysisSimplifyAllSpeciesMergedCategoricalTree.rds")



palette(c("orange", "darkgreen", "black"))


mergedCommonTree = mergedTree
mergedCommonTree$tip.label = ZonomNameConvertVectorCommon(mergedTree$tip.label, tipCol = "ZoonomiaTip")
treeImageFilename = "Output/ComplexDietCentralAnalysisSimplify2/ComplexDietCentralAnalysisSimplify2MergedCategoricalTree.pdf"
pdf(treeImageFilename, height = length(mergedTree$tip.label)/18, width = 10)                     #make a pdf to store the plot, sized based on tree size

plotTreeCategorical(mergedCommonTree, c("Carnivore", "Herbivore", "Omnivore"), master = commonMainTrees$masterTree)
plotTreeCategorical(mergedTree, c("Carnivore", "Herbivore", "Omnivore"), master = mainTrees$masterTree)
dev.off()

treeImageFilename = "Output/ComplexDietCentralAnalysisSimplifyAllSpecies/ComplexDietCentralAnalysisSimplifyAllSpeciesMergedCategoricalTree.pdf"

mergedAllCommonTree = mergedAllTree
mergedAllCommonTree$tip.label = ZonomNameConvertVectorCommon(mergedAllTree$tip.label, tipCol = "ZoonomiaTip")


pdf(treeImageFilename, height = length(mergedTree$tip.label)/18, width = 10)                     #make a pdf to store the plot, sized based on tree size

plotTreeCategorical(mergedAllCommonTree, c("Carnivore", "Herbivore", "Omnivore"), master = commonMainTrees$masterTree)
plotTreeCategorical(mergedAllTree, c("Carnivore", "Herbivore", "Omnivore"), master = mainTrees$masterTree)
dev.off()


#---------------------------------------------------------------------
# --- debug StrictPred 2--- 
# --------------------------------------------------------------------
treeNames = phenotypeTree$tip.label
treeNames = treeNames[order(treeNames)]

setNames = currentSet
setNames = setNames[order(setNames)]

filterNames = speciesFilter
filterNames = filterNames[order(filterNames)]

namesFrame = data.frame(treeNames, setNames, filterNames)

which(!namesFrame[[1]] == namesFrame[[2]])
all.equal(namesFrame[[1]], namesFrame[[2]])
all.equal(namesFrame[[1]], namesFrame[[3]])


phenotypeTree$tip.label %in% currentSet
currentSet %in% phenotypeTree$tip.label

mainTrees$masterTree$tip.label %in% speciesFilter

phenotypeTree$edge.length
currentSet %in% speciesFilter

unique(colnames(RERObject))


table(phenotypeTree$edge.length)
table(pathsObject)


RERObject = getAllResiduals(mainTrees, useSpecies = speciesFilter, plot = F, min.sp=10)  #Calculate the RERs

length(rownames(RERObject))


write.csv(CompareDiets, "Results/CompareDiets.csv")



#---------------------------------------------------------------------
# --- debug StrictPred --- 
# --------------------------------------------------------------------

startingName = "Output/ComplexDietCentralAnalysisSimplifyStrictPred/Alternates/Alternate9ComplexDietCentralAnalysisSimplifyStrictPred"

alternateSets = readRDS("Output/ComplexDietCentralAnalysisSimplifyStrictPred/ComplexDietCentralAnalysisSimplifyStrictPredAlternatePruningSpecies.rds")
alternate = alternateSets[[9]]
speciesFilter = readRDS(paste0(startingName, "SpeciesFilter.rds"))
phenVec = readRDS(paste0(startingName, "CategoricalPhenotypeVector.rds"))
length(phenVec)
tree = readRDS(paste0(startingName, "CategoricalTree.rds"))
length(tree$tip.label)

which(!alternate %in% speciesFilter)
which(!speciesFilter %in% alternate)
which(!speciesFilter %in% names(phenVec))
which(! names(phenVec)%in% speciesFilter)
table(tree$edge.length)

vecFromCluster = c("vs_HLornAna3", "vs_monDom5", "vs_HLdidVir1", "vs_HLgraAgi1",
                   "vs_HLantFla1", "vs_HLphaCin1", "vs_HLvomUrs1", "vs_HLgymLea1",
                   "vs_HLpseCor1", "vs_HLpseCup1", "vs_HLphaGym1", "vs_HLmacGig1",
                   "vs_oryAfe1", "vs_chrAsi1", "vs_HLmicTal1", "vs_HLloxAfr4", "vs_HLtamTet1",
                   "vs_HLnycCou1", "vs_otoGar3", "vs_HLdauMad1", "vs_HLindInd1",
                   "vs_proCoq1", "vs_HLcheMed1", "vs_HLeulMon1", "vs_HLpitPit1",
                   "vs_HLateGeo1", "vs_HLcebAlb1", "vs_HLsagImp1", "vs_HLcalPym1",
                   "vs_HLcalJac4", "vs_ponAbe3", "vs_HLpilTep2", "vs_colAng1", "vs_HLsemEnt1",
                   "vs_HLallNig1", "vs_HLcerNeg1", "vs_HLcerMon1", "vs_HLpapAnu5",
                   "vs_macNem1", "vs_ochPri3", "vs_HLoryCun3", "vs_HLmusAve1", "vs_HLsciVul1",
                   "vs_HLxerIna1", "vs_HLspeDau1", "vs_HLmarFla1", "vs_HLcteGun1",
                   "vs_HLhysCri1", "vs_HLpetTyp1", "vs_hetGla2", "vs_chiLan1", "vs_HLdinBra1",
                   "vs_HLmyoCoy1", "vs_HLcteSoc1", "vs_HLcoePre1", "vs_HLdolPat1",
                   "vs_HLhydHyd1", "vs_HLpedCap1", "vs_HLperLonPac1", "vs_dipOrd2",
                   "vs_HLzapHud1", "vs_HLallBul1", "vs_HLcriGam1", "vs_HLsigHis1",
                   "vs_HLonyTor1", "vs_HLperNas1", "vs_HLperCri1", "vs_HLperEre1",
                   "vs_HLperCal2", "vs_HLperLeu1", "vs_HLperPol1", "vs_HLondZib1",
                   "vs_HLmyoGla2", "vs_HLellLut1", "vs_HLmicArv1", "vs_HLpsaObe1",
                   "vs_HLmerUng1", "vs_HLacoRus1", "vs_HLacoCah1", "vs_HLratRat7",
                   "vs_HLapoSyl1", "vs_HLmasCou1", "vs_HLmusPah1", "vs_HLmusCar1",
                   "vs_HLmusSpr1", "vs_HLsolPar1", "vs_eriEur2", "vs_sorAra2", "vs_HLuroGra1",
                   "vs_HLscaAqu1", "vs_HLtalOcc1", "vs_HLpteGig1", "vs_pteAle1",
                   "vs_HLcynBra1", "vs_HLeidDup1", "vs_HLeidHel2", "vs_HLeonSpe1",
                   "vs_HLrouMad1", "vs_HLcraTho1", "vs_HLmegLyr2", "vs_HLhipGal1",
                   "vs_HLhipArm1", "vs_HLrhiFer5", "vs_HLrhiSin1", "vs_HLnocLep1",
                   "vs_ptePar1", "vs_HLmorBla1", "vs_HLmacCal1", "vs_HLdesRot2",
                   "vs_HLtonSau1", "vs_HLphyDis3", "vs_HLanoCau1", "vs_HLlepYer1",
                   "vs_HLcarPer3", "vs_HLmolMol2", "vs_HLtadBra1", "vs_HLminNat1",
                   "vs_HLmurAurFea1", "vs_HLnycHum2", "vs_HLlasBor1", "vs_eptFus1",
                   "vs_HLpipPip1", "vs_myoDav1", "vs_HLmyoSep1", "vs_HLmyoLuc1",
                   "vs_HLmyoMyo6", "vs_HLmanPen2", "vs_HLparHer1", "vs_HLhyaHya1",
                   "vs_HLcryFer2", "vs_HLhelPar1", "vs_HLmunMug1", "vs_HLaciJub2",
                   "vs_HLpumYag1", "vs_HLpriBen1", "vs_HLlycPic2", "vs_HLailMel2",
                   "vs_HLursThi1", "vs_HLursArc1", "vs_HLarcGaz2", "vs_lepWed1",
                   "vs_HLmirLeo1", "vs_HLeriBar1", "vs_HLphoVit1", "vs_HLhalGry1",
                   "vs_HLspiGra1", "vs_HLailFul2", "vs_HLpotFla1", "vs_HLnasNar1",
                   "vs_HLproLot1", "vs_HLtaxTax1", "vs_HLmelCap1", "vs_HLmarZib1",
                   "vs_HLmusErm1", "vs_HLmusPut1", "vs_HLpteBra1", "vs_HLlonCan1",
                   "vs_HLlutLut1", "vs_equCab3", "vs_HLdicSum1", "vs_HLcerSimCot1",
                   "vs_HLlamGla1", "vs_susScr11", "vs_HLeubGla1", "vs_balAcu1",
                   "vs_HLbalBon1", "vs_HLbalMus1", "vs_HLbalEde1", "vs_HLescRob1",
                   "vs_HLkogBre1", "vs_HLplaMin1", "vs_HLzipCav1", "vs_lipVex1",
                   "vs_HLponBla1", "vs_HLlniGeo1", "vs_orcOrc1", "vs_HLlagObl1",
                   "vs_HLsouChi1", "vs_HLturAdu1", "vs_HLturTru3", "vs_HLtraKan1",
                   "vs_HLokaJoh2", "vs_HLgirCam1", "vs_HLhydIne1", "vs_HLcapPyg1",
                   "vs_HLalcAlc1", "vs_HLmunMun1", "vs_HLmunCri1", "vs_HLelaDav1",
                   "vs_HLmosBer1", "vs_HLtraScr1", "vs_HLbubBub2", "vs_bisBis1",
                   "vs_HLaepMel1", "vs_HLoreOre1", "vs_HLnanGra1", "vs_HLlitWal1",
                   "vs_HLkobLecLec1")





all.equal(vecFromCluster, speciesFilter)
all.equal(alternate, names(phenVec))

tree$tip.label %in% alternate

#---------------------------------------------------------------------
# --- testing diet overlap and predict  --- 
# --------------------------------------------------------------------

mergedData = read.csv("Data/mergedData.csv")

mergeDataInPantheria = which(!is.na(mergedData$panTheriaTrophicLevelCharacter))
mergeDataInWalker = which(!is.na(mergedData$MeyerTrophicLevel))
mergeDataInWalker = which(!is.na(mergedData$Meyer.Lab.Classification.Compressed))
mergeDataInElton = which(!is.na(mergedData$SimplifiedDietConvert))

mergeDataInOtherDiet = union(mergeDataInPantheria, mergeDataInWalker)
mergeDataInBothDiet = intersect(mergeDataInPantheria, mergeDataInWalker)

mergeDataInAllDiet = intersect(mergeDataInBothDiet, mergeDataInElton)

convertToTC3 = function(vector){
  vector = gsub("Carnivore", "C", vector)
  vector = gsub("Piscivore", "C", vector)
  vector = gsub("Planktivore", "C", vector)
  vector = gsub("Insectivore", "C", vector)
  vector = gsub("Hematophagy", "C", vector)
  vector = gsub("Vertivore", "C", vector)
  vector = gsub("Invertivore", "C", vector)
  
  
  vector = gsub("Omnivore", "O", vector)
  vector = gsub("Ambigous", "O", vector)
  vector = gsub("Ambiguous", "O", vector)
  vector = gsub("zMixedPredator", "O", vector)
  
  vector = gsub("Herbivore", "H", vector)
  
  #vector = gsub(NA, "N", vector)
  
}

convertToTC4 = function(vector){
  vector = gsub("Carnivore", "V", vector)
  vector = gsub("Piscivore", "V", vector)
  vector = gsub("Planktivore", "I", vector)
  vector = gsub("Insectivore", "I", vector)
  vector = gsub("Hematophagy", "C", vector)
  vector = gsub("Vertivore", "V", vector)
  vector = gsub("Invertivore", "I", vector)
  
  
  vector = gsub("Omnivore", "O", vector)
  vector = gsub("Ambigous", "O", vector)
  vector = gsub("Ambiguous", "O", vector)
  vector = gsub("zMixedPredator", NA, vector)
  
  vector = gsub("Herbivore", "H", vector)
  
  #vector = gsub(NA, "N", vector)
  
}



allDietOverlap = mergedData[c(1, 2, 24, 18, 83, 82, 79, 78)]


allDietOverlap[,4] = gsub("_", "", allDietOverlap[,4])
allDietOverlap$WalkerTC = allDietOverlap[,4]
allDietOverlap$WalkerTC = convertToTC(allDietOverlap$WalkerTC)

allDietOverlap$PantheriaTC = allDietOverlap[,3]
allDietOverlap$PantheriaTC = convertToTC(allDietOverlap$PantheriaTC)



allDietOverlap$Derek90TC = allDietOverlap$TranslatedDerek
allDietOverlap$Derek90TC = convertToTC(allDietOverlap$Derek90TC)

allDietOverlap$PiscCarnTC = allDietOverlap$SimplifiedDietConvert
allDietOverlap$PiscCarnTC = convertToTC(allDietOverlap$PiscCarnTC)

allDietOverlap$PiscOmniTC = allDietOverlap$SimplifiedDietConvertEqualDrop
allDietOverlap$PiscOmniTC = convertToTC(allDietOverlap$PiscOmniTC)

allDietOverlap$StrictPredTC = allDietOverlap$SimplifiedDietConvertStrictPred
allDietOverlap$StrictPredTC = convertToTC(allDietOverlap$StrictPredTC)


allDietOverlap$WalkerPantheria = (allDietOverlap$WalkerTC == allDietOverlap$PantheriaTC)
allDietOverlap$WalkerDerek = (allDietOverlap$WalkerTC == allDietOverlap$Derek90TC)
allDietOverlap$PantheriaDerek = (allDietOverlap$PantheriaTC == allDietOverlap$Derek90TC)

allDietOverlap$WalkerPiscCarn = (allDietOverlap$WalkerTC == allDietOverlap$PiscCarnTC)
allDietOverlap$PantheriaPiscCarn = (allDietOverlap$PantheriaTC == allDietOverlap$PiscCarnTC)

allDietOverlap$WalkerPiscOmni = (allDietOverlap$WalkerTC == allDietOverlap$PiscOmniTC)
allDietOverlap$PantheriaPiscOmni = (allDietOverlap$PantheriaTC == allDietOverlap$PiscOmniTC)

allDietOverlap$WalkerStrictPred = (allDietOverlap$WalkerTC == allDietOverlap$StrictPredTC)
allDietOverlap$PantheriaStrictPred = (allDietOverlap$PantheriaTC == allDietOverlap$StrictPredTC)


write.csv(allDietOverlap, "Results/OverlapData3.csv")




table(mergedData$Meyer.Diet.Categorical)
table(mergedData$Meyer.Lab.Classification.Clean)
table(mergedData$Meyer.Lab.Classification.Compressed)


dietCompares = mergedData[mergeDataInWalker,]

which(dietCompares$SimplifiedDietConvertEqualDrop == "zMixedPredator")

dietCompares[which(dietCompares$SimplifiedDietConvertEqualDrop == "zMixedPredator"),c(1,2,5,16,18,24, 82)]




#---------------------------------------------------------------------
# --- testing enrichment num genesets  --- 
# --------------------------------------------------------------------


enrichmentResult = fastwilcoxGMTall(rerStats, annotationsList, outputGeneVals = T, num.g =2) #run enrichment analysis 

vals = rerStats
gmt = annotationsList[[1]]
outputGeneVals = T
num.g =2

genes=NULL
simple=T 
use.all=F
order=F
alternative = "two.sided"

fastwilcoxGMT=function(vals, gmt, simple=T, use.all=F, num.g=10,genes=NULL, outputGeneVals=F, order=F,
                       alternative = "two.sided"){
  vals=vals[!is.na(vals)]
  if(is.null(genes)){
    genes=unique(unlist(gmt$genesets))
  }
  out=matrix(nrow=length(gmt$genesets), ncol=5)
  rownames(out)=gmt$geneset.names
  colnames(out)=c("stat", "pval", "p.adj","num.genes", "gene.vals")
  out=as.data.frame(out)
  genes=intersect(genes, names(vals))
  
  valsr=rank(vals[genes])
  numg=length(vals)+1
  valsallr=rank(vals)
  for( i in 1:nrow(out)){
    
    curgenes=intersect(genes,gmt$genesets[[i]])
    
    bkgenes=setdiff(genes, curgenes)
    
    if (length(bkgenes)==0 || use.all){
      bkgenes=setdiff(names(vals), curgenes)
    }
    if(length(curgenes)>=num.g & length(bkgenes)>2){
      if(!simple){
        # change alternative = "greater" for the one-sided test
        res=wilcox.test(x = vals[curgenes], y=vals[bkgenes], exact=F, alternative = alternative)
        
        out[i, 1:2]=c(res$statistic/(as.numeric(length(bkgenes))*as.numeric(length(curgenes))), res$p.value)
      }
      else{
        # add an alternative parameter (can be "greater" or "two.sided")
        out[i, 1:2]=simpleAUCgenesRanks(valsr[curgenes],valsr[bkgenes], alt = alternative)
        
      }
      out[i,"num.genes"]=length(curgenes)
      if(outputGeneVals){
        if (out[i,1]>0.5){
          oo=order(vals[curgenes], decreasing = T)
          granks=numg-valsallr[curgenes]
        }
        else{
          oo=order(vals[curgenes], decreasing = F)
          granks=valsallr[curgenes]
        }
        
        
        nn=paste(curgenes[oo],round((granks[curgenes])[oo],2),sep=':' )
        out[i,"gene.vals"]=paste(nn, collapse = ", ")
      }
    }
    
  }
  # hist(out[,2])
  out[,1]=out[,1]-0.5
  out[, "p.adj"]=p.adjust(out[,2], method="BH")
  
  out=out[!is.na(out[,2]),]
  if(order){
    out=out[order(-abs(out[,1])),]
  }
  out
}




#---------------------------------------------------------------------
# --- write species filters --- 
# --------------------------------------------------------------------




{  # Bracket used for collapsing purposes
  #File Prefix
  if(!is.na(cmdArgImport('r'))){
    filePrefix = cmdArgImport('r')
  }else{
    stop("THIS IS AN ISSUE MESSAGE; SPECIFY FILE PREFIX")
  }
  
  #  Output Directory 
  if(!dir.exists("Output")){                                      #Make output directory if it does not exist
    dir.create("Output")
  }
  outputFolderNameNoSlash = paste("Output/",filePrefix, sep = "") #Set the prefix sub directory
  if(!dir.exists(outputFolderNameNoSlash)){                       #create that directory if it does not exist
    dir.create(outputFolderNameNoSlash)
  }
  outputFolderName = paste("Output/",filePrefix,"/", sep = "")
  
  #  Force update argument
  forceUpdate = FALSE
  if(!is.na(cmdArgImport('v'))){                                 #Import if update being forced with argument 
    forceUpdate = cmdArgImport('v')
    forceUpdate = as.logical(forceUpdate)
  }else{
    message("Force update not specified, not forcing update")
  }
}

alternateSets = readRDS(paste0(outputFolderName, filePrefix, "AlternatePruningSpecies.rds"))

for(i in 1:length(alternateSets)){
  message(i)
  alternateFilePrefix = paste0("/Alternates/Alternate", i)
  currentSet = alternateSets[[i]]
  currentSpeciesFilterFilename =  paste(outputFolderName, alternateFilePrefix, filePrefix, "SpeciesFilter.rds",sep="") #set a filename for the species filter based on the prefix 
  saveRDS(currentSet, currentSpeciesFilterFilename)
}

#---------------------------------------------------------------------
# --- check for alternate mismatch --- 
# --------------------------------------------------------------------

files <- list.files("Output/ComplexDietCentralAnalysisSimplify2/Alternates", full.names = TRUE)

filters = files[grep("SpeciesFilter.rds", files)]
#load in all of the alternates 
filters = filters[order(filters)]

filters <- filters[order(as.numeric(sub(".*Alternate([0-9]+).*", "\\1", filters)))]

filtersSets = list()
for(i in 1:length(filters)){
  currentSet = readRDS(filters[i])
  filtersSets[[i]] = currentSet
}

phenVecs = files[grep("PhenotypeVector.csv", files)]
#load in all of the alternates 
phenVecs <- phenVecs[order(as.numeric(sub(".*Alternate([0-9]+).*", "\\1", phenVecs)))]

PhenSets = list()
for(i in 1:length(phenVecs)){
  currentSet = read.csv(phenVecs[i])
  PhenSets[[i]] = currentSet[,2]
  #PhenSets[[i]] = PhenSets[[i,2]]
}



currentPhenotypeVectorSavingManual = currentPhenotypeVectorSaving
manualVector = currentPhenotypeVectorSaving[,1]

PhenSets[[i]] %in% manualVector

for(i in 1:100){
  message(i)
  print(paste("A-F ", all(alternateSets[[i]] %in% filtersSets[[i]])))
  print(paste("A-P ", all(alternateSets[[i]] %in% PhenSets[[i]])))
  print(paste("F-A ", all(filtersSets[[i]] %in% alternateSets[[i]])))
  print(paste("F-P ", all(filtersSets[[i]] %in% PhenSets[[i]])))
  print(paste("P-A ", all(PhenSets[[i]] %in% alternateSets[[i]])))
        
  
  
}

TreeSets = list()
altTrees = files[grep("CategoricalTree.rds", files)]
for(i in 1:length(altTrees)){
  currentSet = readRDS(altTrees[i])
  TreeSets[[i]] = currentSet
}

masterTree = mainTrees$masterTree
shortsNumbers = vector()
for(i in 1:length(TreeSets)){
  currentTree = TreeSets[[i]]
  currentSp = currentTree$tip.label
  tipsToRemove = masterTree$tip.label[which(!masterTree$tip.label %in% currentSp)]
  
  currentLengthTree = drop.tip(masterTree, tipsToRemove)
  
  
  currentTerminals = which(currentLengthTree$edge[,2] <= length(currentLengthTree$tip.label))
  
  numShorts = length(which(currentLengthTree$edge.length[currentTerminals] < 0.01))
  shortsNumbers[i] = numShorts
}
max(shortsNumbers)

for(i in 1:length(alternateSets)){
  currentTree = TreeSets[[i]]
  currentSp = alternateSets[[i]]
  tipsToRemove = masterTree$tip.label[which(!masterTree$tip.label %in% currentSp)]
  
  currentLengthTree = drop.tip(masterTree, tipsToRemove)
  
  
  currentTerminals = which(currentLengthTree$edge[,2] <= length(currentLengthTree$tip.label))
  
  numShorts = length(which(currentLengthTree$edge.length[currentTerminals] < 0.01))
  shortsNumbers[i] = numShorts
}
#---------------------------------------------------------------------
# --- nebug alternates --- 
# --------------------------------------------------------------------
brokenSet = currentSet
brokenFilter = speciesFilter
brokenTree = phenotypeTree

"Output/ComplexDietCentralAnalysisSimplify/Alternates/Alternate20ComplexDietCentralAnalysisSimplifySpeciesFilter.rds"



removalTips = data.frame(startRemoveTips1, startRemoveTips2, startRemoveTips3)

startRemoveTips1 %in% startRemoveTips2
startRemoveTips1 %in% startRemoveTips3
startRemoveTips = startRemoveTips1

masterShortTipNames = masterTree$tip.label[masterShortTips]


shortPhens = fullDataPhenotype[names(fullDataPhenotype) %in% masterShortTipNames]
table(shortPhens)
removePhens = fullDataPhenotype[names(fullDataPhenotype) %in% startRemoveTips] 
table(removePhens)


length(which(speciesSets[[j]] %in% startRemoveTips))


phenCHeck = fullDataPhenotype[names(fullDataPhenotype) %in% speciesSets[[j]]] 
table(phenCHeck)

allAlternates = unlist(alternateSets)
length(table(allAlternates))



aoverlapsLengths = integer()
alengths = integer()
for(i in 1:100){
  length = length(alternateSets[[i]])
  alengths = append(alengths, length)
  for(j in 2:100){
    if(i != j){
      overlap = length(which(alternateSets[[i]] %in% alternateSets[[j]]))
      aoverlapsLengths = append(aoverlapsLengths, overlap)
    }
    
  }
}

table(alengths)
table(aoverlapsLengths)




files <- list.files("Output/ComplexDietCentralAnalysisSimplify/Alternates", full.names = TRUE)
files <- list.files("Output/ComplexDietCentralAnalysis/Alternates", full.names = TRUE)


filters = files[grep("SpeciesFilter.rds", files)]


#load in all of the alternates 

filtersSets = list()
for(i in 1:length(filters)){
  currentSet = readRDS(filters[i])
  filtersSets[[i]] = currentSet
}

foverlapsLengths = integer()
flengths = integer()
for(i in 1:100){
  length = length(filtersSets[[i]])
  flengths = append(flengths, length)
  for(j in 2:100){
    if(i != j){
      overlap = length(which(filtersSets[[i]] %in% filtersSets[[j]]))
      foverlapsLengths = append(foverlapsLengths, overlap)
    }
    
  }
}

table(flengths)
table(foverlapsLengths)

files <- list.files("Output/ComplexDietCentralAnalysisSimplify/Alternates", full.names = TRUE)
files <- list.files("Output/ComplexDietCentralAnalysis/Alternates", full.names = TRUE)


phenVecs = files[grep("PhenotypeVector.rds", files)]


#load in all of the alternates 

PhenSets = list()
for(i in 1:length(phenVecs)){
  currentSet = readRDS(phenVecs[i])
  PhenSets[[i]] = currentSet
}


poverlapsLengths = integer()
plengths = integer()
for(i in 1:100){
  length = length(PhenSets[[i]])
  plengths = append(plengths, length)
  for(j in 2:100){
    if(i != j){
      overlap = length(which(names(PhenSets[[i]]) %in% names(PhenSets[[j]])))
      poverlapsLengths = append(poverlapsLengths, overlap)
    }
    
  }
}

table(plengths)
table(poverlapsLengths)



table(alengths)
table(flengths)
table(plengths)

table(aoverlapsLengths)
table(foverlapsLengths)
table(poverlapsLengths)

filtersSets[[i]] %in% names(PhenSets[[i]])


phenotypeVectorFilename = paste(outputFolderName, filePrefix, "CategoricalPhenotypeVector.rds",sep="")
phenotypeVector = readRDS(phenotypeVectorFilename)


speciesFilterFilename = paste(outputFolderName, filePrefix, "SpeciesFilter.rds",sep="") #set a filename for the species filter based on the prefix 

relevantSpecieslist = readRDS(speciesFilterFilename)                          #if not, use the existing list 
speciesFilter = relevantSpecieslist                                           #make the speciesFilter object for later 
relevantSpecies = manualAnnots[ manualAnnots[[nameColumn]] %in% relevantSpecieslist,] #and select the manual annotations entries in that list (useful if the list is more restrictive than it would be by default) 
irrelevantSpecies = manualAnnots[! manualAnnots[[nameColumn]] %in% relevantSpecieslist,]

all.equal(speciesFilter, relevantSpecies)
length(speciesFilter)


phenVec = readRDS("")


currentSet %in% speciesFilter
speciesFilter %in%  currentSet

phenotypeTree$tip.label %in% speciesFilter




testTree = mainTrees$masterTree
tipsToDrop = testTree$tip.label[!testTree$tip.label %in% alternateSets[[1]]]
testTree = drop.tip(testTree, tipsToDrop)
i=i+1
if(i %% 10000 == 0){message(i)}
if(min(testTree$edge.length) < pruningCutoff){
  message("Found Valid Alternate")
  message(i)
  alternateTips = testTree$tip.label
  alternateTips = list(alternateTips)
  alternateSets = append(alternateSets, alternateTips)
}

prunedTree$edge.length

numTips = length(testTree$tip.label)

edgeBranches = which(testTree$edge[,2] <= numTips)
edgeBranches2 = which(!testTree$edge[,2] %in% testTree$edge[,1])


all.equal(edgeBranches, edgeBranches2)

edgeBranches2[which(!edgeBranches2 %in% edgeBranches)]

testTree$edge[905,]

min(testTree$edge[,1])
length(testTree$edge)



#---------------------------------------------------------------------
# --- Summarize binary analyses --- 
# --------------------------------------------------------------------




# --- binary results 


BovCor = readRDS("Output/CladeBinaryBovidae/0-1/CladeBinaryBovidae0-1CorrelationFile.rds")
CerCor = readRDS("Output/CladeBinaryCervidae/0-1/CladeBinaryCervidae0-1CorrelationFile.rds")
CriCor = readRDS("Output/CladeBinaryCricetidae/0-1/CladeBinaryCricetidae0-1CorrelationFile.rds")
HysCor = readRDS("Output/CladeBinaryHystricognathi/0-1/CladeBinaryHystricognathi0-1CorrelationFile.rds")
PerCor = readRDS("Output/CladeBinaryPeropdidae/0-1/CladeBinaryPeropdidae0-1CorrelationFile.rds")
VesCor = readRDS("Output/CladeBinaryVespertilionidae/0-1/CladeBinaryVespertilionidae0-1CorrelationFile.rds")


#---------------------------------------------------------------------
# --- Compare diet categorizations  --- 
# --------------------------------------------------------------------

mergedData = read.csv("Data/MergedData.csv")
mainSpeciesSet = readRDS("OUtput/ComplexDietCentralAnalysis/ComplexDietCentralAnalysisSpeciesFilter.rds")

mergedData$inMainAnalysis = mergedData$ZoonomiaTip %in% mainSpeciesSet

dietConversions = mergedData[,c(2,22,26,29,37,78,79,80)]



dietConversions[which(dietConversions$TranslatedDerek != dietConversions$SimplifiedDietConvertNoScav & mergedData$inMainAnalysis),]

dietConversions[which(dietConversions$TranslatedDerek != dietConversions$SimplifiedDietConvert & mergedData$inMainAnalysis),]

dietConversions[which( mergedData$SimplifiedDietConvertNoScav == "MixedPredator" & mergedData$Diet.Vfish >= 50),]



mainPhenv = readRDS("OUtput/ComplexDietCentralAnalysis/ComplexDietCentralAnalysisCategoricalPhenotypeVector.rds")
simpPhenv = readRDS("OUtput/ComplexDietCentralAnalysisSimplify/ComplexDietCentralAnalysisSimplifyCategoricalPhenotypeVector.rds")

names(mainPhenv) == names(simpPhenv)

phenvCompare = data.table(names(mainPhenv), mainPhenv, simpPhenv)
phenvCompare$mainPhenv = gsub("Insectivore", "Invertivore", phenvCompare$mainPhenv)


phenvCompare[which(phenvCompare$mainPhenv != phenvCompare$simpPhenv),]




#---------------------------------------------------------------------
# --- Compare new diet assignment  --- 
# --------------------------------------------------------------------


library(tibble)
library(ggvenn)
library(gridExtra)


convertColnameToIndex = function(name, df){
  output = which(colnames(df)==name)
  output
}
makeVenn = function(i, titleName,  gene = T){
  if(gene){
  venn_list <- list(
    Main = rownames(mainSigDataGene)[which(mainSigDataGene[,i])],
    Simplified = rownames(simpSigDataGene)[which(simpSigDataGene[,i])]
  )
  titleName = paste(titleName, "Gene")
  }else{
    venn_list <- list(
      Main = rownames(mainSigDataGO)[which(mainSigDataGO[,i])],
      Simplified = rownames(simpSigDataGO)[which(simpSigDataGO[,i])]
    )
    titleName = paste(titleName, "GO")
  }
  vennout = ggvenn(venn_list) +  labs(title = titleName)
  vennout
}
makeComboVenn = function( gene = T){
  venn1 = makeVenn(1, "CH", gene)
  venn2 = makeVenn(2, "HI", gene)
  venn3 = makeVenn(3, "HV", gene)
  venn4 = makeVenn(4, "IV", gene)
  
  vennSet = grid.arrange(venn1, venn2, venn3, venn4)
  vennSet
}

mainDataGene = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysiscombinedGeneResultsWithAlternates.rds")
simpDataGene = readRDS("Output/ComplexDietCentralAnalysisSimplify/ComplexDietCentralAnalysisSimplifycombinedGeneResults.rds")
simpDataGene = readRDS("Output/ComplexDietCentralAnalysisSimplifyNoScav/ComplexDietCentralAnalysisSimplifyNoScavcombinedGeneResults.rds")
simpDataGene = readRDS("Output/ComplexDietCentralAnalysisPiscFix/ComplexDietCentralAnalysisPiscFixcombinedGeneResults.rds")
simpDataGene = readRDS("Output/ComplexDietCentralAnalysisNoFishbat/ComplexDietCentralAnalysisNoFishbatcombinedGeneResults.rds")
simpDataGene = readRDS("Output/ComplexDietCentralAnalysisSimplifyEqualDrop/ComplexDietCentralAnalysisSimplifyEqualDropcombinedGeneResultsWithAlternates.rds")
simpDataGene = readRDS("Output/ComplexDietCentralAnalysisSimplify2/ComplexDietCentralAnalysisSimplify2combinedGeneResultsWithAlternates.rds")
simpDataGene = readRDS("Output/ComplexDietCentralAnalysisSimplifyExclusivePred/ComplexDietCentralAnalysisSimplifyExclusivePredcombinedGeneResultsWithAlternates.rds")

mainDataGene = readRDS("Output/ComplexDietCentralAnalysisSimplify2/ComplexDietCentralAnalysisSimplify2combinedGeneResultsWithAlternates.rds")



mainDataGO = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysiscombinedGOResults-KeggReactomeOld.rds")
mainDataGO = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysiscombinedGOResultsWithAlternates-KeggReactome.rds")
simpDataGO = readRDS("Output/ComplexDietCentralAnalysisSimplify/ComplexDietCentralAnalysisSimplifycombinedGOResults-KeggReactome.rds")
simpDataGO = readRDS("Output/ComplexDietCentralAnalysisSimplifyNoScav/ComplexDietCentralAnalysisSimplifyNoScavcombinedGOResults-KeggReactome.rds")
simpDataGO = readRDS("Output/ComplexDietCentralAnalysisPiscFix/ComplexDietCentralAnalysisPiscFixcombinedGOResults-KeggReactome.rds")
simpDataGO = readRDS("Output/ComplexDietCentralAnalysisNoFishbat/ComplexDietCentralAnalysisNoFishbatcombinedGOResults-KeggReactome.rds")
simpDataGO = readRDS("Output/ComplexDietCentralAnalysisSimplify2/ComplexDietCentralAnalysisSimplify2combinedGOResultsWithAlternates-KeggReactome.rds")
simpDataGO = readRDS("Output/ComplexDietCentralAnalysisSimplifyExclusivePred/ComplexDietCentralAnalysisSimplifyExclusivePredcombinedGOResultsWithAlternates-KeggReactome.rds")

simpDataGO = read.csv("Output/ComplexDietCentralAnalysisSimplifyExclusivePred/ComplexDietCentralAnalysisSimplifyExclusivePredcombinedGOResultsWithAlternates-KeggReactomeOLD.csv")
rownames(simpDataGO) = simpDataGO$X
mainDataGO = readRDS("Output/ComplexDietCentralAnalysisSimplify2/ComplexDietCentralAnalysisSimplify2combinedGOResultsWithAlternates-KeggReactome.rds")



#simpData = simpData[,c(1:42)]

mainSigDataGene = mainDataGene[,c(convertColnameToIndex("CH-significant", mainDataGene), convertColnameToIndex("HI-significant", mainDataGene), convertColnameToIndex("HV-significant", mainDataGene), convertColnameToIndex("IV-significant", mainDataGene))]
mainSigDataGene = mainDataGene[,c(convertColnameToIndex("CH-significantRobust", mainDataGene), convertColnameToIndex("HI-significantRobust", mainDataGene), convertColnameToIndex("HV-significantRobust", mainDataGene), convertColnameToIndex("IV-significantRobust", mainDataGene))]

mainSigDataGO = mainDataGO[,c(convertColnameToIndex("CH-significant", mainDataGO), convertColnameToIndex("HI-significant", mainDataGO), convertColnameToIndex("HV-significant", mainDataGO), convertColnameToIndex("IV-significant", mainDataGO))]
mainSigDataGO = mainDataGO[,c(convertColnameToIndex("CH.significantRobust", mainDataGO), convertColnameToIndex("HI.significantRobust", mainDataGO), convertColnameToIndex("HV.significantRobust", mainDataGO), convertColnameToIndex("IV.significantRobust", mainDataGO))]

mainSigDataGO = mainDataGO[,c(convertColnameToIndex("CH-significantRobust", mainDataGO), convertColnameToIndex("HI-significantRobust", mainDataGO), convertColnameToIndex("HV-significantRobust", mainDataGO), convertColnameToIndex("IV-significantRobust", mainDataGO))]


simpSigDataGene = simpDataGene[,c(convertColnameToIndex("CH-significant", simpDataGene), convertColnameToIndex("HI-significant", simpDataGene), convertColnameToIndex("HV-significant", simpDataGene), convertColnameToIndex("IV-significant", simpDataGene))]
simpSigDataGene = simpDataGene[,c(convertColnameToIndex("CH-significantRobust", simpDataGene), convertColnameToIndex("HI-significantRobust", simpDataGene), convertColnameToIndex("HV-significantRobust", simpDataGene), convertColnameToIndex("IV-significantRobust", simpDataGene))]

simpSigDataGO = simpDataGO[,c(convertColnameToIndex("CH-significant", simpDataGO), convertColnameToIndex("HI-significant", simpDataGO), convertColnameToIndex("HV-significant", simpDataGO), convertColnameToIndex("IV-significant", simpDataGO))]
simpSigDataGO = simpDataGO[,c(convertColnameToIndex("CH.significantRobust", simpDataGO), convertColnameToIndex("HI.significantRobust", simpDataGO), convertColnameToIndex("HV.significantRobust", simpDataGO), convertColnameToIndex("IV.significantRobust", simpDataGO))]
simpSigDataGO = simpDataGO[,c(convertColnameToIndex("CH-significantRobust", simpDataGO), convertColnameToIndex("HI-significantRobust", simpDataGO), convertColnameToIndex("HV-significantRobust", simpDataGO), convertColnameToIndex("IV-significantRobust", simpDataGO))]





geneVenn = makeComboVenn(T)
goVenn = makeComboVenn(F)
# megaVenn = grid.arrange(geneVenn, goVenn, ncol =1)
megaVenn



simpIVSig = rownames(simpDataGO)[which(simpDataGO$`IV-significantRobust`)]
mainIVSig = rownames(mainDataGO)[which(mainDataGO$`IV-significantRobust`)]
IvSig = append(simpIVSig, mainIVSig)

IVOUtsMain = mainDataGO[rownames(mainDataGO) %in% IvSig,grep("IV", colnames(mainDataGO))]
IVOUtsSimp = simpDataGO[rownames(simpDataGO) %in% IvSig,grep("IV", colnames(simpDataGO))]

names(IVOUtsSimp) = paste0("simp-", names(IVOUtsSimp))

IVOutCombined = cbind(IVOUtsMain, IVOUtsSimp)


IVsigGO = rownames(simpDataGO)[simpDataGO$`IV-significantRobust`]
IVMainGO = rownames(mainDataGO)[mainDataGO$`IV-significantRobust`]

IvOverlap = IVsigGO[IVsigGO %in% IVMainGO]

CHsigGO = rownames(simpDataGO)[simpDataGO$`CH-significantRobust`]
CHMainGO = rownames(mainDataGO)[mainDataGO$`CH-significantRobust`]

CHOverlap = CHsigGO[CHsigGO %in% CHMainGO]


HIsigGO = rownames(simpDataGO)[simpDataGO$`HI-significantRobust`]
HIMainGO = rownames(mainDataGO)[mainDataGO$`HI-significantRobust`]

HIOverlap = HIsigGO[HIsigGO %in% HIMainGO]



IvOLMain = mainDataGO[IvOverlap,grep("IV", colnames(mainDataGO))]
IvOLSimp = simpDataGO[IvOverlap,grep("IV", colnames(simpDataGO))]

names(IvOLSimp) = paste0("simp-", names(IvOLSimp))

IvOlCombined = cbind(IvOLMain, IvOLSimp)







table(simpDataGO$`IV-PadjNumSignificant`)
#there are 13 of them significant in more than 

which(simpDataGO$`IV-PadjNumSignificant` > 51)
grep("IV", colnames(simpDataGO))



convertColnameToIndex("IV-significant", simpData)
plot(mainData$`IV-p.adj`, simpData$`IV-p.adj`)

which(mainSigData$`CH-significantRobust` & !simpSigData$`CH-significant`)

IvMissingInSimp = which(mainSigData$`IV-significantRobust` & !simpSigData$`IV-significant`)
IvMissingInSimp = which(mainSigData$`IV-significant` & !simpSigData$`IV-significant`)

IvMissingInMain = which(!mainSigData$`IV-significantRobust` & simpSigData$`IV-significant`)

length(which(mainSigData$`IV-significantRobust` & !simpSigData$`IV-significant`))

simpData[IvMissingInSimp,]


simpData[which(mainSigData$`HI-significantRobust` & !simpSigData$`HI-significant`),]

IVData = cbind(mainData[,c(37:42)], simpData[,c(37:42)])

IVData[IvMissingInSimp, c(1,2)]

missingSets = IVData[c(545, 546),]


missingSets[1,c(5,11)]


get_positions <- function(x) {
  pairs <- trimws(unlist(strsplit(x, ",")))
  
  genes <- sub(":.*", "", pairs)
  positions <- as.numeric(sub(".*:", "", pairs))
  
  setNames(positions, genes)
}

sumImportantChange
IVData$sumChange = rep(NA)
IVData$sumImportantChange = rep(NA)

i=33
for(i in 1:1803){
IV_gene_vals = IVData[i,c(5)]
IV_gene_vals_1 = IVData[i,c(11)]
  


pos1 <- get_positions(IV_gene_vals)
pos2 <- get_positions(IV_gene_vals_1)

# Genes present in both lists
common_genes <- intersect(names(pos1), names(pos2))

# Position changes
position_changes <- data.frame(
  Gene = common_genes,
  MainPosition = pos1[common_genes],
  SimpPosition = pos2[common_genes],
  Change = pos2[common_genes] - pos1[common_genes]
)

# Only genes whose position actually changed
position_changes <- position_changes[
  position_changes$MainPosition != position_changes$SimpPosition,
]

position_changes
numImportantPositions = max(length(which(position_changes$MainPosition < 1000)), length(which(position_changes$SimpPosition < 1000)))
IVData$sumChange[i] = sum(position_changes$Change)
IVData$sumImportantChange[i] = sum(head(position_changes$Change,numImportantPositions))

}

IVData$pValChange = (IVData[,3]) - (IVData[,9])

plot(IVData$sumChange, IVData$pValChange)
abline(v = 0)  # vertical line at x = 0
abline(h = 0)  # horizontal line at y = 0

plot(IVData$sumImportantChange, IVData$pValChange, xlim = c(min(IVData$sumImportantChange), -min(IVData$sumImportantChange)))



missingSets = IVData[IvMissingInSimp,]


length(which())

confusingCol = mainSigData$`IV-significant`

all.equal(confusingCol,mainSigData[,4])

#---------------------------------------------------------------------
# --- Debug Alternates --- 
# --------------------------------------------------------------------

files <- list.files("Output/ComplexDietCentralAnalysisSimplify/Alternates", full.names = TRUE)

phenVecs = files[grep("PhenotypeVector.rds", files)]


#load in all of the alternates 

PhenSets = list()
for(i in 1:length(phenVecs)){
  currentSet = readRDS(phenVecs[i])
  PhenSets[[i]] = currentSet
}

phenFullVec = unlist(PhenSets)

table(phenFullVec)



#---------------------------------------------------------------------
# --- Compare new and old trees --- 
# --------------------------------------------------------------------

oldFullTree = readRDS("Output/ComplexDietCentralAnalysisAllSpecies/ComplexDietCentralAnalysisAllSpeciesCategoricalCommonTree.rds")
newFullTree = readRDS("Output/ComplexDietCentralAnalysisSimplifyAllSpecies/ComplexDietCentralAnalysisSimplifyAllSpeciesCategoricalCommonTree.rds")

oldFullTree$tip.label[which(!oldFullTree$tip.label %in% newFullTree$tip.label)]

changedBranches = which(newFullTree$edge.length != oldFullTree$edge.length)

newFullTree$edge[changedBranches,]
edgcols = rep("black", length(newFullTree$edge.length))
edgcols[changedBranches] = "hotpink"
pdf("testTree.pdf", height = 40, width = 14)
plot.phylo(newFullTree, font = 2, edge.color = edgcols)
dev.off()

#Looks all good 

oldMainTree = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysisCategoricalCommonTree.rds")
newMainTree = readRDS("Output/ComplexDietCentralAnalysisSimplify/ComplexDietCentralAnalysisSimplifyCategoricalCommonTree.rds")

all.equal(oldMainTree$tip.label, newMainTree$tip.label)

files <- list.files("Output/ComplexDietCentralAnalysisSimplify/Alternates", full.names = TRUE)

file.rename(
  files,
  file.path(
    dirname(files),
    gsub("alysisSpecies", "alysisSimplifySpecies", basename(files))
  )
)

#load in all of the alternates 

alternateSets = list()
for(i in 1:length(files)){
  currentSet = readRDS(files[i])
  alternateSets[[i]] = currentSet
}


#---------------------------------------------------------------------
# --- Futzing with classification --- 
# --------------------------------------------------------------------

mergedData = read.csv("Data/mergedData.csv")

rawColumn = mergedData$DerekDietClassification90InsVertivoreSorting



substitutions=list(
            c("C-Invertebrate-eater", "Invertivore"), c("C-InsVertivore-Insectivore", "Invertivore"),
            c("C-Herpetivore", "Vertivore"),
            c("C-Piscivore", "Vertivore"), c("C-InsVertivore-Piscivore", "Vertivore"),
            c("C-Endotherm-Carnivore", "Vertivore"), c("C-Scavenger", "Vertivore"), c("C-Nonspecific-Vertebrate-eater", "Vertivore"),
            c("C-Terrestrial-vertebrates-eater", "Vertivore"), c("C-All-vertebrate-eater", "Vertivore"), c("C-InsVertivore-Carnivore", "Vertivore"),
            c("C-InsVertivore-Mixed", "Omnivore"), 
            c("O-For Examination", "Omnivore"), c("O-Scavenger", "Omnivore"),
            c("H-Frugivore", "Herbivore"), 
            c("H-Nectarivore", "Herbivore"), 
            c("H-High-sugar-plants-Eater", "Herbivore"),
            c("H-Granivore", "Herbivore"), c("H-Nonspecific-Herbivore", "Herbivore"), 
            c("H-Low-sugar-plants-Eater", "Herbivore"), c("H-All-plants-Eater", "Herbivore"),
            c("O-Generalist", "Omnivore")
          )

replacedColumn = rawColumn 
for( i in 1:length(substitutions)){
  substitutePhenotypes = substitutions[[i]]
  message(paste("replacing", substitutePhenotypes[1], "with", substitutePhenotypes[2]))
  replacedColumn = gsub(substitutePhenotypes[1], substitutePhenotypes[2], replacedColumn)
}


mergedData$VertAllPlusScav = mergedData$Diet.Scav + mergedData$Diet.VertAll

mergedData$TranslatedDerek = replacedColumn

mergedData$SimplifiedDietConvertNoScav = rep(NA)
mergedData$SimplifiedDietConvertNoScav[which(is.na(mergedData$Diet.Inv))] = "NoData"
mergedData$SimplifiedDietConvertNoScav[which(mergedData$Diet.PlantAll >= 90)] = "Herbivore"
mergedData$SimplifiedDietConvertNoScav[which((mergedData$Diet.VertAll + mergedData$Diet.Inv) >= 90)] = "MixedPredator"
mergedData$SimplifiedDietConvertNoScav[which(mergedData$Diet.Inv >= 90)] = "Invertivore"
mergedData$SimplifiedDietConvertNoScav[which(mergedData$Diet.VertAll >= 90)] = "Vertivore"
mergedData$SimplifiedDietConvertNoScav[is.na(mergedData$SimplifiedDietConvertNoScav)] = "Omnivore"

mergedData$SimplifiedDietConvertNoScav[which( mergedData$SimplifiedDietConvertNoScav == "MixedPredator" & mergedData$Diet.Inv > 50)] = "Invertivore"
mergedData$SimplifiedDietConvertNoScav[which( mergedData$SimplifiedDietConvertNoScav == "MixedPredator" & mergedData$Diet.VertAll > 50)] = "Vertivore"

mergedData$SimplifiedDietConvertNoScav[which( mergedData$SimplifiedDietConvertNoScav == "MixedPredator" & mergedData$Diet.Vfish >= 50)] = "Vertivore"

mergedData$SimplifiedDietConvertNoScav[which( mergedData$SimplifiedDietConvertNoScav == "NoData")] = NA


mergedData$SimplifiedDietConvert = rep(NA)
mergedData$SimplifiedDietConvert[which(is.na(mergedData$Diet.Inv))] = "NoData"
mergedData$SimplifiedDietConvert[which(mergedData$Diet.PlantAll >= 90)] = "Herbivore"
mergedData$SimplifiedDietConvert[which((mergedData$VertAllPlusScav + mergedData$Diet.Inv) >= 90)] = "MixedPredator"
mergedData$SimplifiedDietConvert[which(mergedData$Diet.Inv >= 90)] = "Invertivore"
mergedData$SimplifiedDietConvert[which(mergedData$VertAllPlusScav >= 90)] = "Vertivore"
mergedData$SimplifiedDietConvert[is.na(mergedData$SimplifiedDietConvert)] = "Omnivore"

mergedData$SimplifiedDietConvert[which( mergedData$SimplifiedDietConvert == "MixedPredator" & mergedData$Diet.Inv > 50)] = "Invertivore"
mergedData$SimplifiedDietConvert[which( mergedData$SimplifiedDietConvert == "MixedPredator" & mergedData$VertAllPlusScav > 50)] = "Vertivore"

mergedData$SimplifiedDietConvert[which( mergedData$SimplifiedDietConvert == "MixedPredator" & mergedData$Diet.Vfish >= 50)] = "Vertivore"
#mergedData$SimplifiedDietConvert[which( mergedData$SimplifiedDietConvert == "MixedPredator")] = "Omnivore"


mergedData$SimplifiedDietConvert[which( mergedData$SimplifiedDietConvert == "NoData")] = NA


mergedData$SimplifiedDietConvertEqualDrop = rep(NA)
mergedData$SimplifiedDietConvertEqualDrop[which(is.na(mergedData$Diet.Inv))] = "NoData"
mergedData$SimplifiedDietConvertEqualDrop[which(mergedData$Diet.PlantAll >= 90)] = "Herbivore"
mergedData$SimplifiedDietConvertEqualDrop[which((mergedData$VertAllPlusScav + mergedData$Diet.Inv) >= 90)] = "zMixedPredator"
mergedData$SimplifiedDietConvertEqualDrop[which(mergedData$Diet.Inv >= 90)] = "Invertivore"
mergedData$SimplifiedDietConvertEqualDrop[which(mergedData$VertAllPlusScav >= 90)] = "Vertivore"
mergedData$SimplifiedDietConvertEqualDrop[is.na(mergedData$SimplifiedDietConvertEqualDrop)] = "Omnivore"

mergedData$SimplifiedDietConvertEqualDrop[which( mergedData$SimplifiedDietConvertEqualDrop == "zMixedPredator" & mergedData$Diet.Inv > 50)] = "Invertivore"
mergedData$SimplifiedDietConvertEqualDrop[which( mergedData$SimplifiedDietConvertEqualDrop == "zMixedPredator" & mergedData$VertAllPlusScav > 50)] = "Vertivore"

mergedData$SimplifiedDietConvertEqualDrop[which(mergedData$SimplifiedDietConvertEqualDrop == "NoData")] = NA



mergedData$SimplifiedDietConvertStrictPred = rep(NA)
mergedData$SimplifiedDietConvertStrictPred[which(is.na(mergedData$Diet.Inv))] = "NoData"
mergedData$SimplifiedDietConvertStrictPred[which(mergedData$Diet.PlantAll >= 90)] = "Herbivore"
mergedData$SimplifiedDietConvertStrictPred[which((mergedData$VertAllPlusScav + mergedData$Diet.Inv) >= 90)] = "zMixedPredator"
mergedData$SimplifiedDietConvertStrictPred[which(mergedData$Diet.Inv >= 90)] = "Invertivore"
mergedData$SimplifiedDietConvertStrictPred[which(mergedData$VertAllPlusScav >= 90)] = "Vertivore"
mergedData$SimplifiedDietConvertStrictPred[is.na(mergedData$SimplifiedDietConvertStrictPred)] = "Omnivore"

mergedData$SimplifiedDietConvertStrictPred[which(mergedData$SimplifiedDietConvertStrictPred == "NoData")] = NA



which( mergedData$SimplifiedDietConvert == "MixedPredator")

mergedData[which(!mergedData$SimplifiedDietConvertNoScav == mergedData$SimplifiedDietConvert),c(2,78,79,80)]
mergedData[which(!mergedData$TranslatedDerek == mergedData$SimplifiedDietConvertNoScav),c(2,78,79)]
length(which(!mergedData$TranslatedDerek == mergedData$SimplifiedDietConvertNoScav))
mergedData[15, c(2,77,78)]

#write.csv(mergedData, "Data/mergedData.csv")

#---------------------------------------------------------------------
# --- inclusion columns for supplement --- 
# --------------------------------------------------------------------

supTable = read.csv("Results/S1.csv")




mainSpecies = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysisSpeciesFilter.rds")

inMain = (supTable$Tree.Tip.Name %in% mainSpecies)


alternatesSpecies = list()


files <- list.files(
  pattern = "^Alternate[0-9]+ComplexDietCentralAnalysisCategoricalPhenotypeVector\\.rds$",
  full.names = TRUE
)

i <- sub(
  "^Alternate([0-9]+)ComplexDietCentralAnalysisCategoricalPhenotypeVector\\.rds$",
  "\\1",
  basename(files)
)

folder <- "Output/ComplexdietCentralAnalysis/Alternates"

files <- list.files(
  folder,
  pattern = "^Alternate[0-9]+ComplexDietCentralAnalysisCategoricalPhenotypeVector\\.rds$",
  full.names = TRUE
)

dataList <- lapply(files, readRDS)

names(dataList) <- sub(
  "^Alternate([0-9]+)ComplexDietCentralAnalysisCategoricalPhenotypeVector\\.rds$",
  "\\1",
  basename(files)
)

dataList <- lapply(dataList, names)
combinedVector <- unlist(dataList, use.names = FALSE)
numberofAlternates = table(combinedVector)

nameCounts <- as.data.frame(table(combinedVector))
names(nameCounts) <- c("Name", "Count")


toAdd = supTable$Tree.Tip.Name[which(!supTable$Tree.Tip.Name %in% nameCounts$Name)]
emptyCOunts = rep(0, length(toAdd))

addTable = data.frame(toAdd, emptyCOunts)
names(addTable) <- c("Name", "Count")
nameCounts = rbind (nameCounts, addTable)

nameCounts = nameCounts[order(match(nameCounts$Name, supTable$Tree.Tip.Name)),]



supTable$inMain = inMain
supTable$numAltern = nameCounts$Count

write.csv(supTable, "Results/OUtSup1.csv")
#---------------------------------------------------------------------
# --- Looking into TEs --- 
# --------------------------------------------------------------------


library(xlsx)
transposonData = read.xlsx("Data/science.abn1430_tables_s1_to_s8.xlsx", 8)

mainData = read.csv("Data/MergedData.csv")

length(which(tolower(transposonData$Species) %in% mainData$Scientific_Binomial))

tolower()

mainData$Scientific_Binomial[340]

transposonData$newDiet = rep(NA, 178)
match(tolower(transposonData$Species), mainData$Scientific_Binomial)
transposonData$newDiet = mainData$DerekDietClassification90InsVertivoreSorting[match(tolower(transposonData$Species), mainData$Scientific_Binomial)]
substitutions =list(
            c("C-Invertebrate-eater", "Insectivore"), c("C-InsVertivore-Insectivore", "Insectivore"),
            c("C-Herpetivore", "Vertivore"),
            c("C-Piscivore", "Vertivore"), c("C-InsVertivore-Piscivore", "Vertivore"),
            c("C-Endotherm-Carnivore", "Vertivore"), c("C-Scavenger", "Vertivore"), c("C-Nonspecific-Vertebrate-eater", "Vertivore"),
            c("C-Terrestrial-vertebrates-eater", "Vertivore"), c("C-All-vertebrate-eater", "Vertivore"), c("C-InsVertivore-Carnivore", "Vertivore"),
            c("C-InsVertivore-Mixed", "Omnivore"), 
            c("O-For Examination", "Omnivore"), c("O-Scavenger", "Omnivore"),
            c("H-Frugivore", "Herbivore"), 
            c("H-Nectarivore", "Herbivore"), 
            c("H-High-sugar-plants-Eater", "Herbivore"),
            c("H-Granivore", "Herbivore"), c("H-Nonspecific-Herbivore", "Herbivore"), 
            c("H-Low-sugar-plants-Eater", "Herbivore"), c("H-All-plants-Eater", "Herbivore"),
            c("O-Generalist", "Omnivore")
          )

for( i in 1:length(substitutions)){
  substitutePhenotypes = substitutions[[i]]
  message(paste("replacing", substitutePhenotypes[1], "with", substitutePhenotypes[2]))
  transposonData$newDiet = gsub(substitutePhenotypes[1], substitutePhenotypes[2], transposonData$newDiet)
}

legend = F
colorScale = (c("#1B9E77", "#000000", "#7570B3", "#E7298A"))

plot = ggplot(transposonData, aes(x=newDiet, y=log(Young_DNA_proportion), col=newDiet)) + 
  geom_violin(adjust=1/3, show.legend = legend) +
  geom_jitter(position=position_jitter(0.2), show.legend = legend) +
  theme_classic() +
  scale_color_manual(values=colorScale) +
  theme(text = element_text(size = 20))




#---------------------------------------------------------------------
# --- Comapring result from equal branch length inference--- 
# --------------------------------------------------------------------

paperGenes = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysiscombinedGeneResultsWithAlternatesPermulated.rds")

test = paperGenes[which(paperGenes$`IV-significantRobust`),]





paperAnalysis = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysisPairwiseCorrelationFile.rds")
ELAnalysis = readRDS("Output/ComplexDietCentralAnalysisEqualLengthInference/ComplexDietCentralAnalysisEqualLengthInferencePairwiseCorrelationFile.rds")

paperAnalysis = paperAnalysis[-c(2,3,6,8)]
ELAnalysis = ELAnalysis[-c(2,3,6,8)]

compareAnalyses = function(ogAnalysis, newAnalysis, prefix = NULL, pCuttof = 0.05){
  
  if(!is.null(prefix)){
    outputOGAnalysis = list()
    outputNewAnalysis = list()
    for(i in 1:length(prefix)){
      currentCols = grep(prefix[i], names(ogAnalysis))
      currentCols = currentCols[1:3]
      
      tempOGAnalysis = ogAnalysis[,currentCols]
      tempNewAnalysis = newAnalysis[,currentCols]
      
      tempOGAnalysis = list(tempOGAnalysis)
      names(tempOGAnalysis) = prefix[i]
      
      tempNewAnalysis = list(tempNewAnalysis)
      names(tempNewAnalysis) = prefix[i]
      
      
      outputOGAnalysis = append(outputOGAnalysis, tempOGAnalysis)
      outputNewAnalysis = append(outputNewAnalysis, tempNewAnalysis)
    }
    
    newAnalysis = outputNewAnalysis
    ogAnalysis = outputOGAnalysis
  }
  
  
  
  results <- tibble(
    index = integer(),
    correlation = numeric(),
    rhoRankcorrelation = numeric(),
    numMismatchedNAs = integer(),
    pCorrelation = numeric(),
    numMismatchedPs = integer(),
    padjCorrelation = numeric(),
    numMismatchedPadjs = integer(),
    totalOgSigPadjs = integer(), 
    totalNewSigPadjs = integer(),
    missingFraction = integer(),
    mismatchedNAs = I(list()),
    mismatchedPs = I(list()),
    mismatchedPadjs = I(list())
  )
  
  for (i in seq_along(ogAnalysis)) {
    
    correlation <- cor(ogAnalysis[[i]][[1]], newAnalysis[[i]][[1]], use = "complete.obs")
    
    mismiatchedNAs <- which(!is.na(ogAnalysis[[i]][[1]]) %in% is.na(newAnalysis[[i]][[1]]))
    
    ogAnalysis[[i]] = ogAnalysis[[i]] %>% mutate(rhoRank = rank(ogAnalysis[[i]][[1]])) 
    newAnalysis[[i]] = newAnalysis[[i]] %>% mutate(rhoRank = rank(newAnalysis[[i]][[1]]))
    
    rhoRankcorrelation <- cor(ogAnalysis[[i]][[4]], newAnalysis[[i]][[4]], use = "complete.obs")
    
    pCorrelation <- cor(ogAnalysis[[i]][[2]], newAnalysis[[i]][[2]], use = "complete.obs")
    
    mismiatchedPs <- which(ogAnalysis[[i]][[2]] < pCuttof)[which(!which(ogAnalysis[[i]][[2]] < pCuttof) %in% 
                                                                   which(newAnalysis[[i]][[2]] < pCuttof))]
    
    padjCorrelation <- cor(ogAnalysis[[i]][[3]], newAnalysis[[i]][[3]], use = "complete.obs")
    
    mismiatchedPadjs <- which(ogAnalysis[[i]][[3]] < pCuttof)[which(!which(ogAnalysis[[i]][[3]] < pCuttof) %in% 
                                                                      which(newAnalysis[[i]][[3]] < pCuttof))]
    
    totalOgSigPadjs = length(which(ogAnalysis[[i]][[3]] < pCuttof))
    totalNewSigPadjs = length(which(newAnalysis[[i]][[3]] < pCuttof))
    
    missingFraction = length(mismiatchedPadjs) / length(which(ogAnalysis[[i]][[3]] < pCuttof))
    
    results <- rbind(results, tibble(
      index = i,
      correlation = correlation,
      rhoRankcorrelation = rhoRankcorrelation,
      numMismatchedNAs = length(mismiatchedNAs),
      pCorrelation = pCorrelation,
      numMissingPs = length(mismiatchedPs),
      padjCorrelation = padjCorrelation,
      numMissingPadjs = length(mismiatchedPadjs),
      totalOgSigPadjs = totalOgSigPadjs,
      totalNewSigPadjs = totalNewSigPadjs,
      missingFraction = missingFraction,
      mismatchedNAs = list(mismiatchedNAs),
      missingPs = list(mismiatchedPs),
      missingPadjs = list(mismiatchedPadjs),
    ))
  }
  row.names(results) = names(ogAnalysis)
  return(results)
}




geneMissingInNewAnalysis = compareAnalyses(paperAnalysis, newAnalysis, c("HI", "HV", "IV"))

i=1

cor(paperAnalysis[[i]][[1]], ELAnalysis[[i]][[1]], use = "complete.obs") *  cor(paperAnalysis[[i]][[1]], ELAnalysis[[i]][[1]], use = "complete.obs")


correlation <- cor(paperAnalysis[[i]][[1]], ELAnalysis[[i]][[1]], use = "complete.obs") *  cor(paperAnalysis[[i]][[1]], ELAnalysis[[i]][[1]], use = "complete.obs")




#---------------------------------------------------------------------
# --- set up equal branch length ASR --- 
# --------------------------------------------------------------------

commonMainTreesBackup = commonMainTrees
commonMainTrees$masterTree$edge.length = rep(1, length(commonMainTrees$masterTree$edge.length))

mainTreesBackup = mainTrees
mainTrees$masterTree$edge.length = rep(1, length(mainTrees$masterTree$edge.length))


equalCategoricalTree = categoricalTree
mainCategoricalTree = readRDS("output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysisCategoricalTree.rds")

equalCategoricalTree$edge.length == mainCategoricalTree$edge.length

mismatchBranches = which(!equalCategoricalTree$edge.length == mainCategoricalTree$edge.length)
length(mismatchBranches)

edgelabels(edge = mismatchBranches, cex=0.4, col= "purple", frame="none")

#---------------------------------------------------------------------
# --- figure out gmt--- 
# --------------------------------------------------------------------

c2 = read.gmt("Data/c2.all.v2026.1.Hs.symbols.gmt")

hsSymbols = read.gmt("Data/EnrichmentHsSymbolsFile.gmt")

length(which(hsSymbols$geneset.names %in%c2$geneset.names))
length(hsSymbols$geneset.names)
#---------------------------------------------------------------------
# --- top genes look into --- 
# --------------------------------------------------------------------

geneData = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysiscombinedGeneResultsWithAlternates.rds")

order(geneData$`CH-PadjNumSignificantInclusive`)

importantFirstCols = c(2,3,4,19,21,22,23)
importantSecondCols = importantFirstCols +27 +27
importantThirdCols = importantSecondCols +27 +27
importantFourthCols = importantThirdCols +27 +27


importantCols = c(importantFirstCols, importantSecondCols, importantThirdCols, importantFourthCols)

crits = c(1,3,5,6)

crits+7

criticalCols = importantCols[c(crits, crits+7, crits+14, crits+21)]

geneData[order(geneData$`CH-PadjNumSignificantInclusive`), criticalCols]


importantData = geneData[,criticalCols]
names(importantData) = substr(names(importantData),1, 10)
importantData

importantData = importantData[c(1,5,9,13,3,7,11,15,2,6,10,14,4,8,12,16)]

importantData[order(importantData$`HV-RhoMed`, decreasing = T),]

HVONlyGenes = importantData[which(importantData$`HV-signifi` & !importantData$`CH-signifi`),]
HVONlyGenes[order(HVONlyGenes$`HV-RhoMedi`, decreasing = T),]

geneData[which(rownames(geneData) == "DNAH14"),]


#---------------------------------------------------------------------
# --- invert exclusive lookinto --- 
# --------------------------------------------------------------------

goData = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysiscombinedGOResultsWithAlternatesPermulated-KeggReactome.rds")
goData$ChRhoBoost = abs(goData$`CH-statMedian`) - abs(goData$`HI-statMedian`)
which(goData$`CH-significantRobust` | goData$`HI-significantRobust` | goData$`HV-significantRobust`)

sigGoData = goData[which(goData$`CH-significantRobust` | goData$`HI-significantRobust` | goData$`HV-significantRobust`),]


sigGoData[which(abs(sigGoData$`CH-stat`) < abs(sigGoData$`HI-stat`) & sigGoData$`HI-stat` > 0 & sigGoData$`CH-significantRobust`),c(5, 56)]


goData$ChRhoBoost = abs(goData$`CH-statMedian`) - abs(goData$`HI-statMedian`)
goData$HVRhoBoost = abs(goData$`HI-statMedian`) - abs(goData$`CH-statMedian`) 

sigGoData = sigGoData[order(sigGoData$ChRhoBoost, decreasing = T),]

sigGoData = sigGoData[order(sigGoData$ChRhoBoost, decreasing = T), c(5,6,56,57, 107,108)]


goData[which(rownames(goData) == "REACTOME_FATTY_ACIDS"),]


smGEnes = c("CYP17A1", "HSD11B1", "HSD17B2", "CYP21A2", "HSD17B3", "HSD3B2", "CYP11A1", "COMT", "SULT1E1", "CYP7B1", "CYP1A1", "HSD17B6")
#CYP17A1:11, HSD11B1:16, HSD17B2:37, CYP21A2:269, HSD17B3:287, COMT:342, HSD3B2:359, SULT1E1:479, CYP7B1:568, CYP1A1:576, HSD17B6:684, CYP11A1:795, AKR1D1:838, CYP7A1:1215, SRD5A3:2017, HSD11B2:2198, CYP19A1:2378, SULT2B1:8717.5, UGT1A6:8717.5, UGT1A1:8717.5, STS:8717.5, SRD5A1:8717.5, HSD17B12:8717.5, SRD5A2:8717.5, HSD17B1:8717.5, CYP1B1:8717.5, HSD17B8:8717.5, HSD17B7:14164

smeGeneData = importantData[which(rownames(importantData) %in% smGEnes),]

smeGeneData = smeGeneData[order(smeGeneData$`HV-RhoMedi` - smeGeneData$`HI-RhoMedi`, decreasing = T),]


sigGoData[1,]$`HV-statMedian`/sigGoData[1,]$`CH-statMedian`
#---------------------------------------------------------------------
# --- Update supplement Files --- 
# --------------------------------------------------------------------

permData = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysisPermulationsPValueCorrelations.rds")
mainData = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysiscombinedGeneResultsWithAlternates.rds")

HIPermData = permData[[2]][[1]][,c(4,5)]
names(HIPermData) = c("HI-permP", "HI-permP.adj")
HIPermData$`HI-permSignificant` = HIPermData$`HI-permP.adj` <0.05

HVPermData = permData[[2]][[4]][,c(4,5)]
names(HVPermData) = c("HV-permP", "HV-permP.adj")
HVPermData$`HV-permSignificant` = HVPermData$`HV-permP.adj` <0.05

IVPermData = permData[[2]][[5]][,c(4,5)]
names(IVPermData) = c("IV-permP", "IV-permP.adj")
IVPermData$`IV-permSignificant` = IVPermData$`IV-permP.adj` <0.05


combinedAltAndPermGene = mainData
combinedAltAndPermGene = combinedAltAndPermGene %>% add_column(HIPermData, 
                                                               .after ="HI-significantRobust")

combinedAltAndPermGene = combinedAltAndPermGene %>% add_column(HVPermData, 
                                                               .after ="HV-significantRobust")

combinedAltAndPermGene = combinedAltAndPermGene %>% add_column(IVPermData, 
                                                               .after ="IV-significantRobust")

saveRDS(combinedAltAndPermGene, "Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysiscombinedGeneResultsWithAlternatesPermulated.rds")
write.csv(combinedAltAndPermGene, "Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysiscombinedGeneResultsWithAlternatesPermulated.csv")


mainGOData = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysiscombinedGOResultsWithAlternates-KeggReactome.rds")
HIpermImport = readRDS("Output/ComplexDietCentralAnalysis/Herbivore-Insectivore/ComplexDietCentralAnalysisHerbivore-InsectivoreEnrichment-Permulation-KeggReactome.rds")
HVpermImport = readRDS("Output/ComplexDietCentralAnalysis/Herbivore-Vertivore/ComplexDietCentralAnalysisHerbivore-VertivoreEnrichment-Permulation-KeggReactome.rds")
IVpermImport = readRDS("Output/ComplexDietCentralAnalysis/Insectivore-Vertivore/ComplexDietCentralAnalysisInsectivore-VertivoreEnrichment-Permulation-KeggReactome.rds")

HIPermData = HIpermImport[[1]][,c(2,3)]
names(HIPermData) = c("HI-permP", "HI-permP.adj")
HIPermData$`HI-permSignificant` = HIPermData$`HI-permP.adj` <0.05

HVPermData = HVpermImport[[1]][,c(2,3)]
names(HVPermData) = c("HV-permP", "HV-permP.adj")
HVPermData$`HV-permSignificant` = HVPermData$`HV-permP.adj` <0.05

IVPermData = IVpermImport[[1]][,c(2,3)]
names(IVPermData) = c("IV-permP", "IV-permP.adj")
IVPermData$`IV-permSignificant` = IVPermData$`IV-permP.adj` <0.05


combinedAltAndPermGO = mainGOData
combinedAltAndPermGO = combinedAltAndPermGO %>% add_column(HIPermData, 
                                                           .after ="HI-significantRobust")

combinedAltAndPermGO = combinedAltAndPermGO %>% add_column(HVPermData, 
                                                           .after ="HV-significantRobust")

combinedAltAndPermGO = combinedAltAndPermGO %>% add_column(IVPermData, 
                                                           .after ="IV-significantRobust")

saveRDS(combinedAltAndPermGO, "Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysiscombinedGOResultsWithAlternatesPermulated-KeggReactome.rds")
write.csv(combinedAltAndPermGO, "Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysiscombinedGOResultsWithAlternatesPermulated-KeggReactome.csv")


# ------

mainTrees = readRDS("data/zoonomiaAllMammalsTrees.rds")
masterTree = mainTrees$masterTree

fourCatTree = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysisCategoricalTree.rds")
threeCatTree = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysisMergedCategoricalTree.rds")
equalTree = readRDS("Output/ComplexDietCentralAnalysisEqualLengthInference/ComplexDietCentralAnalysisEqualLengthInferenceCategoricalTree.rds")

BovTree = readRDS("Output/CladeBinaryBovidae/CladeBinaryBovidaeCategoricalTree.rds")
CerTree = readRDS("Output/CladeBinaryCervidae/CladeBinaryCervidaeCategoricalTree.rds")
CriTree = readRDS("Output/CladeBinaryCricetidae/CladeBinaryCricetidaeCategoricalTree.rds")
HysTree = readRDS("Output/CladeBinaryHystricognathi/CladeBinaryHystricognathiCategoricalTree.rds")
PerTree = readRDS("Output/CladeBinaryPeropdidae/CladeBinaryPeropdidaeCategoricalTree.rds")
VesTree = readRDS("Output/CladeBinaryVespertilionidae/CladeBinaryVespertilionidaeCategoricalTree.rds")

writeLines("----Master Tree----", "Output/ComplexDietCentralAnalysis/SupplementaryData5--Trees.txt")
write.tree(masterTree, file = "Output/ComplexDietCentralAnalysis/SupplementaryData5--Trees.txt", append = T)
cat("\n----Four Category Tree----\n", file = "Output/ComplexDietCentralAnalysis/SupplementaryData5--Trees.txt", append = T)
write.tree(fourCatTree, file = "Output/ComplexDietCentralAnalysis/SupplementaryData5--Trees.txt", append = T)
cat("\n----Three Category Tree----\n",  file ="Output/ComplexDietCentralAnalysis/SupplementaryData5--Trees.txt", append = T)
write.tree(threeCatTree, file = "Output/ComplexDietCentralAnalysis/SupplementaryData5--Trees.txt", append = T)

cat("\n----Alternate Reconstruction Tree----\n",  file ="Output/ComplexDietCentralAnalysis/SupplementaryData5--Trees.txt", append = T)
write.tree(equalTree, file = "Output/ComplexDietCentralAnalysis/SupplementaryData5--Trees.txt", append = T)

cat("\n\n----Clade Binary Trees----\n",  file ="Output/ComplexDietCentralAnalysis/SupplementaryData5--Trees.txt", append = T)
cat("\n----Bovidae Tree----\n",  file ="Output/ComplexDietCentralAnalysis/SupplementaryData5--Trees.txt", append = T)
write.tree(BovTree, file = "Output/ComplexDietCentralAnalysis/SupplementaryData5--Trees.txt", append = T)
cat("\n----Cervidae Tree----\n",  file ="Output/ComplexDietCentralAnalysis/SupplementaryData5--Trees.txt", append = T)
write.tree(CerTree, file = "Output/ComplexDietCentralAnalysis/SupplementaryData5--Trees.txt", append = T)
cat("\n----Cricetidae Tree----\n",  file ="Output/ComplexDietCentralAnalysis/SupplementaryData5--Trees.txt", append = T)
write.tree(CriTree, file = "Output/ComplexDietCentralAnalysis/SupplementaryData5--Trees.txt", append = T)
cat("\n----Hystricognathi Tree----\n",  file ="Output/ComplexDietCentralAnalysis/SupplementaryData5--Trees.txt", append = T)
write.tree(HysTree, file = "Output/ComplexDietCentralAnalysis/SupplementaryData5--Trees.txt", append = T)
cat("\n----Peropdidae Tree----\n",  file ="Output/ComplexDietCentralAnalysis/SupplementaryData5--Trees.txt", append = T)
write.tree(PerTree, file = "Output/ComplexDietCentralAnalysis/SupplementaryData5--Trees.txt", append = T)
cat("\n----Vespertilionidae Tree----\n",  file ="Output/ComplexDietCentralAnalysis/SupplementaryData5--Trees.txt", append = T)
write.tree(VesTree, file = "Output/ComplexDietCentralAnalysis/SupplementaryData5--Trees.txt", append = T)



cat("\n\n----Alternate Trees----\n",  file ="Output/ComplexDietCentralAnalysis/SupplementaryData5--Trees.txt", append = T)

for(i in 1:100){
  alternateTree = readRDS(paste0("Output/ComplexDietCentralAnalysis/Alternates/Alternate", i, "ComplexDietCentralAnalysisCategoricalTree.rds"))
  
  cat(paste("\n----Alternate", i, "Trees----\n"),  file ="Output/ComplexDietCentralAnalysis/SupplementaryData5--Trees.txt", append = T)
  write.tree(alternateTree, file = "Output/ComplexDietCentralAnalysis/SupplementaryData5--Trees.txt", append = T)

  
}


alternateTree1 = readRDS(paste0("Output/ComplexDietCentralAnalysis/Alternates/Alternate", 1, "ComplexDietCentralAnalysisCategoricalTree.rds"))
alternateTree2 = readRDS(paste0("Output/ComplexDietCentralAnalysis/Alternates/Alternate", 2, "ComplexDietCentralAnalysisCategoricalTree.rds"))

alternateTree1$tip.label %in% alternateTree2$tip.label




# --- binary results 


BovCor = readRDS("Output/CladeBinaryBovidae/0-1/CladeBinaryBovidae0-1CorrelationFile.rds")
CerCor = readRDS("Output/CladeBinaryCervidae/0-1/CladeBinaryCervidae0-1CorrelationFile.rds")
CriCor = readRDS("Output/CladeBinaryCricetidae/0-1/CladeBinaryCricetidae0-1CorrelationFile.rds")
HysCor = readRDS("Output/CladeBinaryHystricognathi/0-1/CladeBinaryHystricognathi0-1CorrelationFile.rds")
PerCor = readRDS("Output/CladeBinaryPeropdidae/0-1/CladeBinaryPeropdidae0-1CorrelationFile.rds")
VesCor = readRDS("Output/CladeBinaryVespertilionidae/0-1/CladeBinaryVespertilionidae0-1CorrelationFile.rds")


#---------------------------------------------------------------------
# --- Looking into msucle results --- 
# --------------------------------------------------------------------
moveIndexToNames = function(x){
  # split at ":"
  parts <- strsplit(x, ":")
  
  # keep only the part after ":"
  x <- sapply(parts, `[`, 1)
  
  # set names to the part before ":"
  names(x) <- sapply(parts, `[`, 2)
  x
}

CHStraitedGenes = c("DES", "ACTC1", "MYH8", "NEB")
names(CHStraitedGenes) = c(163, 231, 272, 289)
HIStraitedGenes = c("ES", "MYH8", "TCAP", "TMOD4", "NEB")
names(HIStraitedGenes) = c(76, 322, 343, 401, 564)

HVStraitedGenes = c("ACTC1:110", "DMD:368", "NEB:455", "MYH8:745")
HVStraitedGenes = moveIndexToNames(HVStraitedGenes)

CHMuscleGenes = c("CAMK2B:3", "ITPR1:97", "DES:163", "ACTC1:231", "MYH8:272", "CACNB2:280", "NEB:289", "CACNA1G:297")
CHMuscleGenes = moveIndexToNames(CHMuscleGenes)
HIMuscleGenes = c("CAMK2B:5", "DES:76", "ITPR1:85", "KCNJ12:212", "KCNJ14:239", "ATP1A3:245", "ATP2B3:269", "MYH8:322", "TCAP:343", "ATP2A2:350", "CAMK2D:385")
HIMuscleGenes = moveIndexToNames(HIMuscleGenes)
HVMuscleGenes = c("CAMK2B:20", "ACTC1:110", "CACNA1G:135", "KCNH2:211", "CACNB2:263", "ACTG2:298", "SLC8A3:340", "DMD:368", "ACTA2:381", "NEB:455", "SCN5A:495")
HVMuscleGenes = moveIndexToNames(HVMuscleGenes)

CHCardiacGenes = c("CAMK2B:3", "ITPR1:97", "CACNB2:280")
CHCardiacGenes = moveIndexToNames(CHCardiacGenes)
HICardiacGenes = c("CAMK2B:5", "ITPR1:85", "KCNJ12:212", "KCNJ14:239", "ATP1A3:245", "ATP2B3:269", "ATP2A2:350") 
HICardiacGenes = moveIndexToNames(HICardiacGenes)
HVCardiacGenes = c("CAMK2B:20", "KCNH2:211", "CACNB2:263", "SLC8A3:340", "SCN5A:495")
HVCardiacGenes = moveIndexToNames(HVCardiacGenes)


allGenes = list(CHStraitedGenes, HIStraitedGenes, HVStraitedGenes, CHMuscleGenes, HIMuscleGenes, HVMuscleGenes, CHCardiacGenes, HICardiacGenes, HVCardiacGenes)
names(allGenes) = c("CHStraitedGenes", "HIStraitedGenes", "HVStraitedGenes", "CHMuscleGenes", "HIMuscleGenes", "HVMuscleGenes", "CHCardiacGenes", "HICardiacGenes", "HVCardiacGenes")


geneslist = sort(unique(unlist(allGenes)))

geneTable = matrix(nrow = length(geneslist), ncol = length(allGenes))
geneTable = data.frame(geneTable)


rownames(geneTable) = geneslist
colnames(geneTable) = names(allGenes)

for(i in 1:ncol(geneTable)){
  geneTable[i] = names(allGenes[[i]])[match(rownames(geneTable), allGenes[[i]])]
}

geneTable$NumberOfSets <- rowSums(!is.na(geneTable))


geneTable = geneTable[order(geneTable$NumberOfSets, decreasing = T),]

#---------------------------------------------------------------------
# --- looking into permualtions --- 
# --------------------------------------------------------------------
library(tibble)
library(ggvenn)
library(gridExtra)



permData = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysisPermulationsPValueCorrelations.rds")
mainData = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysiscombinedGeneResultsWithAlternates.rds")

HIPermData = permData[[2]][[1]][,c(4,5)]
names(HIPermData) = c("HI-permP", "HI-permP.adj")
HIPermData$`HI-permSignificant` = HIPermData$`HI-permP.adj` <0.05

HVPermData = permData[[2]][[4]][,c(4,5)]
names(HVPermData) = c("HV-permP", "HV-permP.adj")
HVPermData$`HV-permSignificant` = HVPermData$`HV-permP.adj` <0.05

IVPermData = permData[[2]][[5]][,c(4,5)]
names(IVPermData) = c("IV-permP", "IV-permP.adj")
IVPermData$`IV-permSignificant` = IVPermData$`IV-permP.adj` <0.05


combinedAltAndPermGene = mainData
combinedAltAndPermGene = combinedAltAndPermGene %>% add_column(HIPermData, 
                                                       .after ="HI-significantRobust")

combinedAltAndPermGene = combinedAltAndPermGene %>% add_column(HVPermData, 
                                                       .after ="HV-significantRobust")

combinedAltAndPermGene = combinedAltAndPermGene %>% add_column(IVPermData, 
                                                       .after ="IV-significantRobust")



mainGOData = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysiscombinedGOResultsWithAlternates-KeggReactome.rds")
HIpermImport = readRDS("Output/ComplexDietCentralAnalysis/Herbivore-Insectivore/ComplexDietCentralAnalysisHerbivore-InsectivoreEnrichment-Permulation-KeggReactome.rds")
HVpermImport = readRDS("Output/ComplexDietCentralAnalysis/Herbivore-Vertivore/ComplexDietCentralAnalysisHerbivore-VertivoreEnrichment-Permulation-KeggReactome.rds")
IVpermImport = readRDS("Output/ComplexDietCentralAnalysis/Insectivore-Vertivore/ComplexDietCentralAnalysisInsectivore-VertivoreEnrichment-Permulation-KeggReactome.rds")

HIPermData = HIpermImport[[1]][,c(2,3)]
names(HIPermData) = c("HI-permP", "HI-permP.adj")
HIPermData$`HI-permSignificant` = HIPermData$`HI-permP.adj` <0.05

HVPermData = HVpermImport[[1]][,c(2,3)]
names(HVPermData) = c("HV-permP", "HV-permP.adj")
HVPermData$`HV-permSignificant` = HVPermData$`HV-permP.adj` <0.05

IVPermData = IVpermImport[[1]][,c(2,3)]
names(IVPermData) = c("IV-permP", "IV-permP.adj")
IVPermData$`IV-permSignificant` = IVPermData$`IV-permP.adj` <0.05


combinedAltAndPermGO = mainGOData
combinedAltAndPermGO = combinedAltAndPermGO %>% add_column(HIPermData, 
                                                           .after ="HI-significantRobust")

combinedAltAndPermGO = combinedAltAndPermGO %>% add_column(HVPermData, 
                                                           .after ="HV-significantRobust")

combinedAltAndPermGO = combinedAltAndPermGO %>% add_column(IVPermData, 
                                                           .after ="IV-significantRobust")





currentData = permData[[2]][[1]]
currentData = permData[[2]][[4]]
currentData = permData[[2]][[5]]

signedLogBaseP = log10(currentData$P) * sign(currentData$Rho) *-1
signedLogPermP = log10(currentData$permP) * sign(currentData$Rho) *-1
cor(signedLogBaseP, signedLogPermP, method = "spearman", use = "complete.obs")
plot1 = plot(signedLogBaseP, signedLogPermP)
plot2 = plot(signedLogBaseP, signedLogPermP)
plot3 = plot(signedLogBaseP, signedLogPermP)

plot1






combinedAltAndPerm = combinedAltAndPermGene
combinedAltAndPerm = combinedAltAndPermGO


convertColnameToIndex = function(name){
  output = which(colnames(combinedAltAndPerm)==name)
  output
}

prefix = "HI"
prefix = "HV"
prefix = "IV"

prefix = prefix
robustSig = convertColnameToIndex(paste0(prefix, "-significantRobust"))
permSig = convertColnameToIndex(paste0(prefix, "-permSignificant"))
mainSig = convertColnameToIndex(paste0(prefix, "-significant"))
altSig = convertColnameToIndex(paste0(prefix, "-PadjNumSignificantInclusive"))


venn_list <- list(
  Column1 = which(combinedAltAndPerm[,altSig] > 50),
  Column2 = which(combinedAltAndPerm[,permSig]),
  Column3 = which(combinedAltAndPerm[,mainSig])
  
)

makeVenn = function(infix){
  robustSig = convertColnameToIndex(paste0(infix, "-significantRobust"))
  permSig = convertColnameToIndex(paste0(infix, "-permSignificant"))
  mainSig = convertColnameToIndex(paste0(infix, "-significant"))
  altSig = convertColnameToIndex(paste0(infix, "-PadjNumSignificant"))
  
  
  venn_list <- list(
    Alternates = which(combinedAltAndPerm[,altSig] > 50),
    Permulations = which(combinedAltAndPerm[,permSig]),
    Main = which(combinedAltAndPerm[,mainSig])
    
  )
  vennout = ggvenn(venn_list) +  labs(title = infix)
}


venn1 = makeVenn("HI")
venn2 = makeVenn("HV")
venn3 = makeVenn("IV")

grid.arrange(venn1, venn2, venn3, nrow = 2)


combinedAltAndPermGO[which(combinedAltAndPermGO$`IV-significantRobust` & !combinedAltAndPermGO$`IV-permSignificant`),]


convertColnameToIndexGene = function(name){
  output = which(colnames(combinedAltAndPermGene)==name)
  output
}

hiPermCol = convertColnameToIndexGene(paste0("HI-permP.adj"))
hvPermCol = convertColnameToIndexGene(paste0("HV-permP.adj"))
ivPermCol = convertColnameToIndexGene(paste0("IV-permP.adj"))



combinedAltAndPermGene
selectedGene = c("CLDN16", "CPB1", "SLC14A2", "ACADSB", "PNLIP", "SLC13A2") 
selectedRow = which(rownames(combinedAltAndPermGene) %in% selectedGene)

combinedAltAndPermGene[selectedRow,c(hiPermCol, hvPermCol, ivPermCol)]





# all.equal(rownames(mainData), rownames(permData[[1]]))



#---------------------------------------------------------------------
# --- Working with alex data  --- 
# --------------------------------------------------------------------

length(which(raw_data$`IV-significantRobust`))


bforePurne = node_data$pathway
afterPrune = node_data$pathway
which(!bforePurne %in% afterPrune)

bforePurne[21]
#---------------------------------------------------------------------
# --- Debugging making maingRObust venn diagrams   --- 
# --------------------------------------------------------------------



length(which(geneSignificanceResults$`CH-significantRobust`))

which(GoSignificanceResults$`HV-significant` & GoSignificanceResults$`HI-significant` & !GoSignificanceResults$`CH-significant`)

GoSignificanceResults[635,]

length(which(geneSignificanceResults$`IV-significantRobust`))

length(which(GoSignificanceResults$`IV-significantRobust`))

length(which(GoCombinedResults$`IV-PadjNumSignificant` > 50))

length(which(GoCombinedResults$`IV-PadjMedian` <= 0.05))
length(which(GoCombinedResults$`IV-PadjMean` <= 0.05))



convertColnameToIndex = function(name){
  output = which(colnames(combinedResults)==name)
  output
}

convertColnameToIndexGO = function(name){
  output = which(colnames(GoCombinedResults)==name)
  output
}


altOnlyGoIndexes =  which(GoCombinedResults$`IV-PadjNumSignificant` > 50)[!which(GoCombinedResults$`IV-PadjNumSignificant` > 50) %in% which(GoSignificanceResults$`IV-significantRobust`)]

robustGO = GoCombinedResults[which(GoSignificanceResults$`IV-significantRobust`),c(convertColnameToIndexGO("IV-stat"), convertColnameToIndexGO("IV-significantRobust"),convertColnameToIndexGO("IV-PadjNumSignificant"))]

altOnlyGo = GoCombinedResults[altOnlyGoIndexes,c(convertColnameToIndexGO("IV-stat"), convertColnameToIndexGO("IV-significantRobust"),convertColnameToIndexGO("IV-PadjNumSignificant"))]

nrow(robustGO)
#---------------------------------------------------------------------
# --- Making column that includes main analysis as an alternate   ---DFDF
# --------------------------------------------------------------------
combinedGeneDataFilename = "Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysiscombinedGeneResultsWithAlternates.rds"
combinedGeneDataFilename = "Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysiscombinedGOResultsWithAlternates-KeggReactome.rds"

combinedResults = readRDS(combinedGeneDataFilename)


prefixes <- unique(sub("-.*", "", names(combinedResults)))

for (prefix in prefixes) {
  
  p_col <- paste0(prefix, "-P")
  padj_col <- paste0(prefix, "-p.adj")
  
  pnum_col <- paste0(prefix, "-PNumSignificant")
  padjnum_col <- paste0(prefix, "-PadjNumSignificant")
  
  # skip if columns don't exist
  if (!all(c(p_col, padj_col, pnum_col, padjnum_col) %in% names(combinedResults))) next
  
  # new column names
  pnum_new <- paste0(prefix, "-PNumSignificant_Inclusive")
  padjnum_new <- paste0(prefix, "-PadjNumSignificant_Inclusive")
  
  # create new values
  combinedResults[[pnum_new]] <- combinedResults[[pnum_col]] + 
    ifelse(combinedResults[[p_col]] < 0.05, 1, 0)
  
  combinedResults[[padjnum_new]] <- combinedResults[[padjnum_col]] + 
    ifelse(combinedResults[[padj_col]] < 0.05, 1, 0)
  
  # move columns to correct position (after originals)
  p_index <- match(pnum_col, names(combinedResults))
  padj_index <- match(padjnum_col, names(combinedResults))
  
  # reorder for PNum
  combinedResults <- combinedResults[, append(
    names(combinedResults)[-which(names(combinedResults) == pnum_new)],
    pnum_new,
    after = p_index
  )]
  
  # reorder for PadjNum
  combinedResults <- combinedResults[, append(
    names(combinedResults)[-which(names(combinedResults) == padjnum_new)],
    padjnum_new,
    after = padj_index
  )]
}

saveRDS(combinedResults, combinedGeneDataFilename)

#---------------------------------------------------------------------
# --- Looking into the the differences in venn diagrams and generally the difference in alternate genes  --- 
# --------------------------------------------------------------------
library(ggvenn)
combinedGeneDataFilename = "Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysiscombinedGeneResultsWithAlternates.rds"
combinedResults = readRDS(combinedGeneDataFilename)

vespertDataFilename = "Output/CladeBinaryVespertilionidae/CladeBinaryVespertilionidaePairwiseCorrelationFile.rds"
vespertResults = readRDS(vespertDataFilename)
vespertResults = vespertResults[[1]]

significanceColumns = names(combinedResults)[grep("PadjNumSignificant", names(combinedResults))]
for(i in significanceColumns){
  currentCol = combinedResults[,names(combinedResults) == i]
  currentCol = (currentCol > 50)
  combinedResults = cbind(combinedResults, currentCol)
  
  currentName = i 
  currentName = gsub("PadjNumSignificant", "AlternateSignificant", i)
  names(combinedResults)[ncol(combinedResults)] = currentName
}

convertColnameToIndex = function(name){
  output = which(colnames(combinedResults)==name)
  output
}


ChMainSig = convertColnameToIndex(paste0("CH", "-significant"))
ChAltSig = convertColnameToIndex(paste0("CH", "-AlternateSignificant_Inclusive"))

prefix = "HI"

mainSig = convertColnameToIndex(paste0(prefix, "-significant"))
altSig = convertColnameToIndex(paste0(prefix, "-AlternateSignificant_Inclusive"))
pAdj = convertColnameToIndex(paste0(prefix, "-p.adj"))
numAlts = convertColnameToIndex(paste0(prefix, "-PadjNumSignificant_Inclusive"))


{
  DiffGenes = which(combinedResults[,mainSig] != combinedResults[,altSig])

  venn_list <- list(
    Column1 = which(combinedResults[,mainSig]),
    Column2 = which(combinedResults[,altSig])
  )
  ggvenn(venn_list)    
  
  mainOnlyGenes = DiffGenes[which(combinedResults[DiffGenes, mainSig])]
  altOnlyGenes = DiffGenes[which(combinedResults[DiffGenes, altSig])]
  sharedGenes = which(combinedResults[,mainSig] & combinedResults[,altSig] )
  allGenes =    which(combinedResults[,mainSig] | combinedResults[,altSig] )
  allMain =     which(combinedResults[,mainSig])
  allAlt =      which(combinedResults[,altSig] )
  
  combinedResults[mainOnlyGenes, pAdj]

  plot(combinedResults[mainOnlyGenes, pAdj])
  plot(combinedResults[sharedGenes, pAdj])
  plot(combinedResults[allGenes,pAdj], combinedResults[allGenes,numAlts])
  plot(combinedResults[mainOnlyGenes,pAdj], combinedResults[mainOnlyGenes,numAlts])
  plot(combinedResults[,pAdj], combinedResults[,numAlts])
  plot(combinedResults[allMain,pAdj], combinedResults[allMain,numAlts]) 
  
  cor(combinedResults[,pAdj], combinedResults[,numAlts], use = "complete.obs", method = "pearson")
  cor(combinedResults[,pAdj], combinedResults[,numAlts], use = "complete.obs", method = "spearman")
  
  vespertResults$category = rep(NA, nrow(vespertResults))
  vespertResults$category = ifelse(combinedResults[, mainSig] & combinedResults[, altSig], "both",
                                   ifelse(combinedResults[, mainSig], "main_only",
                                          ifelse(combinedResults[, altSig], "alt_only", "neither")))
  
  vespertResults$categoryByCH = rep(NA, nrow(vespertResults))
  vespertResults$categoryByCH = ifelse(combinedResults[, mainSig] & combinedResults[, altSig] & combinedResults[, ChMainSig] & combinedResults[, ChAltSig], "both (+CH)",
                                   ifelse(combinedResults[, mainSig] & combinedResults[, ChMainSig], "main_only (+CH)",
                                          ifelse(combinedResults[, altSig] & combinedResults[, ChAltSig], "alt_only (+CH)", 
                                                 ifelse(combinedResults[, mainSig] & combinedResults[, altSig], "both",
                                                        ifelse(combinedResults[, mainSig], "main_only",
                                                               ifelse(combinedResults[, altSig], "alt_only", 
                                                                      "none"))))))
  

    plot1 = ggplot(vespertResults, aes(x=category, y=p.adj))+
    geom_violin(adjust=1/3) +
    geom_jitter(position=position_jitter(0.2)) +
    theme_classic()+
    ggtitle(prefix)
  
    
    plot2 = ggplot(vespertResults, aes(x=category, y=Rho))+
      geom_violin(adjust=1/3) +
      geom_jitter(position=position_jitter(0.2)) +
      theme_classic()+
      ggtitle(prefix) 
    
    plot3 = ggplot(vespertResults, aes(x=categoryByCH, y=p.adj))+
      geom_violin(adjust=1/3) +
      geom_jitter(position=position_jitter(0.2)) +
      theme_classic()+
      ggtitle(prefix)
    
    
    plot4 = ggplot(vespertResults, aes(x=categoryByCH, y=Rho))+
      geom_violin(adjust=1/3) +
      geom_jitter(position=position_jitter(0.2)) +
      theme_classic()+
      ggtitle(prefix) 
    
    plot5 = ggplot(vespertResults, aes(x=factor(categoryByCH, levels = c("alt_only", "both", "main_only", "alt_only (+CH)", "both (+CH)", "main_only (+CH)", "none", "NA")), y=Rho))+
                     geom_violin(adjust=1/3) +
                     geom_jitter(position=position_jitter(0.2)) +
                     theme_classic()+
                     ggtitle(prefix) 
    
    plot5
}




# -- Look into vespert results --- 








CHDiffGenes = which(combinedResults$`CH-significant` != combinedResults$`CH-AlternateSignificant`)
HIDiffGenes = which(combinedResults$`HI-significant` != combinedResults$`HI-AlternateSignificant`)
HVDiffGenes = which(combinedResults$`HV-significant` != combinedResults$`HV-AlternateSignificant`)





length(CHDiffGenes)
which(combinedResults[CHDiffGenes, convertColnameToIndex("CH-significant")])
which(combinedResults[CHDiffGenes, convertColnameToIndex("CH-alternatesignificant")])



View(combinedResults[CHDiffGenes,])


library(ggvenn)


# Prepare list
venn_list <- list(
  Column1 = which(combinedResults$`CH-significant`),
  Column2 = which(combinedResults$`CH-AlternateSignificant`)
)

# Plot
ggvenn(venn_list)


#seeing if any genes significant in exactly 49 alternates 

prefix = "CH"
numAlts = convertColnameToIndex(paste0(prefix, "-PadjNumSignificant"))
which(combinedResults[numAlts] == 50)

combinedResults[c(116,161,1137,1471), mainSig]
which()

#---------------------------------------------------------------------
# --- Looking into insvertivore sorting   --- 
# --------------------------------------------------------------------

mergedData = read.csv("Data/mergedData.csv")

insVertivores = grep("InsVertivore", mergedData$DerekDietClassification90InsVertivoreSorting)



relevantData = mergedData[insVertivores,c(1:3, 22, 26:35)]

sortCargories = unique(relevantData$DerekDietClassification90InsVertivoreSorting)

sortedInvertivore = grep("InsVertivore-Insec", relevantData$DerekDietClassification90InsVertivoreSorting)
sortedVertivore = grep("InsVertivore-Carn", relevantData$DerekDietClassification90InsVertivoreSorting)
sortedMixed = grep("InsVertivore-Mix", relevantData$DerekDietClassification90InsVertivoreSorting)
sortedPisc = grep("InsVertivore-Pisc", relevantData$DerekDietClassification90InsVertivoreSorting)


relevantData[sortedInvertivore,c(1:2,5:10)]
relevantData[sortedVertivore, c(1:2,5:10)]
relevantData[sortedMixed,c(1:2,5:10)]
relevantData[sortedPisc,c(1:2,5:10)]


View(relevantData)
#---------------------------------------------------------------------
# --- Fast method to run DisplayCategoricalRERTree  --- 
# --------------------------------------------------------------------

source("Src/Reu/ZonomNameConvertVectorCommon.R")
source("Src/Loc/Dev/DisplayCategoricalRERTree.R")




RERFileName = paste(outputFolderName, filePrefix, "RERFile.rds", sep= "")       #Set a filename for the RERs based on the prefix
RERObject = readRDS(RERFileName)

pathsFileName = paste(outputFolderName, filePrefix,"CategoricalPathsFile.rds", sep= "")
pathsObject = readRDS(pathsFileName)


difGene = "CLDN16"

tipColumn = "ZoonomiaTip"
annotationLocation = "Data/MergedData.csv"

treesObj = mainTrees
rermat = RERObject
phenv = pathsObject
index = difGene


rer.cex = 0.7
tip.cex = 0.7
nalab = 'NA'
plot = T
subsetTree = T
equalLengths=T
minWidth=1
maxWidth=6
library(geiger)





pdf("testRERTree.pdf", height = 14, width = 14)
displayCategoricalRERTree(treesObj = mainTrees, rermat = RERObject, index = difGene, phenv = pathsObject, subsetTree = T, tipCol = tipColumn, annotLocation = annotationLocation)
dev.off()

returnRersAsTree(mainTrees, wynnContinousRER, difGene, wynnContinousPath)





#---------------------------------------------------------------------
# --- Making tree for leah  --- 
# --------------------------------------------------------------------

leahTree = read.tree("Data/110_species_tree.newick")

harshalMainTree = read.newick("Data/roadies_v1.1.16b.nwk")
leahFakeMainTrees = list()
leahFakeMainTrees$masterTree = leahTree
saveRDS(leahFakeMainTrees, "Data/LeahFakeMainTrees.rds")



#---------------------------------------------------------------------
# --- looking into ultrametric tree conversion --- 
# --------------------------------------------------------------------

demoTree = read.tree(("../../MiscData/maximum_likelihood_tree.txt"))

plotTree(demoTree)


#---------------------------------------------------------------------
# --- making plots for wynn  --- 
# --------------------------------------------------------------------

# Load in Data 
wynnCategoricalRER = readRDS("Output/HarshalCategoricalRER/HarshalCategoricalRERRERFile.rds")
wynnCategoricalPath = readRDS("Output/HarshalCategoricalRER/HarshalCategoricalRERCategoricalPathsFile.rds")
wynnCategoricalRERCommon = wynnCategoricalRER
colnames(wynnCategoricalRERCommon) = ZonomNameConvertVectorCommon(colnames(wynnCategoricalRERCommon),  annotationLocation = "Data/VGP_Mammals_Diet.csv", tipColumn = "Accession")

wynnContinousRER = readRDS("Output/HarshalContinousRERMod/HarshalContinousRERModRERFile.rds")
wynnContinousPath = readRDS("Output/HarshalContinousRERMod/HarshalContinousRERModContinuousPathsFile.rds")
wynnContinousRERCommon = wynnContinousRER
colnames(wynnContinousRERCommon) = ZonomNameConvertVectorCommon(colnames(wynnContinousRERCommon), annotationLocation = "Data/VGP_Mammals_Diet.csv", tipColumn = "Accession")

difGene = "NCE.ALDH1A1.subset_aln_region3568_start356701_w300.fa.filt"
palette(c( "#cc6677", "#117733","#33bbee", "white"))


#Make the RER Plot 
png("Output/Misc/WynnCategoricalPlot.png", width = 2000, height = 2000)
plotRers(wynnCategoricalRERCommon, difGene, wynnCategoricalPath, sortrers = T)
dev.off()


# Make the Scatter plot 
phenotypeTree = readRDS("Output/HarshalCategoricalRER/HarshalCategoricalRERCategoricalTree.rds")
speciesFilter = readRDS("Output/HarshalContinousRERMod/HarshalContinousRERModSpeciesFilter.rds")

RelativeEvolutionaryRate = wynnContinousRERCommon[which(rownames(wynnContinousRERCommon) == difGene),]
ChangeinTMM = wynnContinousPath
matchedPathsObject = tree2Paths(phenotypeTree, mainTrees, useSpecies = speciesFilter, categorical = TRUE) #do not binarize; the categorical data is already contained in the phenotype tree.

scatterPlotData = data.frame(RelativeEvolutionaryRate, ChangeinTMM, RERNames = names(RelativeEvolutionaryRate), categoricalPath = matchedPathsObject)


scatterPLot = ggplot(data = scatterPlotData, aes(x = ChangeinTMM, y = RelativeEvolutionaryRate, color = factor(categoricalPath))) +
  geom_point() +
  scale_color_manual(values = palette())+ 
  theme_minimal()
  #+ geom_text(aes(label = RERNames))


scatterPlot = ggplot(data = scatterPlotData, aes(x = ChangeinTMM, y = RelativeEvolutionaryRate, color = factor(categoricalPath))) +
  geom_point() +
  scale_color_manual(
    name = "Diet",                                # legend title
    breaks = c(1,2,3,4),                             # only show first 3
    labels = c("Carnivore", "Herbivore", "Omnivore", ""),
    values = setNames(palette()[1:4], c(1,2,3,4))     # match colors to those levels
  )+ 
  theme_minimal()
#+ geom_text(aes(label = RERNames)) #Optional adding of name labels 



png("Output/Misc/WynnScatterPlot.png", width = 1000, height = 1000)
scatterPlot
dev.off()












#---------------------------------------------------------------------
# --- Making RERTree for wynn   --- 
# --------------------------------------------------------------------

source("Src/Reu/ZonomNameConvertVectorCommon.R")


mainTrees = readRDS("data/HarshalFakeMainTrees.rds")
mainTrees = readRDS("data/NCE.ALDH1A1.filt.trees.rds")



RERObject = wynnContinousRER
pathsObject = wynnContinousPath

RERObject = wynnCategoricalRER
pathsObject = wynnCategoricalPath
RERObject = wynnCategoricalRERCommon
RERObject = wynnContinousRERCommon
#source("Src/Reu/RERConvergeFunctions.R")

difGene = "NCE.ALDH1A1.subset_aln_region3568_start356701_w300.fa.filt"



tipColumn = "Accession"
tipCol = "Accession"
annotLocation = "Data/VGP_Mammals_Diet.csv"

treesObj = mainTrees
rermat = RERObject
phenv = pathsObject
index = difGene

displayCategoricalRERTree = function(treesObj, rermat, index, phenv = NULL, subsetTree = T, equalLengths = T, minWidth = 1, maxWidth = 6, tipCol = "tipColumn", annotLocation = "Data/mergedData.csv"){
  treesObj$trees[[index]]$tip.label = ZonomNameConvertVectorCommon(treesObj$trees[[index]]$tip.label, tipColumn = tipCol, annotationLocation = annotLocation)
  treesObj$masterTree$tip.label = ZonomNameConvertVectorCommon(treesObj$masterTree$tip.label, tipColumn = tipCol, annotationLocation = annotLocation)
  colnames(rermat) = ZonomNameConvertVectorCommon(colnames(rermat), tipColumn = tipCol)
  returnRersAsTreeNew(treesObj, rermat, index, phenv, 0.7, 0.7, 'NA', T, subsetTree, equalLengths, minWidth, maxWidth)
  
  
  
  
}

rer.cex = 0.7
tip.cex = 0.7
nalab = 'NA'
plot = T
subsetTree = T
equalLengths=T
minWidth=1
maxWidth=6
library(geiger)

treesObj$trees[[index]]


returnRersAsTreeNew <- function(treesObj, rermat, index, phenv = NULL, rer.cex = 0.7,
                                tip.cex = 0.7, nalab = 'NA', plot = T, subsetTree = F, equalLengths=T, minWidth=1, maxWidth=6){
  trgene <- treesObj$trees[[index]]
  if (subsetTree == TRUE) {
    #Use phenv to identify non-NA paths and then subset gene tree by species
    subsp <- unique(colnames(rermat)[which(!is.na(phenv))])
    trgene <- pruneTree(trgene, subsp) #drops only species not in subsp
  }
  if(equalLengths){trgene$edge.length <- rep(2,nrow(trgene$edge))}
  ee=edgeIndexRelativeMaster(trgene, treesObj$masterTree)
  ii= treesObj$matIndex[ee[, c(2,1)]]
  rertree=rermat[index,ii]
  rertree[is.nan(rertree)]=NA #replace NaNs from C functions
  relativeRER = abs(scale(rertree))
  if(any(relativeRER > 3, na.rm=T)){
    message("WARNING: Potential misannotated gene!")
    for(i in which(relativeRER > 3)){
      if(!is.na(rownames(relativeRER)[i])){
        message(paste("Potential Misannotated gene at branch:", i ))
        message(paste("Species of branch=:", rownames(relativeRER)[i]))
        if(plot){
          message("Branch highlighted in PINK.")
        }
      }
    }
  }
  
  if (plot) {
    clampedRelativeRER = relativeRER
    clampedRelativeRER[order(relativeRER, decreasing = TRUE)[1:3]] <- sort(relativeRER, decreasing = TRUE)[4] #Set the top three high RERs to the fourth highest RER; this matches the calculation ignoring top values and prevents mis-IDed genes from throwing off the scale. 
    rerWidth = rescale(clampedRelativeRER, c(minWidth,maxWidth))
    rerWidth[is.na(rerWidth)] = 1
    
    par(mar = c(1,1,1,0))
    edgcols <- rep('black', nrow(trgene$edge))
    edgwds <- rep(1, nrow(trgene$edge))
    if(!is.null(phenv)){
      edgcols <- rep('black', nrow(trgene$edge))
      edgwds <- rerWidth
      if(length(unique(phenv) < length(palette()))){ # add a catch for continuous phenotypes and not run it in that case
        for(j in unique(phenv)[!is.na(unique(phenv))]){
          edgcols[phenv[ii]==j] <- palette()[j]
        }      
      }
    }
    if(any(relativeRER > 3, na.rm=T)){
      for(i in which(relativeRER > 3)){
        if(!is.na(rownames(relativeRER)[i])){
          edgcols[i] = 'hotpink'
        }
      }
    }
    plot.phylo(trgene, font = 2, edge.color = edgcols, edge.width = edgwds, cex = tip.cex)
    rerlab <- round(rertree,3)
    rerlab[is.na(rerlab)] <- nalab
    if(!is.null(phenv)){ #reset the color so that the RER labels still have the correct palette color even with the pink warning
      for(j in unique(phenv)[!is.na(unique(phenv))]){
        edgcols[phenv[ii]==j] <- palette()[j]
      }
    }
    edgelabels(rerlab, bg = NULL, adj = c(0.5,0.9), col = edgcols, frame = 'none',cex = rer.cex, font =2)
  }
  
  
  trgene$edge.length <- rertree
  return(trgene)
}

{
  edgeIndexRelativeMaster= function(tree, masterTree){
    map=matchAllNodes(tree,masterTree)
    newedge=tree$edge
    newedge[,1]=map[newedge[,1],2]
    newedge[,2]=map[newedge[,2],2]
    newedge
  }
  matchAllNodes = function(tree1, tree2){
    map=matchNodesInject(tree1,tree2)
    map=map[order(map[,1]),]
    map
  }
  matchNodesInject = function (tr1, tr2){
    if(length(tmpsp<-setdiff(tr1$tip.label, tr2$tip.label))>0){
      #stop(paste(paste(tmpsp, ","), "in tree1 do not exist in tree2"))
      stop(c("The following species in tree1 do not exist in tree2: ",paste(tmpsp, ", ")))
    }
    commontiplabels <- intersect(tr1$tip,tr2$tip)
    if(RF.dist(pruneTree(tr1,commontiplabels),pruneTree(tr2,commontiplabels))>0){
      stop("Discordant tree topology detected - gene/trait tree and treesObj$masterTree have irreconcilable topologies")
    }
    #if(RF.dist(tr1,tr2)>0){
    #  stop("Discordant tree topology detected - trait tree and treesObj$masterTree have irreconcilable topologies")
    #}
    
    toRm=setdiff(tr2$tip.label, tr1$tip.label)
    desc.tr1 <- lapply(1:tr1$Nnode + length(tr1$tip), function(x) extract.clade(tr1,
                                                                                x)$tip.label)
    names(desc.tr1) <- 1:tr1$Nnode + length(tr1$tip)
    desc.tr2 <- lapply(1:tr2$Nnode + length(tr2$tip), function(x) extract.clade(tr2,
                                                                                x)$tip.label)
    names(desc.tr2) <- 1:tr2$Nnode + length(tr2$tip)
    Nodes <- matrix(NA, length(desc.tr1), 2, dimnames = list(NULL,
                                                             c("tr1", "tr2")))
    for (i in 1:length(desc.tr1)) {
      Nodes[i, 1] <- as.numeric(names(desc.tr1)[i])
      for (j in 1:length(desc.tr2)) if (all(desc.tr1[[i]] %in%
                                            desc.tr2[[j]]))
        Nodes[i, 2] <- as.numeric(names(desc.tr2)[j])
    }
    
    iim=match(tr1$tip.label, tr2$tip.label)
    Nodes=rbind(cbind(1:length(tr1$tip.label),iim),Nodes)
    if(any(table(Nodes[,2])>1)){
      stop("Incorrect pseudorooting detected - use fixPseudoroot() function to correct trait tree topology")
    }
    
    Nodes
  }
  
}

pdf("WynnCategoricalTree.pdf", height = 14, width = 14)
displayCategoricalRERTree(treesObj = mainTrees, rermat = RERObject, index = difGene, phenv = pathsObject, subsetTree = T, tipCol = tipColumn, annotLocation = annotationLocation)
dev.off()

returnRersAsTree(mainTrees, wynnContinousRER, difGene, wynnContinousPath)

#---------------------------------------------------------------------
# --- helping wynn  --- 
# --------------------------------------------------------------------

harshalMainTree = read.newick("Data/roadies_v1.1.16b.nwk")
HarshalFakeMainTrees = list()
HarshalFakeMainTrees$masterTree = harshalMainTree
saveRDS(HarshalFakeMainTrees, "Data/HarshalFakeMainTrees.rds")

harhsalCategoricalTree = readRDS(categoricalTreeFilename)

plotTree(harhsalCategoricalTree)
plotTree(commonMainTrees$masterTree)
nodelabels(cex=0.6, col= "green", frame="none")

ZoonomTreeNameToCommon(commonMainTrees$masterTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)

TmmDataOg = read.table("Data/VGP_ALDH1A1_Accession.tsv")
tmmData = as.data.frame(t(TmmDataOg))
names(tmmData) = c("Accession", "ALDH1A1")
tmmData = tmmData[-1,]
write.csv(tmmData, "Data/VGP_ALDH1A1_Accession.csv")
#---------------------------------------------------------------------
# --- overlap counting  --- 
# --------------------------------------------------------------------

geneResults = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysisCombinedGeneResults.rds")
length(which(geneResults$`CH-significant`))
length(which(geneResults$`CH-HI-UnpermOverlap`))
length(which(geneResults$`CH-HV-UnpermOverlap`))

length(which(geneResults$`HI-significant`))
length(which(geneResults$`HV-significant`))
length(which(geneResults$`HI-HV-UnpermOverlap`))


GoResults = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysisCombinedGoResults-KEggReactome.rds")
length(which(GoResults$`CH-significant`))
length(which(GoResults$`CH-HI-UnpermOverlap`))
length(which(GoResults$`CH-HV-UnpermOverlap`))

#---------------------------------------------------------------------
# --- ceteacean counting  --- 
# --------------------------------------------------------------------

mergedData = read.csv("Data/mergedData.csv")

psiciSpecies = which(mergedData$DerekDietClassification90InsVertivoreSorting == "C-InsVertivore-Piscivore")

psiciSpeciesNames = mergedData$ZoonomiaTip[psiciSpecies]


speciesFilter = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysisSpeciesFilter.rds")

length(which(psiciSpeciesNames %in% speciesFilter))


which(as.logical(mergedData$isCricetidae))

cricedidaeVector = readRDS("Output/CladeBinaryCricetidae/CladeBinaryCricetidaeCategoricalPhenotypeVector.rds")
which(cricedidaeVector ==1)

#---------------------------------------------------------------------
# --- looking into geneset distribution  --- 
# --------------------------------------------------------------------

gmtAnnotations = read.gmt("Data/KeggReactome.gmt")


signalingWords = c("Signaling", "Methylation", "Ubiquitination")
metabolismWords = c("Metabolism", "Catabolism", "Anabolism", "Biosynthesis")

genesetNames = gmtAnnotations$geneset.names

metabolismPathways = c()

metabolismPathways = append(metabolismPathways, grep(metabolismWords[1], genesetNames, ignore.case = T))
metabolismPathways = append(metabolismPathways,grep(metabolismWords[2], genesetNames, ignore.case = T))
metabolismPathways = append(metabolismPathways,grep(metabolismWords[3], genesetNames, ignore.case = T))
metabolismPathways = append(metabolismPathways,grep(metabolismWords[4], genesetNames, ignore.case = T))
metabolismPathways = unique(metabolismPathways)

signalingPathways = c()
signalingPathways = append(signalingPathways, grep(signalingWords[1], genesetNames, ignore.case = T))
signalingPathways = append(signalingPathways,grep(signalingWords[2], genesetNames, ignore.case = T))
signalingPathways = append(signalingPathways,grep(signalingWords[3], genesetNames, ignore.case = T))
signalingPathways = unique(signalingPathways)

#---------------------------------------------------------------------
# --- alternate visualization --- 
# --------------------------------------------------------------------



combinedCorrelations = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysiscombinedGOResultsWithAlternates-KeggReactome.rds")
combinedCorrelations = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysiscombinedGeneResultsWithAlternates.rds")

combinedCorrelations = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysisAlternatesCombinedCorrelations.rds")


names(combinedCorrelations)
IndexesToRun = c(7, 1, 4, 5)


for(i in IndexesToRun){
  currentDataset = combinedCorrelations[i]
  currentName = names(currentDataset)
  
  currentDataset = combinedCorrelations[[i]]
  
  names(currentDataset)
  
}





targetCol = "-PadjMedian"

alternateTable = dataframe()

for (p in prefixes) {
  
  ogValCol = paste0(p, "-p.adj")
  sig_col <- paste0(p, "-significant")
  testCol <- paste0(p, targetCol)
  
  which(names(combinedCorrelations) %in% c(ogValCol, sig_col, testCol))
  
  alternateTable = cbind(alternateTable, combinedCorrelations[,which(names(combinedCorrelations) %in% c(ogValCol, sig_col, testCol))])
  

  ggplot()  
  
  
  
}







# initialize result table
resultAlternateTable <- data.frame(
  gene = rownames(df)
)

# loop through prefixes
for (p in prefixes) {
  
  sig_col <- paste0(p, "-significant")
  testCol <- paste0(p, targetCol)
  
  
  
  
  # safety check (in case some prefixes are missing columns)
  if (!all(c(sig_col, testCol) %in% names(df))) next
  
  resultAlternateTable[[paste0(p, "_testCol")]] <-
    ifelse(df[[sig_col]] == TRUE,
           df[[testCol]],
           0)
}
rownames(resultAlternateTable) = resultAlternateTable$gene

importantCols = c(2,4,6,8)

par(mfrow= c(2,2))
for(i in importantCols){
  imporantRows = which(resultAlternateTable[i] > 0 & !is.na(resultAlternateTable[i]))
  hist(resultAlternateTable[imporantRows,c(i)], main = paste0(substr(colnames(resultAlternateTable)[i],0, 2), targetCol, " Num sig:", length(which(resultAlternateTable[imporantRows,c(i)] < 0.05)), " Fraction sig:", (round(length(which(resultAlternateTable[imporantRows,c(i)] < 0.05))/length(imporantRows),3 ))), breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.60, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0), xlim = c(0,1))
}





#---------------------------------------------------------------------
# --- Look into semantic similarity  --- 
# --------------------------------------------------------------------
require(devtools)

install_github("TranslationalBioinformaticsUnit/GeneSetCluster")







#---------------------------------------------------------------------
# --- Make code to split significance results by direction  --- 
# --------------------------------------------------------------------

combinedResults

statColumns = names(combinedResults)[grep("Rho", names(combinedResults))]
GoStatColumns = names(GoCombinedResults)[grep("stat", names(GoCombinedResults))]

geneSignificanceResultsDirection = combinedResults[, names(combinedResults) %in% c(significanceColumns, statColumns)]

GoSignificanceResultsDirection = GoCombinedResults[, names(GoCombinedResults) %in% c(GoSignificanceColumns, GoStatColumns)]



sigDirectionData  =geneSignificanceResultsDirection
sigDirectionData  =GoSignificanceResultsDirection

usedStatCols = GoStatColumns
usedStatCols = statColumns

directionData = GoSignificanceResultsDirection



makeDirectionalResults = function(directionData, usedStatCols){
  sigDirectionDataPositive = directionData
  sigDirectionDataNegative = directionData
  for(i in usedStatCols){
    currentPrefix = substr(i, 1, 2)
    
    # get matching significant column
    statCol = names(directionData)[grep(currentPrefix, names(directionData))][1]
    sigCol = names(directionData)[grep(currentPrefix, names(directionData))][2]
    
    # set significance to NA where Rho > 0
    
    sigDirectionDataPositive[[sigCol]][sigDirectionDataPositive[[statCol]] > 0] <- NA
    
    sigDirectionDataNegative[[sigCol]][sigDirectionDataNegative[[statCol]] < 0] <- NA
    
  } 
  sigDirectionDataNegative = sigDirectionDataNegative[, -which(names(sigDirectionDataNegative) %in% usedStatCols)]
  sigDirectionDataPositive = sigDirectionDataPositive[, -which(names(sigDirectionDataPositive) %in% usedStatCols)]
  directionalResults = list(sigDirectionDataPositive, sigDirectionDataNegative)
  return(directionalResults)
}

testOut = makeDirectionalResults(geneSignificanceResultsDirection, statColumns)

#---------------------------------------------------------------------
# --- Misc results review as part of paper refresh   --- 
# --------------------------------------------------------------------

fullPredatorPhenotypes = readRDS("Output/ComplexDietCentralAnalysisAllSpecies/ComplexDietCentralAnalysisAllSpeciesCategoricalPhenotypeVector.rds")

table(fullPredatorPhenotypes)


mainPhenotypes = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysisCategoricalPhenotypeVector.rds")
length(mainPhenotypes)
table(mainPhenotypes)


mainCorrelations = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysiscombinedGeneResultsWithAlternates.rds")
mainCorrelations[(which(mainCorrelations$`IV-significant`)),]


mainEnrichments = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysiscombinedGOResultsWithAlternates-KeggReactome.rds")
length(which(mainEnrichments$`IV-significant`))

length(which(mainEnrichments$`IV-p.adj` <0.1))


oldEnrichments = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferencecombinedGOResults-KeggReactome.rds")
length(which(oldEnrichments$`IV-p.adj` <0.1))


mainCorrelations[,grep("significant", names(mainCorrelations))]

colSums(mainCorrelations[,grep("significant", names(mainCorrelations))], na.rm = T)

#---------------------------------------------------------------------
# --- Look into continuous results    --- 
# --------------------------------------------------------------------
protienEnrichments = readRDS("Output/ContinousCdcaCp/ContinousCdcaCpEnrichment-KeggReactome.rds")[[1]]



#---------------------------------------------------------------------
# --- look into alterante affecting gene rank    --- 
# --------------------------------------------------------------------
combinedCorrelations = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysiscombinedGeneResults.rds")


readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysisCategoricalTree.rds")
readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysisSpeciesFilter.rds")





#---------------------------------------------------------------------
# --- look into alterante results   --- 
# --------------------------------------------------------------------
combinedCorrelations = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysiscombinedGeneResultsWithAlternates.rds")

combinedCorrelations = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysiscombinedGOResultsWithAlternates-KeggReactome.rds")



prefixes <- unique(sub("-.*", "", colnames(combinedCorrelations)))

for (p in prefixes) {
  mean_col <- paste0(p, "-PadjMean")
  sd_col   <- paste0(p, "-PadjSD")
  new_col  <- paste0(p, "-PadjMeanMinusSD")
  
  # only create if both columns exist
  if (mean_col %in% colnames(combinedCorrelations) & sd_col %in% colnames(combinedCorrelations)) {
    combinedCorrelations[[new_col]] <- combinedCorrelations[[mean_col]] - combinedCorrelations[[sd_col]]
  }
}



combinedCorrelations

str(combinedCorrelations)

df <- combinedCorrelations

# extract prefixes (CH, CO, HI, etc.)
prefixes <- unique(sub("-.*", "", grep("-", names(df), value = TRUE)))

targetCol = "-PadjMean"
targetCol = "-PadjMeanMinusSD"
targetCol = "-PadjMedian"
targetCol = "-PadjNumSignificant"


# initialize result table
resultAlternateTable <- data.frame(
  gene = rownames(df)
)

# loop through prefixes
for (p in prefixes) {
  
  sig_col <- paste0(p, "-significant")
  testCol <- paste0(p, targetCol)
  
  
  # safety check (in case some prefixes are missing columns)
  if (!all(c(sig_col, testCol) %in% names(df))) next
  
  resultAlternateTable[[paste0(p, "_testCol")]] <-
    ifelse(df[[sig_col]] == TRUE,
           df[[testCol]],
           0)

}
rownames(resultAlternateTable) = resultAlternateTable$gene

importantCols = c(2,4,6,8)

sigInfo = data.frame(row.names = c("numTotal", "numSig", "FractionSig", "numSigPlus", "FractionSigPlus"))

par(mfrow= c(2,2))
for(i in importantCols){
  imporantRows = which(resultAlternateTable[i] != 0 & !is.na(resultAlternateTable[i]))
  hist(resultAlternateTable[imporantRows,c(i)], main = paste0(substr(colnames(resultAlternateTable)[i],0, 2), targetCol, " Num sig:", length(which(resultAlternateTable[imporantRows,c(i)] < 0.05)), " Fraction sig:", (round(length(which(resultAlternateTable[imporantRows,c(i)] < 0.05))/length(imporantRows),3 ))), breaks = c(-1, 0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.60, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0), xlim = c(0,1))
  #hist(resultAlternateTable[imporantRows,c(i)], main = paste0(substr(colnames(resultAlternateTable)[i],0, 2), targetCol, " Over 50:", round(length(which(resultAlternateTable[imporantRows,c(i)] > 50))/length(imporantRows),2), " Mean sig:", round((mean(resultAlternateTable[imporantRows,c(i)])), 0)))
  
  
  numTotal = length(imporantRows)
  numSig =  length(which(resultAlternateTable[imporantRows,c(i)] < 0.05))
  FractionSig = length(which(resultAlternateTable[imporantRows,c(i)] < 0.05))/length(imporantRows)
  numSigPlus = length(which(resultAlternateTable[imporantRows,c(i)] < 0.1 & resultAlternateTable[imporantRows,c(i)] > 0.05))
  fractionSigPlus = length(which(resultAlternateTable[imporantRows,c(i)] < 0.1 & resultAlternateTable[imporantRows,c(i)] > 0.05))/length(imporantRows)
  sigInfo = cbind(sigInfo, c(numTotal, numSig, FractionSig, numSigPlus, fractionSigPlus))
}


sum(sigInfo[2,])/sum(sigInfo[1,])
sum(sigInfo[4,])/sum(sigInfo[1,])


importantCols = c(2,4,6,8)
importantCols = c(2,4,6)
outVals = numeric()
par(mfrow= c(2,2))
for(i in importantCols){
  imporantRows = which(resultAlternateTable[i] != 0 & !is.na(resultAlternateTable[i]))
  outVals = append(outVals, resultAlternateTable[imporantRows,c(i)])
  hist(resultAlternateTable[imporantRows,c(i)], main = paste0(substr(colnames(resultAlternateTable)[i],0, 2), targetCol, " Over 50:", round(length(which(resultAlternateTable[imporantRows,c(i)] > 50))/length(imporantRows),2), " Mean sig:", round((mean(resultAlternateTable[imporantRows,c(i)])), 0)))
  
}
sum(outVals)/length(outVals)
length(which(outVals > 50))/length(outVals)

sum(sigInfo[2,])/sum(sigInfo[1,])
sum(sigInfo[4,])/sum(sigInfo[1,])


library(dplyr)
library(ggplot2)
library(ggpmisc)
library(gridExtra)




df = combinedCorrelations

prefixes <- c("CH", "HI", "HV", "IV")

plots <- lapply(prefixes, function(pref) {
  
  x_col <- paste0(pref, "-p.adj")
  y_col <- paste0(pref, "-PadjMedian")
  
  df_clean <- df %>%
    filter(!is.na(.data[[x_col]]), !is.na(.data[[y_col]]))
  
  p_all <- ggplot(df_clean, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", se = TRUE) +
    ggtitle(paste(pref, "- All Points")) +
    stat_poly_eq(
      aes(label = paste(..rr.label.., ..eq.label.., sep = "~~~")),
      formula = y ~ x,
      parse = TRUE
    ) +
    theme_minimal()+ scale_x_log10()+ scale_y_log10()
  
  p_sig <- df_clean %>%
    filter(.data[[x_col]] < 0.05) %>%
    ggplot(aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE) +
    ggtitle(paste(pref, "- p.adj < 0.05")) +
    stat_poly_eq(
      aes(label = paste(..rr.label.., ..eq.label.., sep = "~~~")),
      formula = y ~ x,
      parse = TRUE
    ) +
    theme_minimal()+ scale_x_log10() + scale_y_log10()
  
  list(all = p_all, sig = p_sig)
})

j=4
grid.arrange(plots[[j]][[1]], plots[[j]][[2]], nrow=1)



# ?hist

#plot pvalues 

# initialize result table
resultAlternateTable <- data.frame(
  gene = rownames(df)
)

# loop through prefixes
for (p in prefixes) {
  
  sig_col <- paste0(p, "-significant")
  mainCol <- paste0(p, "-p.adj")
  alternateCol <- paste0(p, "-PadjMean")
  
  
  
  
  # safety check (in case some prefixes are missing columns)
  if (!all(c(sig_col, alternateCol) %in% names(df))) next
  resultAlternateTable[[paste0(p, "mainCol")]] <-
    ifelse(df[[sig_col]] == TRUE,
           df[[mainCol]],
           0)
  resultAlternateTable[[paste0(p, "alternateCol")]] <-
    ifelse(df[[sig_col]] == TRUE,
           df[[alternateCol]],
           0)
}


rownames(resultAlternateTable) = resultAlternateTable$gene

importantCols = c(2,6,10,14)

par(mfrow= c(2,2))
for(i in importantCols){
  
  
  imporantRows = which(resultAlternateTable[i] > 0 & !is.na(resultAlternateTable[i]))
  plot(resultAlternateTable[imporantRows, i], resultAlternateTable[imporantRows,(i+1)])
  
  #hist(resultAlternateTable[imporantRows,c(i)], main = paste(colnames(resultAlternateTable)[i], "Mean:", mean(resultAlternateTable[imporantRows,c(i)])))
}



#OG HIstogram 
# initialize result table
sigResultAlternateTable <- data.frame(
  gene = rownames(df)
)

# loop through prefixes
for (p in prefixes) {
  
  sig_col <- paste0(p, "-significant")
  padj_col <- paste0(p, "-PadjNumSignificant")
  
  
  
  
  # safety check (in case some prefixes are missing columns)
  if (!all(c(sig_col, padj_col) %in% names(df))) next
  
  sigResultAlternateTable[[paste0(p, "_padjSigCount")]] <-
    ifelse(df[[sig_col]] == TRUE,
           df[[padj_col]],
           0)
}

rownames(sigResultAlternateTable) = sigResultAlternateTable$gene

importantCols = c(2,4,6,8)

par(mfrow= c(2,2))
for(i in importantCols){
  imporantRows = which(sigResultAlternateTable[i] > 0 & !is.na(sigResultAlternateTable[i]))
  sigResultAlternateTable[imporantRows,c(1,i)]
  mean(sigResultAlternateTable[imporantRows,c(i)])
  length(sigResultAlternateTable[imporantRows,c(i)] >50)
  length(sigResultAlternateTable[imporantRows,c(i)])
  hist(sigResultAlternateTable[imporantRows,c(i)], main = paste(colnames(sigResultAlternateTable)[i], "Mean:", mean(sigResultAlternateTable[imporantRows,c(i)])))
}


for(i in importantCols){
  imporantRows = which(sigResultAlternateTable[i] > 0 & !is.na(sigResultAlternateTable[i]))
  sigResultAlternateTable[imporantRows,c(1,i)]
  mean(sigResultAlternateTable[imporantRows,c(i)])
  length(sigResultAlternateTable[imporantRows,c(i)] >50)
  length(sigResultAlternateTable[imporantRows,c(i)])
  hist(sigResultAlternateTable[imporantRows,c(i)], main = paste(colnames(sigResultAlternateTable)[i], "Mean:", mean(sigResultAlternateTable[imporantRows,c(i)])))
}



?hist
#---------------------------------------------------------------------
# --- code to add alternates to main    --- 
# --------------------------------------------------------------------
library(purrr)


CorrAlternates = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysisAlternatesCombinedCorrelations.rds")
EnrichAlternates = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysisAlternatesCombinedEnrichments-KeggReactome.rds")


list = EnrichAlternates$`Carnivore-Herbivore`$p.adj
counts <- rowSums(sapply(list, function(x) x < 0.05))


CorrSummaries = map(CorrAlternates, ~ .x[, (ncol(.x) - 4):ncol(.x), drop = FALSE])
EnrichSummaries = map(EnrichAlternates, ~ .x[, (ncol(.x) - 4):ncol(.x), drop = FALSE])



CorSumSingle <- imap(CorrSummaries, function(df, nm) {
  parts <- strsplit(nm, "-")[[1]]
  prefix <- paste0(substr(parts[1], 1, 1),
                   substr(parts[2], 1, 1))
  
  df %>% rename_with(~ paste0(prefix, "-", .x))
}) %>%
  bind_cols()

EnrichSumSingle <- imap(CorrSummaries, function(df, nm) {
  parts <- strsplit(nm, "-")[[1]]
  prefix <- paste0(substr(parts[1], 1, 1),
                   substr(parts[2], 1, 1))
  
  df %>% rename_with(~ paste0(prefix, "-", .x))
}) %>%
  bind_cols()



CorrMain = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysisPairwiseCorrelationFile.rds")


CorrComare = Map(cbind, CorrMain, CorrSummaries)

CorrComare <- lapply(CorrComare, function(df) {
  df$significant <- df$p.adj < 0.05
  df
})

lapply(CorrComare, function(dataframe){
  significantGenes = which(dataframe$significant)
  numSigAlternates = dataframe$PadjNumSignificant[significantGenes]
})


combinedCorrelations = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysiscombinedGeneResults.rds")

combinedCorrelationsWithAlternates = cbind(combinedCorrelations, CorSumSingle)








#---------------------------------------------------------------------
# --- Writring code to combine the alternates together  --- 
# --------------------------------------------------------------------

CorrTest = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysisPairwiseCorrelationFile.rds")
CorrTest2 = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysisPairwiseCorrelationFile.rds")


combined <- lapply(names(runs[[1]]), function(nm) {
  
  # Extract list of data frames for this comparison
  dfs <- lapply(runs, function(x) x[[nm]])
  
  n <- nrow(dfs[[1]])
  
  data.frame(
    Rho   = I(lapply(seq_len(n), function(i) sapply(dfs, function(df) df$Rho[i]))),
    P     = I(lapply(seq_len(n), function(i) sapply(dfs, function(df) df$P[i]))),
    p.adj = I(lapply(seq_len(n), function(i) sapply(dfs, function(df) df$p.adj[i])))
  )
})
names(combined) <- names(runs[[1]])




combined <- lapply(combined, function(df) {
  
  df$P_lt_0.05_count <- sapply(df$P, function(x) sum(x < 0.05, na.rm = TRUE))
  
  df$Padj_lt_0.05_count <- sapply(df$p.adj, function(x) sum(x < 0.05, na.rm = TRUE))
  
  df$Rho_max_diff <- sapply(df$Rho, function(x) {
    if (all(is.na(x))) return(NA_real_)
    max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  })
  
  df$P_max_diff <- sapply(df$P, function(x) {
    if (all(is.na(x))) return(NA_real_)
    max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  })
  
  df$Padj_max_diff <- sapply(df$p.adj, function(x) {
    if (all(is.na(x))) return(NA_real_)
    max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  })
  
  df
})




combined <- lapply(names(runs[[1]]), function(nm) {
  
  dfs <- lapply(runs, function(x) x[[nm]])
  
  list(
    Rho   = do.call(cbind, lapply(dfs, `[[`, "Rho")),
    P     = do.call(cbind, lapply(dfs, `[[`, "P")),
    p.adj = do.call(cbind, lapply(dfs, `[[`, "p.adj"))
  )
})

names(combined) <- names(runs[[1]])



#---------------------------------------------------------------------
# --- Comapring GO results from the histriognathi versions --- 
# --------------------------------------------------------------------
noYeastGOspecies = readRDS("Output/CategoricalInsVertivoreTreeNoYeastLiamInference/CategoricalInsVertivoreTreeNoYeastLiamInferencespeciesFilter.rds")

cat(ZonomNameConvertVectorCommon(noYeastGOspecies, tipColumn = "ZoonomiaTip"), sep = ", ")


noRodentspecies = readRDS("Output/ComplexDietCentralAnalysisNoRodent/ComplexDietCentralAnalysisNoRodentspeciesFilter.rds")

noYeastGOspecies %in% noRodentspecies
noRodentspecies %in% noYeastGOspecies
noRodentspecies[!noRodentspecies %in% noYeastGOspecies]

noRodentPhen = readRDS("Output/ComplexDietCentralAnalysisNoRodent/ComplexDietCentralAnalysisNoRodentcategoricalphenotypeVector.rds")

noRodentPhen[names(noRodentPhen) %in% noRodentspecies[!noRodentspecies %in% noYeastGOspecies]]


paperGOAnalysis = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferencecombinedGOResults-KeggReactome.rds")
noYeastGOAnalysis = readRDS("Output/CategoricalInsVertivoreTreeNoYeastLiamInference/CategoricalInsVertivoreTreeNoYeastLiamInferencecombinedGOResults-KeggReactome.rds")
noRodentGOAnalysis = readRDS("Output/ComplexDietCentralAnalysisNoRodent/ComplexDietCentralAnalysisNoRodentcombinedGOResults-KeggReactome.rds")








paperAnalysis = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferencePairwiseCorrelationFile.rds")
noYeastAnalysis = readRDS("Output/CategoricalInsVertivoreTreeNoYeastLiamInference/CategoricalInsVertivoreTreeNoYeastLiamInferencePairwiseCorrelationFile.rds")

noRodentAnalysis = readRDS("Output/ComplexDietCentralAnalysisNoRodent/ComplexDietCentralAnalysisNoRodentPairwiseCorrelationFile.rds")

paperAnalysis = paperAnalysis[-c(2,3,6,8)]
noYeastAnalysis = noYeastAnalysis[-c(2,3,6,8)]
noRodentAnalysis = noRodentAnalysis[-c(2,3,6,8)]

compareAnalyses = function(ogAnalysis, newAnalysis, prefix = NULL, pCuttof = 0.05){
  
  if(!is.null(prefix)){
    outputOGAnalysis = list()
    outputNewAnalysis = list()
    for(i in 1:length(prefix)){
      currentCols = grep(prefix[i], names(ogAnalysis))
      currentCols = currentCols[1:3]
      
      tempOGAnalysis = ogAnalysis[,currentCols]
      tempNewAnalysis = newAnalysis[,currentCols]
      
      tempOGAnalysis = list(tempOGAnalysis)
      names(tempOGAnalysis) = prefix[i]
      
      tempNewAnalysis = list(tempNewAnalysis)
      names(tempNewAnalysis) = prefix[i]
      
      
      outputOGAnalysis = append(outputOGAnalysis, tempOGAnalysis)
      outputNewAnalysis = append(outputNewAnalysis, tempNewAnalysis)
    }
    
    newAnalysis = outputNewAnalysis
    ogAnalysis = outputOGAnalysis
  }
  
  
  
  results <- tibble(
    index = integer(),
    correlation = numeric(),
    rhoRankcorrelation = numeric(),
    numMismatchedNAs = integer(),
    pCorrelation = numeric(),
    numMismatchedPs = integer(),
    padjCorrelation = numeric(),
    numMismatchedPadjs = integer(),
    totalOgSigPadjs = integer(), 
    totalNewSigPadjs = integer(),
    missingFraction = integer(),
    mismatchedNAs = I(list()),
    mismatchedPs = I(list()),
    mismatchedPadjs = I(list())
  )
  
  for (i in seq_along(ogAnalysis)) {
    
    correlation <- cor(ogAnalysis[[i]][[1]], newAnalysis[[i]][[1]], use = "complete.obs")
    
    mismiatchedNAs <- which(!is.na(ogAnalysis[[i]][[1]]) %in% is.na(newAnalysis[[i]][[1]]))
    
    ogAnalysis[[i]] = ogAnalysis[[i]] %>% mutate(rhoRank = rank(ogAnalysis[[i]][[1]])) 
    newAnalysis[[i]] = newAnalysis[[i]] %>% mutate(rhoRank = rank(newAnalysis[[i]][[1]]))
    
    rhoRankcorrelation <- cor(ogAnalysis[[i]][[4]], newAnalysis[[i]][[4]], use = "complete.obs")
    
    pCorrelation <- cor(ogAnalysis[[i]][[2]], newAnalysis[[i]][[2]], use = "complete.obs")
    
    mismiatchedPs <- which(ogAnalysis[[i]][[2]] < pCuttof)[which(!which(ogAnalysis[[i]][[2]] < pCuttof) %in% 
                                                                   which(newAnalysis[[i]][[2]] < pCuttof))]
    
    padjCorrelation <- cor(ogAnalysis[[i]][[3]], newAnalysis[[i]][[3]], use = "complete.obs")
    
    mismiatchedPadjs <- which(ogAnalysis[[i]][[3]] < pCuttof)[which(!which(ogAnalysis[[i]][[3]] < pCuttof) %in% 
                                                                      which(newAnalysis[[i]][[3]] < pCuttof))]
    
    totalOgSigPadjs = length(which(ogAnalysis[[i]][[3]] < pCuttof))
    totalNewSigPadjs = length(which(newAnalysis[[i]][[3]] < pCuttof))
    
    missingFraction = length(mismiatchedPadjs) / length(which(ogAnalysis[[i]][[3]] < pCuttof))
    
    results <- rbind(results, tibble(
      index = i,
      correlation = correlation,
      rhoRankcorrelation = rhoRankcorrelation,
      numMismatchedNAs = length(mismiatchedNAs),
      pCorrelation = pCorrelation,
      numMissingPs = length(mismiatchedPs),
      padjCorrelation = padjCorrelation,
      numMissingPadjs = length(mismiatchedPadjs),
      totalOgSigPadjs = totalOgSigPadjs,
      totalNewSigPadjs = totalNewSigPadjs,
      missingFraction = missingFraction,
      mismatchedNAs = list(mismiatchedNAs),
      missingPs = list(mismiatchedPs),
      missingPadjs = list(mismiatchedPadjs),
    ))
  }
  row.names(results) = names(ogAnalysis)
  return(results)
}


length(which(noYeastGOAnalysis$`HI-significant`))
length(which(noRodentGOAnalysis$`HI-significant`))




newAnalysis = fullTreeAnalysis
newGOAnalysis = fullTreeGOAnalysis



geneMissingInNewAnalysis = compareAnalyses(noYeastAnalysis, noRodentAnalysis)
geneMissingInOldAnalysis = compareAnalyses(noRodentAnalysis, noYeastAnalysis)


noRodentAnalysis[[3]]$Rho
noYeastAnalysis[[3]]$Rho


rhoDifference = noRodentAnalysis[[3]]$Rho - noYeastAnalysis[[3]]$Rho

which(rhoDifference == max(rhoDifference, na.rm = T))

rownames(noYeastAnalysis[[3]])[which(rhoDifference == max(rhoDifference, na.rm = T))]
maxDifGene = rownames(noRodentAnalysis[[3]])[which(rhoDifference == max(rhoDifference, na.rm = T))]


noRodentRER =  readRDS("Output/ComplexDietCentralAnalysisNoRodent/ComplexDietCentralAnalysisNoRodentRERFile.rds")
noYeastRER =  readRDS("Output/CategoricalInsVertivoreTreeNoYeastLiamInference/CategoricalInsVertivoreTreeNoYeastLiamInferenceRERFile.rds")
noRodentPaths =  readRDS("Output/ComplexDietCentralAnalysisNoRodent/ComplexDietCentralAnalysisNoRodentCategoricalPathsFile.rds")
noYeastPaths =  readRDS("Output/CategoricalInsVertivoreTreeNoYeastLiamInference/CategoricalInsVertivoreTreeNoYeastLiamInferenceCategoricalPathsFile.rds")


plotRers(noRodentRER, maxDifGene, phenv = noRodentPaths)
plotRers(noYeastRER, maxDifGene, phenv = noYeastPaths)




source("Src/Loc/Dev/DisplayCategoricalRERTree.R")

pdf(height = 10, width = 15)
par(mfrow=c(1,2))
displayCategoricalRERTree(mainTrees, noRodentRER, maxDifGene, noRodentPaths, subsetTree = T)
displayCategoricalRERTree(mainTrees, noYeastRER, maxDifGene, noYeastPaths, subsetTree = T)
dev.off()

displayCategoricalRERTree(mainTrees, noRodentRER, maxDifGene, noRodentPaths, subsetTree = F)



goMissingInNewAnalysis = compareAnalyses(noYeastGOAnalysis, noRodentGOAnalysis, c("HI", "HV", "IV", "CH"), pCuttof = 0.1)
goMissingInOldAnalysis = compareAnalyses(noRodentGOAnalysis, noYeastGOAnalysis, c("HI", "HV", "IV", "CH"), pCuttof = 0.1)

goMissingInNewAnalysis = compareAnalyses(noYeastGOAnalysis, noRodentGOAnalysis, c("IV"), pCuttof = 0.1)
goMissingInOldAnalysis = compareAnalyses(noRodentGOAnalysis, noYeastGOAnalysis, c("IV"), pCuttof = 0.1)



#---------------------------------------------------------------------
# --- Looking into Hystriognathi results  --- 
# --------------------------------------------------------------------

histiCorrleaitons = readRDS("Output/CladeBinaryHystricognathi/0-1/CladeBinaryHystricognathi0-1CorrelationFile.rds")

histiCorrleaitons[order(histiCorrleaitons$p.adj),]


length(which(histiCorrleaitons$p.adj < 0.05))

histiEnrich = readRDS("Output/CladeBinaryHystricognathi/0-1/CladeBinaryHystricognathi0-1Enrichment-KeggReactome.rds")
histiEnrich = histiEnrich[[1]]

which(histiEnrich$p.adj < 0.05)


# --- -other binaries

clade= "Hystricognathi"
clade= "Bovidae"
clade= "Cricetidae"
clade= "Vespertilionidae"
clade= "Peropdidae"


binaryCorrelations = readRDS(paste0("Output/CladeBinary",clade,"/0-1/CladeBinary", clade, "0-1CorrelationFile.rds"))
binaryEnrich = readRDS(paste0("Output/CladeBinary",clade,"/0-1/CladeBinary", clade, "0-1Enrichment-KeggReactome.rds"))
binaryEnrich = binaryEnrich[[1]]


length(which(binaryCorrelations$p.adj < 0.05))
which(binaryEnrich$p.adj < 0.05)
length(which(binaryEnrich$p.adj < 0.05))


#binaryCorrelations[order(binaryCorrelations$p.adj),]


#---------------------------------------------------------------------
# --- Setting up clade listing in mergeData  --- 
# --------------------------------------------------------------------

inMergeData = read.csv("Data/mergedData.csv")

inMergeData$isBovidae = rep(0)
inMergeData$isCricetidae = rep(0)
inMergeData$isCervidae = rep(0)
inMergeData$isHystricognathi = rep(0)
inMergeData$isVespertilionidae = rep(0)
inMergeData$isPeropdidae = rep(0)
inMergeData$isInAllSpeciesAnalysis = rep(0)
inMergeData$isInAnalysisWithFullFamilies = rep(0)


bovidaeList = c("vs_HLoryGaz1", "vs_HLbeaHun1", "vs_HLkobLecLec1", "vs_HLkobLecLec1", "vs_HLmadKir1", "vs_HLneoPyg1", "vs_HLphiMax1", "vs_HLoreOre1", "vs_HLneoMos1", "vs_HLaepMel1", "vs_HLtraImb1",  "vs_bisBis1", "vs_HLoviNivLyd1", "vs_HLproPrz1", "Bovidae", "13Herb")

CricetidaeList = c("vs_HLellTal1", "vs_HLellLut1", "vs_HLarvAmp1","vs_HLmicAgr2", "vs_HLmyoGla2", "vs_HLondZib1", "voleClade", "6Herb")
CervidaeList = c("vs_HLhydIne1", "vs_HLmunMun1", "vs_HLodoHem1", "vs_HLantAme1", "vs_HLgirCam1", "Cervidae", "5Herb")

HystricognathiList = c("vs_HLhysCri1", "vs_HLthrSwi1", "vs_HLpetTyp1", "vs_hetGla2", "vs_chiLan1", "vs_HLdinBra1", "vs_HLcteSoc1", "vs_octDeg1", "vs_HLcoePre1", "vs_HLdasPun1", "vs_HLdolPat1", "vs_HLmyoCoy1", "vs_HLhydHyd1", "vs_HLcavTsc1", "gundiGuineaPigClade", "14Herb")

VespertilionidaeList = c("vs_HLmurAurFea1", "Murina", "vs_HLmyoLuc1", "Nearctic", "vs_myoDav1", "vs_HLmyoMyo6", "vs_HLmyoSep1", "vs_HLmyoLuc1", "Myotis", "vs_HLpipPip1", "vs_HLlasBor1", "vs_HLnycHum2", "vs_eptFus1", "Vespertilioninae","Vespertilionidae", "10Inse")
Peropodidae = c("vs_HLmacSob1", "vs_HLpteGig1", "FoxLongTounge", "vs_HLeidHel2", "vs_HLcynBra1", "outerPeropodidae", "vs_HLeonSpe1", "vs_HLrouLes1", "Roussetinae", "Peropodidae", "6Herb")

complexDietAllSpecies = readRDS("Output/ComplexDietCentralAnalysisAllSpecies/ComplexDietCentralAnalysisAllSpeciesSpeciesFilter.rds")
complexDietFullFamily = readRDS("Output/ComplexDietCentralAnalysisNoFamilyPrune/ComplexDietCentralAnalysisNoFamilyPruneSpeciesFilter.rds")




inMergeData[inMergeData$ZoonomiaTip %in% complexDietAllSpecies,]$isInAllSpeciesAnalysis =1
inMergeData[inMergeData$ZoonomiaTip %in% complexDietFullFamily,]$isInAnalysisWithFullFamilies =1
inMergeData[inMergeData$ZoonomiaTip %in% bovidaeList,]$isBovidae =1
inMergeData[inMergeData$ZoonomiaTip %in% CricetidaeList,]$isCricetidae =1
inMergeData[inMergeData$ZoonomiaTip %in% CervidaeList,]$isCervidae =1
inMergeData[inMergeData$ZoonomiaTip %in% HystricognathiList,]$isHystricognathi =1
inMergeData[inMergeData$ZoonomiaTip %in% VespertilionidaeList,]$isVespertilionidae =1
inMergeData[inMergeData$ZoonomiaTip %in% Peropodidae,]$isPeropdidae =1


inMergeData[inMergeData$isBovidae ==1,]


test3 = inMergeData[inMergeData$isHystricognathi ==1,]

cat(test3$MSWC_Family)

#write.csv(inMergeData, "Data/mergedData.csv")

inMergeData[inMergeData$ZoonomiaTip %in% bovidaeList,]$MSWC_Family


test1 = inMergeData[inMergeData$ZoonomiaTip %in% bovidaeList,]
test2 = inMergeData[inMergeData$MSWC_Family == "Bovidae" & !is.na(inMergeData$MSWC_Family) & inMergeData$isInAllSpeciesAnalysis,]

test2 = inMergeData[inMergeData$MSWC_Family == "Bovidae" & !is.na(inMergeData$MSWC_Family) & inMergeData$isInMainAnalysis,]



testTree = categoricalTree
annotColumn = "DerekDietClassification90InsVertivoreSorting"

all.equal(test1, test2)


noYeastSpecies = readRDS("Output/CategoricalInsVertivoreTreeNoYeastLiamInference/CategoricalInsVertivoreTreeNoYeastLiamInferenceSpeciesFilter.rds")
mainAnalysisSpecies = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysisSpeciesFilter.rds")
mainAnalysisPhenotypes = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysiscategoricalphenotypevector.rds")



which(!mainAnalysisSpecies %in% noYeastSpecies)
which(!noYeastSpecies %in% mainAnalysisSpecies)

mainAnalysisSpecies[which(!mainAnalysisSpecies %in% noYeastSpecies)]


mainAnalysisPhenotypes[names(mainAnalysisPhenotypes) %in% mainAnalysisSpecies[which(!mainAnalysisSpecies %in% noYeastSpecies)]]

noYeastTree = readRDS("Output/CategoricalInsVertivoreTreeNoYeastLiamInference/CategoricalInsVertivoreTreeNoYeastLiamInferenceCategoricalTree.rds")
mainAnalysisTree = readRDS("Output/ComplexDietCentralAnalysis/ComplexDietCentralAnalysisCategoricalTree.rds")

all.equal(noYeastTree, mainAnalysisTree)

#---------------------------------------------------------------------
# --- Looking into categorical pruning  --- 
# --------------------------------------------------------------------





alternateSets = readRDS(paste0(outputFolderName, filePrefix, "AlternatePruningSpecies.rds"))
which(names(currentCommonPhenotypeVector) %in% "Bushbaby")

currentCommonPhenotypeVector[77]


which(currentCommonCategoricalTree$tip.label %in% "Bushbaby")

which(currentCommonCategoricalTree$edge[,2] == 16)

currentCommonCategoricalTree$edge.length[56]
currentCommonCategoricalTree$edge.length[57]
currentCommonCategoricalTree$edge.length[55]
currentCommonCategoricalTree$edge[57,]
currentCommonCategoricalTree$edge[56,]
currentCommonCategoricalTree$edge[55,]


which(currentCommonCategoricalTree$edge[,1] == 225)
currentCommonCategoricalTree$edge[58,]
currentCommonCategoricalTree$edge[133,]


allSpeciesTree = readRDS("Output/DuplicatePredatorFullTree/DuplicatePredatorFullTreeCategoricalTree.rds")

testPhenVec = readRDS("Output/DuplicatePredatorFullTree/DuplicatePredatorFullTreeCategoricalPhenotypeVector.rds")
write.csv(phenotypeVectorSaving, "Test.csv")

phenotypeVector = testPhenVec

?drop.tip
UseMethod

getS3method("drop.tip", "phylo")

fullTree = readRDS("Output/PredatorFullTree/PredatorFullTreeCategoricalTree.rds")
dupFullTree = readRDS("Output/DuplicatePredatorFullTree/DuplicatePredatorFullTreeCategoricalTree.rds")

fullPhenVec = readRDS("Output/PredatorFullTree/PredatorFullTreeCategoricalPhenotypeVector.rds")
ogPhenVec = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCategoricalPhenotypeVector.rds")
ogTips = names(ogPhenVec)
dropTips = fullTree$tip.label[!fullTree$tip.label %in% ogTips]

phy = fullTree
tip = dropTips
trim.internal = TRUE
subtree = FALSE
root.edge = 0
rooted = is.rooted(phy)
collapse.singles = TRUE
interactive = FALSE

phy = categoricalDropTip(fullTree, dropTips)

tree$edge[which(tree$edge[,2] == 510),]
tree$edge[which(tree$edge[,1] == 510),]

tree$edge[which(tree$edge[,1] == 501),]
tree$edge[which(tree$edge[,2] == 501),]

tree$edge[which(tree$edge[,1] == 510),]
tree$edge[which(tree$edge[,2] == 510),]

table(e1)

ii[1]

tree$edge[502,]
tree$edge[501,]

which(e1 == 207)
which(e2 == 510)

str(fullTree)

mainTrees = readRDS('data/zoonomiaAllMammalsTrees.rds')

palette(c( "darkgreen", "darkblue","black", "red"))


source("Src/Reu/ZoonomTreeNameToCommon.R")
nameColumn = "ZoonomiaTip"

commonMainTrees = mainTrees
commonMainTrees$masterTree = ZoonomTreeNameToCommon(commonMainTrees$masterTree, tipCol = nameColumn)
commonFullTree = ZoonomTreeNameToCommon(fullTree, tipCol = nameColumn)
commonPhy = ZoonomTreeNameToCommon(phy, tipCol = nameColumn)





commonCategoricalTree$tip.label
commonMainTrees$masterTree$tip.label

fullTree$tip.label
dupFullTree$tip.label
mainTrees$masterTree$tip.label
phy$tip.label






commonPhyFlipped = commonPhy

ntips = length(commonPhyFlipped$tip.label)
commonPhyFlipped$edge[c(commonPhyFlipped$edge[,2] < ntips), ]

phy$edge[,2]

tipVal = "vs_monDom5"
which(fullTree$tip.label == tipVal)
which(fullTree$edge[,2] == which(fullTree$tip.label == tipVal))
fullTree$edge.length[which(fullTree$edge[,2] == which(fullTree$tip.label == tipVal))]


which(phy$tip.label == tipVal)
which(phy$edge[,2] == which(phy$tip.label == tipVal))
phy$edge.length[which(phy$edge[,2] == which(phy$tip.label == tipVal))]

tipVal = "Opossum"
which(commonPhy$tip.label == tipVal)
which(commonPhy$edge[,2] == which(commonPhy$tip.label == tipVal))
commonPhy$edge.length[which(commonPhy$edge[,2] == which(commonPhy$tip.label == tipVal))]




which(commonPhy$tip.label == "Opossum")
which(commonMainTrees$masterTree$tip.label == "Opossum")
which(commonFullTree$tip.label == "Opossum")


length(commonFullTree$tip.label)
mainTrees$masterTree
commonMainTrees$masterTree


tree = commonPhy
master = commonMainTrees$masterTree


cm = intersect(master$tip.label, tree$tip.label)
tipsToDrop = master$tip.label[!master$tip.label %in% cm]
master1 = drop.tip(master, tipsToDrop)
master2 = pruneTree(master, cm)


cm2 = intersect(master$tip.label, commonFullTree$tip.label)
master3 = pruneTree(master, cm2)
cm3 = intersect(master3$tip.label, tree$tip.label)
master4 = master3 = pruneTree(master3, cm3)

all.equal(master1, master2)

tipVal = "Opossum"
which(tree$tip.label == tipVal)
which(tree$edge[,2] == which(tree$tip.label == tipVal))

which(master$tip.label == tipVal)
which(master$edge[,2] == which(master$tip.label == tipVal))

which(master1$tip.label == tipVal)
which(master1$edge[,2] == which(master1$tip.label == tipVal))

which(master4$tip.label == tipVal)
which(master4$edge[,2] == which(master1$tip.label == tipVal))

length(master4$edge[,2])
length(tree$edge[,2])


tab2 = tab <- tabulate(e2)






pdf()
plotTree(commonPhy)
nodelabels(cex=0.6, col= "green", frame="none")
tiplabels(cex=0.6, col= "green", frame="none")
plotTreeCategorical(commonPhy, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = commonMainTrees$masterTree)
nodelabels(cex=0.6, col= "green", frame="none")
tiplabels(cex=0.6, col= "green", frame="none")
plotTreeCategorical(phy, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = mainTrees$masterTree)
nodelabels(cex=0.6, col= "green", frame="none")
tiplabels(cex=0.6, col= "green", frame="none")
dev.off()


ogTree = ogPhenVec = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCategoricalTree.rds")
oldPhy = phy

oldPhy = categoricalDropTip(fullTree, dropTips)
phy = categoricalDropTip(dupFullTree, dropTips)

all.equal(phy$tip.label, ogTree$tip.label)
all.equal(phy$edge.length, ogTree$edge.length)
all.equal(phy$edge[,1], ogTree$edge[,1])
all.equal(phy$edge[,2], ogTree$edge[,2])
phy$edge.length == ogTree$edge.length

disagreeingBranches = which(!phy$edge.length == ogTree$edge.length)

commonPhy = ZoonomTreeNameToCommon(phy, tipCol = nameColumn)


pdf(width = 15)  
par(mfrow = c(1,2))
plotTreeCategorical(commonPhy, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = commonMainTrees$masterTree)
edgelabels(edge = disagreeingBranches, cex=0.4, col= "purple", frame="none")
plotTreeCategorical(ogTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = mainTrees$masterTree)
edgelabels(edge = disagreeingBranches, cex=0.4, col= "purple", frame="none")

dev.off()


disagreements = data.frame(disagreeingBranches, phy$edge.length[disagreeingBranches], ogTree$edge.length[disagreeingBranches])
names(disagreements) = c("branch", "newPhen", "OgPhen")
disagreements



phy$edge[98,]

phy$edge[which(phy$edge[,1] == 243),]



phy=dupFullTree

pdf()
plotTreeCategorical(phyUnmerged, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = mainTrees$masterTree)
nodelabels(cex=0.2, col= "purple", frame="none")
dev.off()
singles


closestSpecies = which(phyUnmerged$tip.label == "vs_mesAur1")

edgeFrame$index
edgeFrame[98,]


phyUnmerged$edge[which(phyUnmerged$edge[,1] == 269),]
phyUnmerged$edge[which(phyUnmerged$edge[,1] == 278),]
phyUnmerged$edge[which(phyUnmerged$edge[,1] == 279),]
phyUnmerged$edge[which(phyUnmerged$edge[,1] == 280),]

all.equal(tree$tip.label, ogTree$tip.label)
all.equal(tree$edge.length, ogTree$edge.length)
all.equal(tree$edge[,1], ogTree$edge[,1])
all.equal(tree$edge[,2], ogTree$edge[,2])
tree$edge.length == ogTree$edge.length

disagreeingBranches = which(!tree$edge.length == ogTree$edge.length)


tree$edge[98,]
tree$edge.length[98]
ogTree$edge.length[98]

firstBranch = which(phyUnmerged$edge[,2] == closestSpecies)
firstParent = phyUnmerged$edge[firstBranch,1]

firstParent %in% singles 

which(phyUnmerged$edge[,1] == 277)


phyUnmerged$edge[which(phyUnmerged$edge[,2] == 255),]
phyUnmerged$edge[which(phyUnmerged$edge[,2] == 254),]

phyUnmerged$edge[which(phyUnmerged$edge[,1] == phyUnmerged$edge[,2])]

phyUnmerged$edge.length[c(82,83)]

problemBranch = 98
problemDaughter = phy$edge[problemBranch,2]

problemii = which(ii == problemDaughter)
problemjj = jj[problemii]
phyUnmerged$edge[problemjj,]

which(phyUnmerged$edge[,1] == 239)


all.equal(e1, edgeFrame$e1)

edgeFrame[1,]

singles
jj[120]
ii[120]

phyUnmerged$edge[15,]
phyUnmerged$edge[16,]

phyUnmerged$edge.length[16]
phyUnmerged$tip.label[2]

fullPhenVec[names(fullPhenVec) %in% phyUnmerged$tip.label[2]]

#Okay, that one's really weird -- it seems to be a bug in the originalPredatorFullTree? Because it's not there on the duplicate. 

fullTree = readRDS("Output/PredatorFullTree/PredatorFullTreeCategoricalTree.rds")
dupFullTree = readRDS("Output/DuplicatePredatorFullTree/DuplicatePredatorFullTreeCategoricalTree.rds")

all.equal(fullTree, dupFullTree)

all.equal(fullTree$tip.label, dupFullTree$tip.label)
all.equal(fullTree$edge.length, dupFullTree$edge.length)
all.equal(fullTree$edge[,1], dupFullTree$edge[,1])
all.equal(fullTree$edge[,2], dupFullTree$edge[,2])
fullTree$edge.length == dupFullTree$edge.length

disagreeingBranchesFull = which(!fullTree$edge.length == dupFullTree$edge.length)

pdf(width = 15)  
par(mfrow = c(1,2))
plotTreeCategorical(fullTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = mainTrees$masterTree)
edgelabels(edge = disagreeingBranchesFull, cex=0.4, col= "purple", frame="none")
plotTreeCategorical(dupFullTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = mainTrees$masterTree)
edgelabels(edge = disagreeingBranchesFull, cex=0.4, col= "purple", frame="none")

dev.off()


all.equal(args1[2:10], args2[2:10])

#those in thoery SHOULD Be identical, because the args are identical, but they... aren't. Which is very odd. 


beforeCollapeEdgeLengths = phy$edge.length


phy$edge.length

unique(table(phy$edge[,2]))

collapse.singles


rm(collapse.singles)

oldTree = tree
all.equal(oldTree, tree)

table(e1)

#---------------------------------------------------------------------
# --- Making cytoscape plots  --- 
# --------------------------------------------------------------------


GOOutput = readRDS("Output/CategoricalInsVertivoreTreeNoYeastLiamInference/Carnivore-Herbivore/CategoricalInsVertivoreTreeNoYeastLiamInferenceCarnivore-HerbivoreEnrichment-KeggReactome.rds")

which(GoSignificanceResults$`HV-significant` & GoSignificanceResults$`HI-significant` &! GoSignificanceResults$`CH-significant`)
GoSignificanceResults[635,]


#Looking into overrepresented family ratios of total mammals

mammalDatabase = read.csv("../../MiscData/MDD/MDD/Species_Syn_v2.4.csv")

familySizes = table(mammalDatabase$MDD_family)

bigFamilies = familySizes[names(familySizes) %in% c("Bovidae", "Cervidae", "Pteropodidae", "Cercopithecidae", "Cricetidae")]

sum(bigFamilies)/sum(familySizes)

FullHerbivoreData = table(familyByDiet$MSWC_Family[familyByDiet$diet == "Herbivore"])
HerbvioreFamilies = familyByDiet[familyByDiet$MSWC_Family %in% names(FullHerbivoreData),]
HerbvioreFamiliySize = table(HerbvioreFamilies$MSWC_Family)

as.integer(FullHerbivoreData)/as.integer(HerbvioreFamiliySize)

herbivoreRatio = (FullHerbivoreData/HerbvioreFamiliySize)

FullHerbivoreData[names(FullHerbivoreData) %in% c("Megalonychidae",  "Myocastoridae")]

names(FullHerbivoreData)

herbFamilies = familySizes[names(familySizes) %in% names(FullHerbivoreData)]

herbFamilies * herbivoreRatio
names(herbFamilies) %in% names(herbivoreRatio)
herbivoreRatio[!names(herbivoreRatio) %in% names(herbFamilies)]

herbivoreRatio = herbivoreRatio[names(herbivoreRatio) %in% names(herbFamilies)]

all.equal(names(herbivoreRatio), names(herbFamilies))

length(herbFamilies)
length(herbivoreRatio)

estimatedHerbivores = herbFamilies * herbivoreRatio
estimatedHerbivores = estimatedHerbivores + 2 # Adding back the species removed ddue ot the family mismatch

sum(bigFamilies)/ sum(estimatedHerbivores)





table(HerbvioreFamilies$MSWC_Family)
table(FullHerbivoreData)

for(i in 1:length(herbFamilies)){
  currentFamily = herbFamilies[i]
  
  
}

#---------------------------------------------------------------------
# --- Comapring GO results from various pehnotype versions --- 
# --------------------------------------------------------------------
paperGOAnalysis = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferencecombinedGOResults-KeggReactome.rds")
familyAgnosticGOAnalysis =readRDS("Output/CategoricalInsVertivoreTreeFamilyAgnostictLiamInference/CategoricalInsVertivoreTreeFamilyAgnostictLiamInferencecombinedGOResults-KeggReactome.rds")
noYeastGOAnalysis = readRDS("Output/CategoricalInsVertivoreTreeNoYeastLiamInference/CategoricalInsVertivoreTreeNoYeastLiamInferencecombinedGOResults-KeggReactome.rds")
noManualGOAnalysis = readRDS("Output/CategoricalInsVertivoreTreeNoManualLiamInference/CategoricalInsVertivoreTreeNoManualLiamInferencecombinedGOResults-KeggReactome.rds")
fullTreeGOAnalysis = readRDS("Output/PredatorFullTree/PredatorFullTreecombinedGOResults-KeggReactome.rds")
noRodentGOAnalysis = readRDS("Output/ComplexDietCentralAnalysisNoRodent/ComplexDietCentralAnalysisNoRodentcombinedGOResults-KeggReactome.rds")


ComplexDietCentralAnalysisNoRodent


paperAnalysis = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferencePairwiseCorrelationFile.rds")
noYeastAnalysis = readRDS("Output/CategoricalInsVertivoreTreeNoYeastLiamInference/CategoricalInsVertivoreTreeNoYeastLiamInferencePairwiseCorrelationFile.rds")
noManualAnalysis = readRDS("Output/CategoricalInsVertivoreTreeNoManualLiamInference/CategoricalInsVertivoreTreeNoManualLiamInferencePairwiseCorrelationFile.rds")
familyAgnosticAnalysis =readRDS("Output/CategoricalInsVertivoreTreeFamilyAgnostictLiamInference/CategoricalInsVertivoreTreeFamilyAgnostictLiamInferencePairwiseCorrelationFile.rds")
fullTreeAnalysis = readRDS("Output/PredatorFullTree/PredatorFullTreePairwiseCorrelationFile.rds")
noRodentAnalysis = readRDS("Output/ComplexDietCentralAnalysisNoRodent/ComplexDietCentralAnalysisNoRodentPairwiseCorrelationFile.rds")

paperAnalysis = paperAnalysis[-c(2,3,6,8)]
noYeastAnalysis = noYeastAnalysis[-c(2,3,6,8)]
noManualAnalysis = noManualAnalysis[-c(2,3,6,8)]
familyAgnosticAnalysis = familyAgnosticAnalysis[-c(2,3,6,8)]
fullTreeAnalysis = fullTreeAnalysis[-c(2,3,6,8)]
noRodentAnalysis = noRodentAnalysis[-c(2,3,6,8)]

compareAnalyses = function(ogAnalysis, newAnalysis, prefix = NULL, pCuttof = 0.05){
  
  if(!is.null(prefix)){
    outputOGAnalysis = list()
    outputNewAnalysis = list()
    for(i in 1:length(prefix)){
      currentCols = grep(prefix[i], names(ogAnalysis))
      currentCols = currentCols[1:3]
      
      tempOGAnalysis = ogAnalysis[,currentCols]
      tempNewAnalysis = newAnalysis[,currentCols]
      
      tempOGAnalysis = list(tempOGAnalysis)
      names(tempOGAnalysis) = prefix[i]
      
      tempNewAnalysis = list(tempNewAnalysis)
      names(tempNewAnalysis) = prefix[i]
      
      
      outputOGAnalysis = append(outputOGAnalysis, tempOGAnalysis)
      outputNewAnalysis = append(outputNewAnalysis, tempNewAnalysis)
    }
    
    newAnalysis = outputNewAnalysis
    ogAnalysis = outputOGAnalysis
  }
  
  
  
  results <- tibble(
    index = integer(),
    correlation = numeric(),
    rhoRankcorrelation = numeric(),
    numMismatchedNAs = integer(),
    pCorrelation = numeric(),
    numMismatchedPs = integer(),
    padjCorrelation = numeric(),
    numMismatchedPadjs = integer(),
    totalOgSigPadjs = integer(), 
    totalNewSigPadjs = integer(),
    missingFraction = integer(),
    mismatchedNAs = I(list()),
    mismatchedPs = I(list()),
    mismatchedPadjs = I(list())
  )
  
  for (i in seq_along(ogAnalysis)) {
    
    correlation <- cor(ogAnalysis[[i]][[1]], newAnalysis[[i]][[1]], use = "complete.obs")
    
    mismiatchedNAs <- which(!is.na(ogAnalysis[[i]][[1]]) %in% is.na(newAnalysis[[i]][[1]]))
    
    ogAnalysis[[i]] = ogAnalysis[[i]] %>% mutate(rhoRank = rank(ogAnalysis[[i]][[1]])) 
    newAnalysis[[i]] = newAnalysis[[i]] %>% mutate(rhoRank = rank(newAnalysis[[i]][[1]]))
    
    rhoRankcorrelation <- cor(ogAnalysis[[i]][[4]], newAnalysis[[i]][[4]], use = "complete.obs")
    
    pCorrelation <- cor(ogAnalysis[[i]][[2]], newAnalysis[[i]][[2]], use = "complete.obs")
    
    mismiatchedPs <- which(ogAnalysis[[i]][[2]] < pCuttof)[which(!which(ogAnalysis[[i]][[2]] < pCuttof) %in% 
                                                                which(newAnalysis[[i]][[2]] < pCuttof))]
    
    padjCorrelation <- cor(ogAnalysis[[i]][[3]], newAnalysis[[i]][[3]], use = "complete.obs")
    
    mismiatchedPadjs <- which(ogAnalysis[[i]][[3]] < pCuttof)[which(!which(ogAnalysis[[i]][[3]] < pCuttof) %in% 
                                                                   which(newAnalysis[[i]][[3]] < pCuttof))]
    
    totalOgSigPadjs = length(which(ogAnalysis[[i]][[3]] < pCuttof))
    totalNewSigPadjs = length(which(newAnalysis[[i]][[3]] < pCuttof))
    
    missingFraction = length(mismiatchedPadjs) / length(which(ogAnalysis[[i]][[3]] < pCuttof))
    
    results <- rbind(results, tibble(
      index = i,
      correlation = correlation,
      rhoRankcorrelation = rhoRankcorrelation,
      numMismatchedNAs = length(mismiatchedNAs),
      pCorrelation = pCorrelation,
      numMissingPs = length(mismiatchedPs),
      padjCorrelation = padjCorrelation,
      numMissingPadjs = length(mismiatchedPadjs),
      totalOgSigPadjs = totalOgSigPadjs,
      totalNewSigPadjs = totalNewSigPadjs,
      missingFraction = missingFraction,
      mismatchedNAs = list(mismiatchedNAs),
      missingPs = list(mismiatchedPs),
      missingPadjs = list(mismiatchedPadjs),
    ))
  }
  row.names(results) = names(ogAnalysis)
  return(results)
}




length(which(noYeastGOAnalysis$`HI-significant`))
length(which(noRodentGOAnalysis$`HI-significant`))




newAnalysis = fullTreeAnalysis
newGOAnalysis = fullTreeGOAnalysis




goMissingInNewAnalysis = compareAnalyses(noYeastGOAnalysis, noRodentGOAnalysis, c("HI", "HV", "IV", "CH"), pCuttof = 0.1)
goMissingInOldAnalysis = compareAnalyses(noRodentGOAnalysis, noYeastGOAnalysis, c("HI", "HV", "IV", "CH"), pCuttof = 0.1)
geneMissingInNewAnalysis = compareAnalyses(paperAnalysis, newAnalysis)


goMissingInNewAnalysis = compareAnalyses(noYeastGOAnalysis, noRodentGOAnalysis, c("HI"), pCuttof = 0.1)

goMissingInNewAnalysis$missingPadjs

rownames(noRodentGOAnalysis)[goMissingInNewAnalysis$missingPadjs[[1]]]
rownames(noRodentGOAnalysis)[goMissingInOldAnalysis$missingPadjs[[1]]]


goMissingInNewAnalysis = compareAnalyses(paperGOAnalysis, newGOAnalysis, c("HI", "HV", "IV", "CH"), pCuttof = 0.1)
goMissingInOldAnalysis = compareAnalyses(newGOAnalysis, paperGOAnalysis, c("HI", "HV", "IV", "CH"), pCuttof = 0.1)
geneMissingInNewAnalysis = compareAnalyses(paperAnalysis, newAnalysis)






test$mismatchedPadjs[[1]]


rownames(ogGOAnalysis)[test$mismatchedPadjs[[1]]]

sum(newAnalysis$`HV-significant`)
sum(ogGOAnalysis$`HV-significant`)

sum(ogAnalysis[[4]][,3] < 0.05, na.rm = T)
sum(newAnalysis[[4]][,3] < 0.05, na.rm = T)

compare = compareAnalyses(ogAnalysis, newAnalysis)







difGenes = rownames(ogAnalysis[[1]])[compare$mismatchedPadjs[[1]]]


ogRERs = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceRERFile.rds")
commonOgRERs = ogRERs
colnames(commonOgRERs) = ZonomNameConvertVectorCommon(colnames(commonOgRERs), tipColumn = "ZoonomiaTip")
ogPaths = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCategoricalPathsFile.rds")



familyAgnosticRERs = readRDS("Output/CategoricalInsVertivoreTreeFamilyAgnostictLiamInference/CategoricalInsVertivoreTreeFamilyAgnostictLiamInferenceRERFile.rds")
commonFaRERs = familyAgnosticRERs
colnames(commonFaRERs) = ZonomNameConvertVectorCommon(colnames(commonFaRERs), tipColumn = "ZoonomiaTip")
familyAgnosticPaths = readRDS("Output/CategoricalInsVertivoreTreeFamilyAgnostictLiamInference/CategoricalInsVertivoreTreeFamilyAgnostictLiamInferenceCategoricalPathsFile.rds")

palette(c( "darkgreen", "darkblue","black", "red"))
library(gridExtra)




i=1
i = i+1
{
a = plotRers(commonOgRERs, difGenes[i], ogPaths)
b = plotRers(commonFaRERs, difGenes[i], familyAgnosticPaths)
grid.arrange(a,b, ncol=2)

}
ogAnalysis[[1]][compare$mismatchedPadjs[[1]][i],]
familyAgnosticAnalysis[[1]][compare$mismatchedPadjs[[1]][i],]




{
  a = grid.arrange(plotRers(commonOgRERs, difGenes[i], ogPaths),plotRers(commonFaRERs, difGenes[i], familyAgnosticPaths), ncol=2)
  b = grid.arrange(plotRers(ogRERs, difGenes[i], ogPaths),plotRers(familyAgnosticRERs, difGenes[i], familyAgnosticPaths), ncol=2)
  grid.arrange(a,b, nrow=2)
}

significanceResults = geneSignificanceResults

#---------------------------------------------------------------------
# --- Cluster debugging--- 
# --------------------------------------------------------------------
library(xlsx)

test = enrichmentResult[[1]]

testdf = data.frame(c(1,2,3), c("a","b","c"))

write.xlsx(enrichmentResult[[1]], file=enrichmentCsvName, sheetName=enrichmentListName, row.names=T)


args = commandArgs()
marker= "s"                                                         #send the marker value
markerWhole = paste("^", marker, "=", sep='')                    #convert marker to grep format
commandLineValue = grep(markerWhole, args, value = TRUE)         #get a string based on the identifier   
commandLineValue


#---------------------------------------------------------------------
# --- Setup for split-both  // finding out no specific herbivory data --- 
# --------------------------------------------------------------------

treeImageFilename = paste(outputFolderName, filePrefix, "CategoricalTree.pdf", sep="") #make a filename based on the prefix
palette(c(  "darkblue","black", "lightgreen", "darkgreen", "red"))

pdf(treeImageFilename, height = length(phenotypeVector)/18, width = 10)                     #make a pdf to store the plot, sized based on tree size
commonCategoricalTree = char2TreeCategorical(commonPhenotypeVector, commonMainTrees, commonSpeciesFilter, model = modelType, anctrait = ancestralTrait, plot = F)
categoricalTree = char2TreeCategorical(phenotypeVector, mainTrees, speciesFilter, model = modelType, anctrait = ancestralTrait, plot = F) #use the phenotype vector to make a tree

commonCategoricalTreeExtraTip = commonCategoricalTree
categoricalTreeExtraTip = categoricalTree

commonCategoricalTree = drop.tip(commonCategoricalTree, names(nodesToAdd))
categoricalTree = drop.tip(categoricalTree, names(nodesToAdd))
mainTrees$masterTree = drop.tip(mainTrees$masterTree, names(nodesToAdd))
commonMasterAdded = commonMainTrees$masterTree
commonMainTrees$masterTree = drop.tip(commonMainTrees$masterTree, names(nodesToAdd))

plotTreeCategorical(commonCategoricalTree, c("Insectivore", "Omnivore", "SugarHigh", "SugarLow", "Vertivore"), master = commonMainTrees$masterTree)
plotTreeCategorical(categoricalTree, c("Insectivore", "Omnivore", "SugarHigh", "SugarLow", "Vertivore"), master = mainTrees$masterTree)

plotTreeCategorical(commonCategoricalTreeExtraTip, c("Insectivore", "Omnivore", "SugarHigh", "SugarLow", "Vertivore"), master = commonMasterAdded)
plotTreeCategorical(categoricalTreeExtraTip, c("Insectivore", "Omnivore", "SugarHigh", "SugarLow", "Vertivore"), master = masterTreeAdded)
dev.off()  


mergedData = read.csv("Data/mergedData.csv")
fullTree = readRDS("Output/PredatorFullTree/PredatorFullTreeCategoricalTree.rds")

phenData = mergedData[c(mergedData$ZoonomiaTip %in% fullTree$tip.label),]

fullPhenVec = readRDS("Output/PredatorFullTree/PredatorFullTreeCategoricalPhenotypeVector.rds")
dataPhenVec = fullPhenVec[names(fullPhenVec) %in% fullTree$tip.label]
phenData = phenData %>% mutate(diet = dataPhenVec)

herbData = phenData[phenData$diet == "Herbivore",]
nrow(herbData)
length(which(herbData$Diet.PlantO > 99))

length(which(herbData$Diet.PlantO > 99))/nrow(herbData)

length(which(herbData$Diet.PlantO < 10))
otherHerbs = herbData[which(herbData$Diet.PlantO < 10),]

table(phenData$DerekDietClassification90InsVertivoreSorting)

table(phenData$DerekDietClassification90InsVertivoreSorting)
#---------------------------------------------------------------------
# --- Comapring the percentages of species in the pruned over-large families --- 
# --------------------------------------------------------------------

plot(categoricalTree)

demoTree = categoricalTree


demoTree = drop.tip(demoTree, demoTree$tip.label[-seq(from = 2, to = 196, by = 2)])

plotTreeCategorical(demoTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = mainTrees$masterTree)



demoTree$edge.length


highFamilyTips = c("vs_HLellTal1", "vs_HLellLut1", "vs_HLarvAmp1","vs_HLmicAgr2", "vs_HLmyoGla2", "vs_HLondZib1", "voleClade", "6Herb",
  "vs_HLhysCri1", "vs_HLthrSwi1", "vs_HLpetTyp1", "vs_hetGla2", "vs_chiLan1", "vs_HLdinBra1", "vs_HLcteSoc1", "vs_octDeg1", "vs_HLcoePre1", "vs_HLdasPun1", "vs_HLdolPat1", "vs_HLmyoCoy1", "vs_HLhydHyd1", "vs_HLcavTsc1", "gundiGuineaPigClade", "14Herb",
  "vs_HLoryGaz1", "vs_HLbeaHun1", "vs_HLkobLecLec1", "vs_HLmadKir1", "vs_HLneoPyg1", "vs_HLphiMax1", "vs_HLoreOre1", "vs_HLneoMos1", "vs_HLaepMel1", "vs_HLtraImb1",  "vs_bisBis1", "vs_HLoviNivLyd1", "vs_HLproPrz1", "Bovidae", "13Herb",
  "vs_HLhydIne1", "vs_HLmunMun1", "vs_HLodoHem1", "vs_HLantAme1", "vs_HLgirCam1", "Cervidae", "5Herb",
  "vs_HLmacSob1", "vs_HLpteGig1", "FoxLongTounge", "vs_HLeidHel2", "vs_HLcynBra1", "outerPeropodidae", "vs_HLeonSpe1", "vs_HLrouLes1", "Roussetinae", "Peropodidae", "6Herb")


categoricalTree

length(which(highFamilyTips %in% categoricalTree$tip.label))

mergedData = read.csv("Data/mergedData.csv")

mergedData$ZoonomiaTip %in% mainTrees$masterTree$tip.label

fullTree = readRDS("Output/PredatorFullTree/PredatorFullTreeCategoricalTree.rds")
fullPhenVec = readRDS("Output/PredatorFullTree/PredatorFullTreeCategoricalPhenotypeVector.rds")
length(fullPhenVec)

ogPhenVec = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCategoricalPhenotypeVector.rds")
noFamPhenVec  = readRDS("Output/CategoricalInsVertivoreTreeFamilyAgnostictLiamInference/CategoricalInsVertivoreTreeFamilyAgnostictLiamInferenceCategoricalPhenotypeVector.rds")


dataPhenVec = fullPhenVec[names(fullPhenVec) %in% fullTree$tip.label]

phenData = mergedData[c(mergedData$ZoonomiaTip %in% fullTree$tip.label),]


phenData = phenData %>% mutate(diet = dataPhenVec)

familyByDiet = phenData[,c(1,2,5,8,9)]
familyByDiet = familyByDiet %>% mutate(diet = dataPhenVec)

table(familyByDiet$diet)
table(familyByDiet$MSWC_Family[familyByDiet$diet == "Herbivore"])[order(table(familyByDiet$MSWC_Family[familyByDiet$diet == "Herbivore"]))]
length(table(familyByDiet$MSWC_Family[familyByDiet$diet == "Herbivore"])[order(table(familyByDiet$MSWC_Family[familyByDiet$diet == "Herbivore"]))])

table(familyByDiet$MSWC_Family[familyByDiet$diet == "Insectivore"])[order(table(familyByDiet$MSWC_Family[familyByDiet$diet == "Insectivore"]))]
length(table(familyByDiet$MSWC_Family[familyByDiet$diet == "Insectivore"])[order(table(familyByDiet$MSWC_Family[familyByDiet$diet == "Insectivore"]))])


table(familyByDiet$MSWC_Family)[order(table(familyByDiet$MSWC_Family))]

which(familyByDiet$ZoonomiaTip == "vs_HLellTal1")
familyByDiet[62,]

set = which(familyByDiet$ZoonomiaTip %in% c("vs_HLellTal1", "vs_HLellLut1", "vs_HLarvAmp1","vs_HLmicAgr2", "vs_HLmyoGla2", "vs_HLondZib1", "voleClade", "6Herb"))
set = which(familyByDiet$ZoonomiaTip %in% c(  "vs_HLhysCri1", "vs_HLthrSwi1", "vs_HLpetTyp1", "vs_hetGla2", "vs_chiLan1", "vs_HLdinBra1", "vs_HLcteSoc1", "vs_octDeg1", "vs_HLcoePre1", "vs_HLdasPun1", "vs_HLdolPat1", "vs_HLmyoCoy1", "vs_HLhydHyd1", "vs_HLcavTsc1", "gundiGuineaPigClade", "14Herb"))
set = which(familyByDiet$ZoonomiaTip %in% c(          "vs_HLmusSpi1", "vs_HLmusCar1", "vs_HLmasCou1", "vs_HLmusPah1", "vs_HLratNor7", "vs_HLarvNil1", "vs_mm10", "mouseClade", "7Omni,1Herb"))
set = which(familyByDiet$ZoonomiaTip %in% c("vs_HLmusSpi1", "vs_HLmusCar1", "vs_HLmasCou1", "vs_HLmusPah1", "mouseClade"))

familyByDiet[set,]

phenData[which(phenData$ZoonomiaTip %in% c("vs_HLhysCri1", "vs_HLthrSwi1", "vs_HLpetTyp1", "vs_hetGla2", "vs_chiLan1", "vs_HLdinBra1", "vs_HLcteSoc1", "vs_octDeg1", "vs_HLcoePre1", "vs_HLdasPun1", "vs_HLdolPat1", "vs_HLmyoCoy1", "vs_HLhydHyd1", "vs_HLcavTsc1", "gundiGuineaPigClade", "14Herb")),]

prunedFamilyByDiet = familyByDiet
prunedFamilyByDiet = prunedFamilyByDiet[prunedFamilyByDiet$ZoonomiaTip %in% names(ogPhenVec),]


table(prunedFamilyByDiet$diet)
table(prunedFamilyByDiet$MSWC_Family[prunedFamilyByDiet$diet == "Herbivore"])[order(table(prunedFamilyByDiet$MSWC_Family[prunedFamilyByDiet$diet == "Herbivore"]))]
length(table(prunedFamilyByDiet$MSWC_Family[prunedFamilyByDiet$diet == "Herbivore"])[order(table(prunedFamilyByDiet$MSWC_Family[prunedFamilyByDiet$diet == "Herbivore"]))])

table(prunedFamilyByDiet$MSWC_Family[prunedFamilyByDiet$diet == "Insectivore"])[order(table(prunedFamilyByDiet$MSWC_Family[prunedFamilyByDiet$diet == "Insectivore"]))]
length(table(prunedFamilyByDiet$MSWC_Family[prunedFamilyByDiet$diet == "Insectivore"])[order(table(prunedFamilyByDiet$MSWC_Family[prunedFamilyByDiet$diet == "Insectivore"]))])


table(prunedFamilyByDiet$MSWC_Family)[order(table(prunedFamilyByDiet$MSWC_Family))]

noFamPruningFamilyByDiet = familyByDiet
noFamPruningFamilyByDiet = noFamPruningFamilyByDiet[noFamPruningFamilyByDiet$ZoonomiaTip %in% names(noFamPhenVec),]
table(noFamPruningFamilyByDiet$diet)
table(noFamPruningFamilyByDiet$MSWC_Family[noFamPruningFamilyByDiet$diet == "Herbivore"])[order(table(noFamPruningFamilyByDiet$MSWC_Family[noFamPruningFamilyByDiet$diet == "Herbivore"]))]
table(noFamPruningFamilyByDiet$MSWC_Family[noFamPruningFamilyByDiet$diet == "Insectivore"])[order(table(noFamPruningFamilyByDiet$MSWC_Family[noFamPruningFamilyByDiet$diet == "Insectivore"]))]


table(familyByDiet$MSWC_Family[familyByDiet$diet == "Omnivore"])[order(table(familyByDiet$MSWC_Family[familyByDiet$diet == "Omnivore"]))]

table(noFamPruningFamilyByDiet$MSWC_Family[noFamPruningFamilyByDiet$diet == "Omnivore"])[order(table(noFamPruningFamilyByDiet$MSWC_Family[noFamPruningFamilyByDiet$diet == "Omnivore"]))]
length(table(noFamPruningFamilyByDiet$MSWC_Family[noFamPruningFamilyByDiet$diet == "Omnivore"])[order(table(noFamPruningFamilyByDiet$MSWC_Family[noFamPruningFamilyByDiet$diet == "Omnivore"]))])







#---------------------------------------------------------------------
# --- Writing code to compare phenotypes --- 
# --------------------------------------------------------------------


ogAnalysis = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferencePairwiseCorrelationFile.rds")
newAnalysis = readRDS("Output/CategoricalInsVertivoreTreeNoYeastLiamInference/CategoricalInsVertivoreTreeNoYeastLiamInferencePairwiseCorrelationFile.rds")


ogAnalysis = ogAnalysis[-c(7,8)]

ogAnalysis = ogAnalysis[-c(2,3,6)]
newAnalysis = newAnalysis[-c(2,3,6)]

out = compareAnalyses(ogAnalysis, newAnalysis)

compareAnalyses = function(ogAnalysis, newAnalysis){

results <- tibble(
  index = integer(),
  correlation = numeric(),
  rhoRankcorrelation = numeric(),
  numMismatchedNAs = integer(),
  mismatchedNAs = I(list()),
  pCorrelation = numeric(),
  numMismatchedPs = integer(),
  mismatchedPs = I(list()),
  padjCorrelation = numeric(),
  numMismatchedPadjs = integer(),
  mismatchedPadjs = I(list()),
  totalOgSigPadjs = integer(), 
  totalNewSigPadjs = integer()
)

for (i in seq_along(ogAnalysis)) {
  
  correlation <- cor(ogAnalysis[[i]][[1]], newAnalysis[[i]][[1]], use = "complete.obs")
  
  mismiatchedNAs <- which(!is.na(ogAnalysis[[i]][[1]]) %in% is.na(newAnalysis[[i]][[1]]))
  
  ogAnalysis[[i]] = ogAnalysis[[i]] %>% mutate(rhoRank = rank(ogAnalysis[[i]][[1]])) 
  newAnalysis[[i]] = newAnalysis[[i]] %>% mutate(rhoRank = rank(newAnalysis[[i]][[1]]))
  
  rhoRankcorrelation <- cor(ogAnalysis[[i]][[4]], newAnalysis[[i]][[4]], use = "complete.obs")
  
  pCorrelation <- cor(ogAnalysis[[i]][[2]], newAnalysis[[i]][[2]], use = "complete.obs")
  
  mismiatchedPs <- which(!which(ogAnalysis[[i]][[2]] < 0.05) %in% 
                           which(newAnalysis[[i]][[2]] < 0.05))
  
  padjCorrelation <- cor(ogAnalysis[[i]][[3]], newAnalysis[[i]][[3]], use = "complete.obs")
  
  mismiatchedPadjs <- which(!which(ogAnalysis[[i]][[3]] < 0.05) %in% 
                              which(newAnalysis[[i]][[3]] < 0.05))
  
  totalOgSigPadjs = length(which(ogAnalysis[[i]][[3]] < 0.05))
  totalNewSigPadjs = length(which(newAnalysis[[i]][[3]] < 0.05))
  
  results <- rbind(results, tibble(
    index = i,
    correlation = correlation,
    rhoRankcorrelation = rhoRankcorrelation,
    numMismatchedNAs = length(mismiatchedNAs),
    mismatchedNAs = list(mismiatchedNAs),
    pCorrelation = pCorrelation,
    numMismatchedPs = length(mismiatchedPs),
    mismatchedPs = list(mismiatchedPs),
    padjCorrelation = padjCorrelation,
    numMismatchedPadjs = length(mismiatchedPadjs),
    mismatchedPadjs = list(mismiatchedPadjs),
    totalOgSigPadjs = totalOgSigPadjs,
    totalNewSigPadjs = totalNewSigPadjs
  ))
}
row.names(results) = names(ogAnalysis)
return(results)
}



#---------------------------------------------------------------------
# --- comparing new and old carnivory  --- 
# --------------------------------------------------------------------

mergedPairwiseCorrelation = readRDS(paste0(pairwiseCorrelationFileName, ".rds"))
mergedPairwiseCorrelation = pairwiseCategorical
ogCarn = mainPairwiseCategorical[c(7,8,2)]

names(ogCarn)
names(mergedPairwiseCorrelation)

carnCompare = compareAnalyses(ogCarn, mergedPairwiseCorrelation)

all.equal(ogCarn[[1]][[1]], mergedPairwiseCorrelation[[1]][[1]])


ogCarnTree = readRDS("Output/CategoricalInsVertivoreTreeCarnivoreLiamInference/CategoricalInsVertivoreTreeCarnivoreLiamInferenceCategoricalTree.rds")
newCarnTree = readRDS("Output/CategoricalInsVertivoreTreeDuplicateLiamInference/CategoricalInsVertivoreTreeDuplicateLiamInferenceMergedCategoricalTree.rds")

ogLiamTree  = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCategoricalTree.rds")
dupLiamTree = readRDS("Output/CategoricalInsVertivoreTreeDuplicateLiamInference/CategoricalInsVertivoreTreeDuplicateLiamInferenceCategoricalTree.rds")

ogLiamTree$tip.label[which(!ogLiamTree$tip.label %in% dupLiamTree$tip.label)]


all.equal(ogLiamTree, categoricalTree)

phenotypeTree$edge.length
#---------------------------------------------------------------------
# --- Looking into the number of species in various trees--- 
# --------------------------------------------------------------------



testTree = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCategoricalTree.rds")

length(testTree$tip.label)

#---------------------------------------------------------------------
# --- test code during creation of generate alternates --- 
# --------------------------------------------------------------------

ogAnnotCOl = manualAnnots[[annotColumn]] 
newAnnotCol = manualAnnots[[annotColumn]]
new2AnnotCOl = manualAnnots[[annotColumn]]

all.equal(new2AnnotCOl, newAnnotCol)

nrow(relevantSpecies)
nrow(manualAnnots)


nonrelSpecies = manualAnnots[!manualAnnots[[annotColumn]] %in% categoryList,]
nrow(nonrelSpecies)


"vs_NA" %in% categoricalTree$tip.label
"vs_NA" %in% droppedTips
"vs_NA" %in% names(fullDataPhenotype)
"vs_NA" %in% mainTrees$masterTree$tip.label

savedSpeciesSet = randomizedSpeciesSet
saved2 = randomizedSpeciesSet
saved3 = randomizedSpeciesSet

all.equal(saved2, saved3)


length(which(alternateSets[[1]] %in% alternateSets[[6]]))/length(alternateSets[[1]])


names(phenotypeVector)[which(!names(phenotypeVector) %in% categoricalTree$tip.label)]
#---------------------------------------------------------------------
# --- Try to debug errors with partially pruned liam tree --- 
# --------------------------------------------------------------------
c(339,347,351,353,387,388,389,393,509,513,514,515,584,586,587,608,609,611,614,617,620)

c(339,347,351,353,387,388,389,393,509,513,514,515,584,586,587,608,609,611,614,617,620)



prunedSpecies

commonPrunedSpecies = ZonomNameConvertVectorCommon(prunedSpecies, annotationLocation = spreadSheetLocation, tipColumn = nameColumn)

length(commonPrunedSpecies)
length(commonSpeciesFilter)

commonSpeciesFilter %in% commonPrunedSpecies

commonCategoricalTree$tip.label %in% commonPrunedSpecies


plot(commonCategoricalTree)

ZonomNameConvertVectorCommon("vs_HLpanOnc1", annotationLocation = spreadSheetLocation, tipColumn = nameColumn)

commonCategoricalTree

"vs_HLpumYag1"

report= mainTrees$report
speciesGeneNumber = colSums(report)

speciesGeneNumber[which(names(speciesGeneNumber) == "vs_HLpumYag1")]

commonCategoricalTree$edge.length
edgelabels(bg = NULL, adj = c(0.5,0.9), frame = 'none', font =2)


test = plotTreeCategorical(commonCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = commonMainTrees$masterTree)



#---------------------------------------------------------------------
# --- Associate drivingbrahces with eltontraits values  --- 
# --------------------------------------------------------------------
source("Src/Loc/Dev/DisplayCategoricalRERTree.R")
mainTrees = readRDS("data/zoonomiaAllMammalsTrees.rds")
RERObject = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceRERFile.rds")
pathsObject = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCategoricalPathsFile.rds")
palette(c( "darkgreen", "darkblue","black", "red"))

indexVal = "IDO2"


rerTree = displayCategoricalRERTree(mainTrees, RERObject, indexVal, phenv = pathsObject)

terminalBranches = which(rerTree$edge[,2] < length(rerTree$tip.label))
names(terminalBranches) = rerTree$tip.label[rerTree$edge[,2][which(rerTree$edge[,2] < length(rerTree$tip.label))]]
terminalBranches = rerTree$edge.length[terminalBranches]

mergedData = read.csv("Data/MergedData.csv")
mainData = mergedData
maindata =  mainData[which(mainData$ZoonomiaTip %in% names(terminalBranches)),]
terminalBranches = terminalBranches[match(maindata$ZoonomiaTip, names(terminalBranches))]

maindata = cbind(maindata, terminalBranches)

ggplot(maindata, aes(x = terminalBranches, y = Diet.Scav)) +
  geom_point()

#---------------------------------------------------------------------
# --- useDrivingBranches script  --- 
# --------------------------------------------------------------------
source("Src/Loc/Dev/DisplayCategoricalRERTree.R")
mainTrees = readRDS("data/zoonomiaAllMammalsTrees.rds")
RERObject = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceRERFile.rds")
pathsObject = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCategoricalPathsFile.rds")
palette(c( "darkgreen", "darkblue","black", "red"))

indexVal = "SDS"



displayCategoricalRERTree(mainTrees, RERObject, indexVal, phenv = pathsObject, tipCol = "ZoonomiaTip")
pdf(height = 30, width = 15)
displayCategoricalRERTree(mainTrees, RERObject, indexVal, phenv = pathsObject, tipCol = "ZoonomiaTip")
dev.off()

#---------------------------------------------------------------------
# --- Making script ot look into driving branches in result --- 
# --------------------------------------------------------------------

mainTrees = readRDS("data/zoonomiaAllMammalsTrees.rds")
RERObject = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceRERFile.rds")
pathsObject = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCategoricalPathsFile.rds")
#source("Src/Reu/RERConvergeFunctions.R")

source("Src/Reu/ZonomNameConvertVectorCommon.R")

tipCol = "ZoonomiaTip"


treesObj = mainTrees
rermat = RERObject
phenv = pathsObject
index = "SDS"

displayCategoricalRERTree = function(treesObj, rermat, index, phenv = NULL, subsetTree = T, equalLengths = T, minWidth = 1, maxWidth = 6, tipCol = "tipColumn"){
  treesObj$trees[[index]]$tip.label = ZonomNameConvertVectorCommon(treesObj$trees[[index]]$tip.label, tipColumn = tipCol)
  treesObj$masterTree$tip.label = ZonomNameConvertVectorCommon(treesObj$masterTree$tip.label, tipColumn = tipCol)
  colnames(rermat) = ZonomNameConvertVectorCommon(colnames(rermat), tipColumn = tipCol)
  returnRersAsTreeNew(treesObj, rermat, index, phenv, 0.7, 0.7, 'NA', T, subsetTree, equalLengths, minWidth, maxWidth)
  
  
  
  
}



returnRersAsTreeNew <- function(treesObj, rermat, index, phenv = NULL, rer.cex = 0.7,
                             tip.cex = 0.7, nalab = 'NA', plot = T, subsetTree = F, equalLengths=T, minWidth=1, maxWidth=6){
  trgene <- treesObj$trees[[index]]
  if (subsetTree == TRUE) {
    #Use phenv to identify non-NA paths and then subset gene tree by species
    subsp <- unique(colnames(rermat)[which(!is.na(phenv))])
    trgene <- pruneTree(trgene, subsp) #drops only species not in subsp
  }
  if(equalLengths){trgene$edge.length <- rep(2,nrow(trgene$edge))}
  ee=edgeIndexRelativeMaster(trgene, treesObj$masterTree)
  ii= treesObj$matIndex[ee[, c(2,1)]]
  rertree=rermat[index,ii]
  rertree[is.nan(rertree)]=NA #replace NaNs from C functions
  relativeRER = abs(scale(rertree))
  if(any(relativeRER > 3)){
    message("WARNING: Potential misannotated gene!")
    for(i in which(relativeRER > 3)){
      if(!is.na(rownames(relativeRER)[i])){
        message(paste("Potential Misannotated gene at branch:", i ))
        message(paste("Species of branch=:", rownames(relativeRER)[i]))
        if(plot){
          message("Branch highlighted in PINK.")
        }
      }
    }
  }
  
  if (plot) {
    clampedRelativeRER = relativeRER
    clampedRelativeRER[order(relativeRER, decreasing = TRUE)[1:3]] <- sort(relativeRER, decreasing = TRUE)[4] #Set the top three high RERs to the fourth highest RER; this matches the calculation ignoring top values and prevents mis-IDed genes from throwing off the scale. 
    rerWidth = rescale(clampedRelativeRER, c(minWidth,maxWidth))
    rerWidth[is.na(rerWidth)] = 1
    
    par(mar = c(1,1,1,0))
    edgcols <- rep('black', nrow(trgene$edge))
    edgwds <- rep(1, nrow(trgene$edge))
    if(!is.null(phenv)){
      edgcols <- rep('black', nrow(trgene$edge))
      edgwds <- rerWidth
      if(length(unique(phenv) < length(palette()))){ # add a catch for continuous phenotypes and not run it in that case
        for(j in unique(phenv)[!is.na(unique(phenv))]){
          edgcols[phenv[ii]==j] <- palette()[j]
        }      
      }
    }
    if(any(relativeRER > 3)){
      for(i in which(relativeRER > 3)){
        if(!is.na(rownames(relativeRER)[i])){
          edgcols[i] = 'hotpink'
        }
      }
    }
    plot.phylo(trgene, font = 2, edge.color = edgcols, edge.width = edgwds, cex = tip.cex)
    rerlab <- round(rertree,3)
    rerlab[is.na(rerlab)] <- nalab
    if(!is.null(phenv)){ #reset the color so that the RER labels still have the correct palette color even with the pink warning
      for(j in unique(phenv)[!is.na(unique(phenv))]){
        edgcols[phenv[ii]==j] <- palette()[j]
      }
    }
    edgelabels(rerlab, bg = NULL, adj = c(0.5,0.9), col = edgcols, frame = 'none',cex = rer.cex, font =2)
  }

  
  trgene$edge.length <- rertree
  return(trgene)
}

{
  edgeIndexRelativeMaster= function(tree, masterTree){
    map=matchAllNodes(tree,masterTree)
    newedge=tree$edge
    newedge[,1]=map[newedge[,1],2]
    newedge[,2]=map[newedge[,2],2]
    newedge
  }
  matchAllNodes = function(tree1, tree2){
    map=matchNodesInject(tree1,tree2)
    map=map[order(map[,1]),]
    map
  }
  matchNodesInject = function (tr1, tr2){
    if(length(tmpsp<-setdiff(tr1$tip.label, tr2$tip.label))>0){
      #stop(paste(paste(tmpsp, ","), "in tree1 do not exist in tree2"))
      stop(c("The following species in tree1 do not exist in tree2: ",paste(tmpsp, ", ")))
    }
    commontiplabels <- intersect(tr1$tip,tr2$tip)
    if(RF.dist(pruneTree(tr1,commontiplabels),pruneTree(tr2,commontiplabels))>0){
      stop("Discordant tree topology detected - gene/trait tree and treesObj$masterTree have irreconcilable topologies")
    }
    #if(RF.dist(tr1,tr2)>0){
    #  stop("Discordant tree topology detected - trait tree and treesObj$masterTree have irreconcilable topologies")
    #}
    
    toRm=setdiff(tr2$tip.label, tr1$tip.label)
    desc.tr1 <- lapply(1:tr1$Nnode + length(tr1$tip), function(x) extract.clade(tr1,
                                                                                x)$tip.label)
    names(desc.tr1) <- 1:tr1$Nnode + length(tr1$tip)
    desc.tr2 <- lapply(1:tr2$Nnode + length(tr2$tip), function(x) extract.clade(tr2,
                                                                                x)$tip.label)
    names(desc.tr2) <- 1:tr2$Nnode + length(tr2$tip)
    Nodes <- matrix(NA, length(desc.tr1), 2, dimnames = list(NULL,
                                                             c("tr1", "tr2")))
    for (i in 1:length(desc.tr1)) {
      Nodes[i, 1] <- as.numeric(names(desc.tr1)[i])
      for (j in 1:length(desc.tr2)) if (all(desc.tr1[[i]] %in%
                                            desc.tr2[[j]]))
        Nodes[i, 2] <- as.numeric(names(desc.tr2)[j])
    }
    
    iim=match(tr1$tip.label, tr2$tip.label)
    Nodes=rbind(cbind(1:length(tr1$tip.label),iim),Nodes)
    if(any(table(Nodes[,2])>1)){
      stop("Incorrect pseudorooting detected - use fixPseudoroot() function to correct trait tree topology")
    }
    
    Nodes
  }
  
}

pdf(height = 30, width = 15)
returnRersAsTreeNew(treesObj = mainTrees, rermat = RERObject, index = "SDS", phenv = pathsObject, subsetTree = T)
dev.off()

displayCategoricalRERTree(treesObj = mainTrees, rermat = RERObject, index = "SDS", phenv = pathsObject, subsetTree = T, tipCol = "ZoonomiaTip")



#---------------------------------------------------------------------
# --- work on making new p value caluclation script --- 
# --------------------------------------------------------------------

realCors = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCombinedCategoricalCorrelationFile.rds")
intermediateList = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCategoricalPermulationsIntermediates7.rds")

CategoricalCollectIntermediateResults = function(realCors, intermediateList, initial = F, start=1, end=NULL, report=F){
  {totalStart = Sys.time()}
  corsMatEffSize = intermediateList[[1]]
  Peffsize = intermediateList[[2]]
  corsMatPvals = intermediateList[[3]]
  Ppvals = intermediateList[[4]]
  message("Obtaining permulations p-values")
  N = nrow(realCors[[1]]) #
  if(initial){ #Only do this if start = 1, because otherwise it's already made and you'll overwrite the old script's results 
    realCors[[1]]$permP = rep(NA, N) #Make a column for permP values in all of the dataframes 
    realCors[[1]]$numMoreExtremePerms = rep(0, N) #Make a column for permP values in all of the dataframes 
    realCors[[1]]$numTotalPerms = rep(0, N) #Make a column for permP values in all of the dataframes 
    for (j in 1:length(realCors[[2]])) {
      realCors[[2]][[j]]$permP = rep(NA, N) #Make a column for permP values in all of the dataframes 
      realCors[[2]][[j]]$numMoreExtremePerms = rep(0, N) #Make a column for permP values in all of the dataframes 
      realCors[[2]][[j]]$numTotalPerms = rep(0, N) #Make a column for permP values in all of the dataframes 
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
  return(list(res = realCors, pvals = list(corsMatPvals, Ppvals), effsize = list(corsMatEffSize, Peffsize)))
  
}

test = CategoricalCollectIntermediateResults(realCors, intermediateList, report = T)

CollectedIntermediates = test

CategoricalCalculatePValueFromCollectedIntermediates = function(CollectedIntermediates, start=1, end=NULL, report=F){
  realCors = CollectedIntermediates[[1]]
  N = nrow(realCors[[1]]) #
  totalStart = Sys.time()

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
  {totalEnd = Sys.time(); totalDuration = totalEnd - totalStart;message(paste("Completed p-Values; Duration", totalDuration, attr(totalDuration, "units")))}
  
  CollectedIntermediates[[1]] = realCors
  
  return(CollectedIntermediates)
}


test2 = CategoricalCalculatePValueFromCollectedIntermediates(CollectedIntermediates)
test2$res[[1]]

test[[1]]

if(calulateValue){
  source("Src/Reu/CategoricalPermulationsParallelFunctions.R")
  
  #Correlations
  correlationFileName = paste(outputFolderName, filePrefix, "CombinedCategoricalCorrelationFile.rds", sep= "") #Make a correlation filename based on the prefix
  correlationsObject = readRDS(correlationFileName) 
  
  if(onlyCalulateValue){
    if(metacombineValue == F){
      combinedDataFileName = paste(outputFolderName, filePrefix, "Combined", permulationPrefix,"PermulationsIntermediates", runInstanceValue, ".rds", sep="")
    }else{
      combinedDataFileName = paste(outputFolderName, filePrefix, "MetaCombined", permulationPrefix, "PermulationsIntermediates", runInstanceValue, ".rds", sep="")
    }
    combinedPermulationsData = readRDS(combinedDataFileName)
  }
  permulationsPValues = CategoricalCalculatePermulationPValues(correlationsObject, combinedPermulationsData)
  permulationsPValuesOutput = permulationsPValues$res
  
  #Give pairwise outputs descriptive names
  phenotypeVectorFilename = paste(outputFolderName, filePrefix, "CategoricalPhenotypeVector.rds",sep="") #select the phenotype vector based on prefix
  phenotypeVector = readRDS(phenotypeVectorFilename)                            #load in the phenotype vector 
  categories = map_to_state_space(phenotypeVector)                              #and use it to connect branch lengths to phenotype name
  categoryNames = categories$name2index                                         #store the length-phenotype connection
  
  pairwiseTableNames = names(permulationsPValuesOutput[[2]])                               #Prepare to replace the number-number titles with phenotype-phenotype titles
  for(i in 1:length(categoryNames)){                                            #for each phenotype
    pairwiseTableNames= gsub(i, names(categoryNames)[i], pairwiseTableNames)    #replace the number with the phenotype name  
  }
  names(permulationsPValuesOutput[[2]]) = pairwiseTableNames                               #update the dataframe titles
  
  permulationsPValuesFilename = paste(outputFolderName, filePrefix, "PermulationsPValueCorrelations.rds", sep= "")
  saveRDS(permulationsPValuesOutput, permulationsPValuesFilename)
  
  outputSubdirectoryNoslash = paste(outputFolderName, "Overall", sep = "")
  if(!dir.exists(outputSubdirectoryNoslash)){                       #create that directory if it does not exist
    dir.create(outputSubdirectoryNoslash)
  }
  outputSubdirectory = paste(outputSubdirectoryNoslash, "/", sep="")
  
  permulationsPValuesOverallFilename = paste(outputSubdirectory, filePrefix, "OverallPermulationsCorrelationFile.rds", sep= "")
  saveRDS(permulationsPValuesOutput[[1]], permulationsPValuesOverallFilename)
  
  for(i in 1:length(pairwiseTableNames)){
    pairwiseTableNames= gsub(" ", "", pairwiseTableNames)
    
    outputSubdirectoryNoslash = paste(outputFolderName, pairwiseTableNames[i], sep = "")
    if(!dir.exists(outputSubdirectoryNoslash)){                       #create that directory if it does not exist
      dir.create(outputSubdirectoryNoslash)
    }
    outputSubdirectory = paste(outputSubdirectoryNoslash, "/", sep="")
    
    permulationsPValuesPairFilename = paste(outputSubdirectory, filePrefix, pairwiseTableNames[i], "PermulationsCorrelationFile",".rds", sep= "")
    saveRDS(permulationsPValuesOutput[[2]][[i]], permulationsPValuesPairFilename)
  }
}


#---------------------------------------------------------------------
# --- Encorperate permulation P values --- 
# --------------------------------------------------------------------

topPerms = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferencePermulationsPValueCorrelations.rds")

 
topPerms[[1]]$index = seq_len(nrow(topPerms[[1]]))
topPerms[[2]][[1]]$index = seq_len(nrow(topPerms[[2]][[1]]))
topPerms[[2]][[2]]$index = seq_len(nrow(topPerms[[2]][[2]]))
topPerms[[2]][[3]]$index = seq_len(nrow(topPerms[[2]][[3]]))
topPerms[[2]][[4]]$index = seq_len(nrow(topPerms[[2]][[4]]))
topPerms[[2]][[5]]$index = seq_len(nrow(topPerms[[2]][[5]]))
topPerms[[2]][[6]]$index = seq_len(nrow(topPerms[[2]][[6]]))


overallPerms = topPerms[[1]]
pairwaisePerms = topPerms[[2]]


overallPerms = overallPerms[order(overallPerms$p.adj),]
overallPerms$baseRank = seq_len(nrow(overallPerms))
overallPerms = overallPerms[order(overallPerms$permP),]
overallPerms$permRank = seq_len(nrow(overallPerms))

for(i in 1:6){
  pairwaisePerms[[i]] = pairwaisePerms[[i]][order(pairwaisePerms[[i]]$p.adj),]
  pairwaisePerms[[i]]$baseRank = seq_len(nrow(pairwaisePerms[[i]]))
  pairwaisePerms[[i]] = pairwaisePerms[[i]][order(pairwaisePerms[[i]]$permP),]
  pairwaisePerms[[i]]$permRank = seq_len(nrow(pairwaisePerms[[i]]))
  plot(pairwaisePerms[[i]]$baseRank, pairwaisePerms[[i]]$permRank)
}

plot(pairwaisePerms[[i]]$p.adj, pairwaisePerms[[i]]$permP)



# saveRDS(pairwaisePerms, "Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferencePermulationPairwiseCorrelationFile.rds")


permGO = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferencecombinedGOResults-KeggReactomePermulations.rds")
nonpermGO = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferencecombinedGOResults-KeggReactome.rds")

permGO$`CH-stat` - nonpermGO$`CH-stat`


plot(permGO$`HV-p.adj`, nonpermGO$`HV-p.adj`)

test1 = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/Insectivore-Vertivore/CategoricalInsVertivoreTreeLiamInferenceInsectivore-VertivoreCorrelationFile.rds")
test2 = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/Insectivore-Vertivore/CategoricalInsVertivoreTreeLiamInferenceInsectivore-VertivorePermulationsCorrelationFile.rds")

# ------------------------------------------------------------------
# ---  Make New SUpplementary File 3----- 
# ------------------------------------------------------------------

MGI = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsvertivoreTreeLiamInferencecombinedGOResults-MGI_Mammalian_Phenotype_Level_4.rds")
GO = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsvertivoreTreeLiamInferencecombinedGOResults-GO_Biological_Process_2023.rds")
DisGen = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsvertivoreTreeLiamInferencecombinedGOResults-DisGeNet.rds")
tissue = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsvertivoreTreeLiamInferencecombinedGOResults-tissue_specific.rds")
Enrichment = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsvertivoreTreeLiamInferencecombinedGOResults-EnrichmentHsSymbolsFile2.rds")


MGI$GeneSet = rep("MGI Mammalian Phenotype", nrow(MGI))
GO$GeneSet = rep("GO Biological Process", nrow(GO))
DisGen$GeneSet = rep("DisGeNET", nrow(DisGen))
tissue$GeneSet = rep("Tissue Specific", nrow(tissue))
Enrichment$GeneSet = rep("EnrichmentHsSymbols", nrow(Enrichment))


combinedGenesets = rbind(MGI, GO, DisGen, tissue, Enrichment)

combinedGenesets = combinedGenesets[,c(1,296, 2:295)]
write.csv(combinedGenesets, "Output/CategoricalInsVertivoreTreeLiamInference/CombinedGenesetOutput.csv")

# ------------------------------------------------------------------
# ---  Make newick files for trees ----- 
# ------------------------------------------------------------------

phenotypeTree3Diet = readRDS("Output/CategoricalInsVertivoreTreeCarnivoreLiamInference/CategoricalInsVertivoreTreeCarnivoreLiamInferenceCategoricalCommonTree.rds")
phenotypeTree4Diet = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCategoricalCommonTree.rds")
mainTrees = readRDS("Data/zoonomiaAllMammalsTrees.rds")
masterTree = mainTrees$masterTree


write.tree(phenotypeTree3Diet,"Results/PhenTree3.txt")
write.tree(phenotypeTree4Diet,"Results/PhenTree4.txt")
write.tree(masterTree,"Results/masterTree.txt")



# ------------------------------------------------------------------
# ---  Make Liam Radial Tree ----- 
# ------------------------------------------------------------------

liamTree = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCategoricalTree.rds")
nonliamTree = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalTree.rds")


match(liamTree$tip.label, nonliamTree$tip.label)

liamNonliamTipConversionIndex = match(nonliamTree$tip.label, liamTree$tip.label)


{
collapsedClades = data.frame()
collapsedClades[1,] = NA

collapsedClades$Platypus = MRCA(commonCategoricalTree, 195)
collapsedClades$Opossums = MRCA(commonCategoricalTree, c(194,192,193))
collapsedClades$Koala = MRCA(commonCategoricalTree, c(190,191))
collapsedClades$Kangaroos = MRCA(commonCategoricalTree, c(186,187,189,188))
collapsedClades$Anteaters = MRCA(commonCategoricalTree, c(183,182))
collapsedClades$Sloths = MRCA(commonCategoricalTree, c(181,180))
collapsedClades$Elephant= MRCA(commonCategoricalTree, c(178,177,176))
collapsedClades$Aardvark = MRCA(commonCategoricalTree, c(175,174,171,172,173))
collapsedClades$Strepsirrhini = MRCA(commonCategoricalTree, c(28,29,27,26,24,25,23,21,22,20,19))
collapsedClades$Atelidae = MRCA(commonCategoricalTree, c(6,5,4,3,2,1))
collapsedClades$Chimpanze = MRCA(commonCategoricalTree, c(18,17,16,15,14,12,13,10,9,8,7,11))
collapsedClades$Hares = MRCA(commonCategoricalTree, c(71,70))
collapsedClades$Squirrels = MRCA(commonCategoricalTree, c(69,67,68,61,66,65,63,62,64))
collapsedClades$Capybara = MRCA(commonCategoricalTree, c(33,32,31,30))
collapsedClades$Beaver = MRCA(commonCategoricalTree, c(36,35,34))
collapsedClades$Jerboa = MRCA(commonCategoricalTree, c(59,58,57))
collapsedClades$Deomyinae = MRCA(commonCategoricalTree, c(41,39,40,38,37,44,43,42))
collapsedClades$Vole = MRCA(commonCategoricalTree, c(48,47,46,45))
collapsedClades$Neotominae = MRCA(commonCategoricalTree, c(54,49,50,52,51,53))
collapsedClades$`African Hedgehogs` = MRCA(commonCategoricalTree, c(168,169))
collapsedClades$`Talpa europaea` = MRCA(commonCategoricalTree, c(167,165,164,166))
collapsedClades$`Flying Fox` = MRCA(commonCategoricalTree, c(163,161,162))
collapsedClades$Rhinolophidae = MRCA(commonCategoricalTree, c(155,156,158,157,159,160))
collapsedClades$`Big Brown Bat` = MRCA(commonCategoricalTree, c(143,142,139,140,141))
collapsedClades$Phyllostomidae = MRCA(commonCategoricalTree, c(144,145,147,146))
collapsedClades$Noctilio = MRCA(commonCategoricalTree, c(154))
collapsedClades$Horse = MRCA(commonCategoricalTree, c(100,99,98))
collapsedClades$Pig = MRCA(commonCategoricalTree, c(95,96))
collapsedClades$`bos bison` = MRCA(commonCategoricalTree, c(87,89,88))
collapsedClades$`Humpback Whale` = MRCA(commonCategoricalTree, c(75,74,72,73))
collapsedClades$`Dolphins` = MRCA(commonCategoricalTree, c(82,81,80,77,76))
collapsedClades$`Pangolin` = MRCA(commonCategoricalTree, c(138,137))
collapsedClades$`Lion` = MRCA(commonCategoricalTree, c(131,130,129))
collapsedClades$`Meerkat` = MRCA(commonCategoricalTree, c(135,134))
collapsedClades$`Dog` = MRCA(commonCategoricalTree, c(101,102))
collapsedClades$`Brown Bear` = MRCA(commonCategoricalTree, c(125,128,127))
collapsedClades$`Odobenus rosmarus` = MRCA(commonCategoricalTree, c(120,119,118,117))
collapsedClades$`Phocidae` = MRCA(commonCategoricalTree, c(124,123,121,122))
collapsedClades$`Procyon lotor` = MRCA(commonCategoricalTree, c(114,113,112))
collapsedClades$`Lontra provocax` = MRCA(commonCategoricalTree, c(108,107,106,105,104))
collapsedClades$`Tasmanian Devil` = MRCA(commonCategoricalTree, c(184,185))
}




{
  collapsedClades = data.frame()
  collapsedClades[1,] = NA
  
  collapsedClades$Platypus = MRCA(commonCategoricalTree, c(1))
  collapsedClades$Opossums = MRCA(commonCategoricalTree, c(3,4,5))
  collapsedClades$Koala = MRCA(commonCategoricalTree, c(8,9))
  collapsedClades$Kangaroos = MRCA(commonCategoricalTree, c(10,11,12,13))
  collapsedClades$Anteaters = MRCA(commonCategoricalTree, c(23,24))
  collapsedClades$Sloths = MRCA(commonCategoricalTree, c(25,26))
  collapsedClades$Elephant= MRCA(commonCategoricalTree, c(19,20,21))
  collapsedClades$Aardvark = MRCA(commonCategoricalTree, c(14,15,16,17,18))
  collapsedClades$Strepsirrhini = MRCA(commonCategoricalTree, c(27,28,29,30,31,32,33,34,35,36,37))
  collapsedClades$Atelidae = MRCA(commonCategoricalTree, c(38,39,40,41,42,43))
  collapsedClades$Chimpanze = MRCA(commonCategoricalTree, c(44:55))
  collapsedClades$Hares = MRCA(commonCategoricalTree, c(56,57))
  collapsedClades$Squirrels = MRCA(commonCategoricalTree, c(58:66))
  collapsedClades$Capybara = MRCA(commonCategoricalTree, c(67:70))
  collapsedClades$Beaver = MRCA(commonCategoricalTree, c(72:74))
  collapsedClades$Jerboa = MRCA(commonCategoricalTree, c(75:77))
  collapsedClades$Deomyinae = MRCA(commonCategoricalTree, c(90:97))
  collapsedClades$Vole = MRCA(commonCategoricalTree, c(86:89))
  collapsedClades$Neotominae = MRCA(commonCategoricalTree, c(80:85))
  collapsedClades$`African Hedgehogs` = MRCA(commonCategoricalTree, c(99,100))
  collapsedClades$`Talpa europaea` = MRCA(commonCategoricalTree, c(101:104))
  collapsedClades$`Flying Fox` = MRCA(commonCategoricalTree, c(105:107))
  collapsedClades$Rhinolophidae = MRCA(commonCategoricalTree, c(108:113))
  collapsedClades$`Big Brown Bat` = MRCA(commonCategoricalTree, c(125:129))
  collapsedClades$Phyllostomidae = MRCA(commonCategoricalTree, c(121:124))
  collapsedClades$Noctilio = MRCA(commonCategoricalTree, c(114))
  collapsedClades$Horse = MRCA(commonCategoricalTree, c(168:170))
  collapsedClades$Pig = MRCA(commonCategoricalTree, c(172:173))
  collapsedClades$`bos bison` = MRCA(commonCategoricalTree, c(194:196))
  collapsedClades$`Humpback Whale` = MRCA(commonCategoricalTree, c(175:178))
  collapsedClades$`Dolphins` = MRCA(commonCategoricalTree, c(182:186))
  collapsedClades$`Pangolin` = MRCA(commonCategoricalTree, c(130:131))
  collapsedClades$`Lion` = MRCA(commonCategoricalTree, c(137:139))
  collapsedClades$`Meerkat` = MRCA(commonCategoricalTree, c(135:136))
  collapsedClades$`Dog` = MRCA(commonCategoricalTree, c(140:141))
  collapsedClades$`Brown Bear` = MRCA(commonCategoricalTree, c(142:144))
  collapsedClades$`Odobenus rosmarus` = MRCA(commonCategoricalTree, c(146:149))
  collapsedClades$`Phocidae` = MRCA(commonCategoricalTree, c(150:153))
  collapsedClades$`Procyon lotor` = MRCA(commonCategoricalTree, c(156:158))
  collapsedClades$`Lontra provocax` = MRCA(commonCategoricalTree, c(162:166))
  collapsedClades$`Tasmanian Devil` = MRCA(commonCategoricalTree, c(6:7))
  
  #collapsedClades$Cats = MRCA(commonCategoricalTree, c("Jaguar", "Lion", "Cheetah"))
}


#collapsedClades$Cats = MRCA(commonCategoricalTree, c("Jaguar", "Lion", "Cheetah"))

# ------------------------------------------------------------------
# ---  check seize ----- 
# ------------------------------------------------------------------

phenoTree = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCategoricalTree.rds")

table(phenoTree$edge.length)


data = read.csv("Data/mergedData.csv")
table(data$DerekDietClassification90InsVertivoreSorting)
grep("Piscivore")

# ------------------------------------------------------------------
# ---  Examine cortisol results ----- 
# ------------------------------------------------------------------

which(rownames(GoMetaCombined) == "REACTOME_METABOLISM_OF_STEROID_HORMONES")

which(gmt$geneset.names == "REACTOME_METABOLISM_OF_STEROID_HORMONES")

steroidGenes = gmt$genesets[348][[1]]

steriodCHSignifiance = geneMetaCombined$`CH-significant-Liam`[match(steroidGenes, rownames(geneMetaCombined))]
steriodHISignifiance = geneMetaCombined$`HI-significant-Liam`[match(steroidGenes, rownames(geneMetaCombined))]
steriodHVSignifiance = geneMetaCombined$`HV-significant-Liam`[match(steroidGenes, rownames(geneMetaCombined))]



steriodSignificance = data.frame(steroidGenes, steriodCHSignifiance, steriodHISignifiance, steriodHVSignifiance)

genesAllwaysUnsignificant = which(apply(steriodSignificance[,2:4], 1, function(x) all(x == FALSE)))
genesNA = which(apply(steriodSignificance[,2:4], 1, function(x) all(is.na(x))))


steriodSignificant = steriodSignificance[-c(genesAllwaysUnsignificant, genesNA),]


CHOnlyGenes = steriodSignificant[c(which(apply(steriodSignificant[,3:4], 1, function(x) all(x == FALSE))),which(apply(steriodSignificant[,3:4], 1, function(x) all(x == TRUE)))),]
steriodDifferences = steriodSignificant[-c(which(apply(steriodSignificant[,3:4], 1, function(x) all(x == FALSE))),which(apply(steriodSignificant[,3:4], 1, function(x) all(x == TRUE)))),]




# ------------------------------------------------------------------
# ---  Make new venn diagram ----- 
# ------------------------------------------------------------------

vennInputDataframe = vennGoSignificanceResults

makeVennPlot = function(vennInputDataframe, mainTitle, plot = T){
  # Build logical vectors for each set
  set1 <- vennInputDataframe[1] == TRUE
  set2 <- vennInputDataframe[2] == TRUE
  set3 <- vennInputDataframe[3] == TRUE
  
  # Create the Venn counts for each region
  vennCounts = c(
    "Column One" = sum(set1 & !set2 & !set3, na.rm = T),
    "Column Two" = sum(!set1 & set2 & !set3, na.rm = T),
    "Column Three" = sum(!set1 & !set2 & set3, na.rm = T),
    "Column One&Column Two" = sum(set1 & set2 & !set3, na.rm = T),
    "Column One&Column Three" = sum(set1 & !set2 & set3, na.rm = T),
    "Column Two&Column Three" = sum(!set1 & set2 & set3, na.rm = T),
    "Column One&Column Two&Column Three" = sum(set1 & set2 & set3, na.rm = T)
  )
  
  comparisonPrefixes = gsub("-significant", "", names(vennInputDataframe))
  comparisonPrefixes = gsub(commonBackground, "", comparisonPrefixes)
  comparisonNames = sapply(comparisonPrefixes, replacePrefixWithName)
  names(vennCounts)=c(comparisonNames[1], comparisonNames[2], comparisonNames[3], paste(comparisonNames[1], comparisonNames[2], sep="&"),paste(comparisonNames[1], comparisonNames[3], sep="&"),paste(comparisonNames[2], comparisonNames[3], sep="&"), paste(comparisonNames[1], comparisonNames[2], comparisonNames[3], sep="&"))
  
  # - make venn diagram labels, with the combination section on two lines and the solo sections on one 
  totalVennValues = sum(vennCounts)
  vennLabels = paste0(
    vennCounts, "\n (", round(vennCounts / totalVennValues * 100, 1), "%)"
  )
  for(i in 1:3){
    vennLabels[i] = paste0(
      vennCounts[i], " (", round(vennCounts[i] / totalVennValues * 100, 1), "%)"
    )
  }
  
  
  vennCountsCustom = vennCounts
  
  vennCountsCustom[5] = 80
  vennCountsCustom[1] = 40 
  vennCountsCustom[2] = 15
  
  vennLabelsCustom = vennLabels
  vennLabelsCustom[4] = 3
  
  # Create Euler diagram
  fit = euler(vennCountsCustom)
  outPlot = plot(fit,
                 fills = list(fill = vennColorset, alpha = 0.5),
                 labels = list(font = 4),
                 quantities = list(labels = vennLabelsCustom, font = 3),
                 main = mainTitle, theme)
  if(plot){print(outPlot)}
  return(outPlot)
}

GoVennCustom = outPlot

png(width = 1120, height = 560, file = "Output/CategoricalInsvertivoreTreeLiamInference/VennDiagramFigure.png")
combinedPlot = grid.arrange(geneVenn, GoVennCustom, nrow = 1)
dev.off()

# ------------------------------------------------------------------
# --- Add cytoscape index column ----- 
# ------------------------------------------------------------------

cytoscapeNodes1 = read.csv("Output/CategoricalInsVertivoreTreeLiamInference/Herbivore-Insectivore/Cytoscape/CarnivoreConserveddefaultnode.csv")
cytoscapeNodes2 = read.csv("Output/CategoricalInsVertivoreTreeLiamInference/Herbivore-Insectivore/Cytoscape/Bidirectionalnodes.csv")
cytoscapeNodes3 = read.csv("Output/CategoricalInsVertivoreTreeLiamInference/Herbivore-Insectivore/Cytoscape/Herbivorenodes.csv")
cytoscapeNodes4 = read.csv("Output/CategoricalInsVertivoreTreeLiamInference/Herbivore-Insectivore/Cytoscape/BidectionalNewNodes.csv")
cytoscapeNodes5 = read.csv("Output/CategoricalInsVertivoreTreeLiamInference/Insectivore-Vertivore/Cytoscape/IVNodes.csv")
cytoscapeNodes6 = read.csv("Output/CategoricalInsVertivoreTreeLiamInference/Herbivore-Insectivore/Cytoscape/NewPredatorNodes.csv")
cytoscapeNodes7 = read.csv("Output/CategoricalInsVertivoreTreeLiamInference/Herbivore-Insectivore/Cytoscape/NewBidirectionalNodes.csv")


cytoscapeNodes = cytoscapeNodes7

GOOutput = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsvertivoreTreeLiamInferencecombinedGOResults-KeggReactome.rds")

rownames(GOOutput)

cytoscapeNodes$shared.name

test = match(cytoscapeNodes$shared.name, rownames(GOOutput))

which(rownames(GOOutput) == "KEGG_BETA_ALANINE_METABOLISM")
which(rownames(GOOutput) == "REACTOME_DNA_REPAIR")



rownames(GOOutput)[436]

nodeIdexes = data.frame(cytoscapeNodes$shared.name, (match(cytoscapeNodes$shared.name, rownames(GOOutput)))+1)

write.csv(nodeIdexes, "Output/CategoricalInsVertivoreTreeLiamInference/Herbivore-Insectivore/Cytoscape/fixedBidirectionalIndex.csv")


length(which(GOOutput$`IV-significant`))

length(which(geneMetaCombined$`IV-significant-Liam`))

repairVsH = c("PID_FANCONI_PATHWAY","REACTOME_DISEASES_OF_DNA_REPAIR","REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR","REACTOME_DNA_REPAIR","REACTOME_FANCONI_ANEMIA_PATHWAY","REACTOME_HDR_THROUGH_HOMOLOGOUS_RECOMBINATION_HRR","REACTOME_HDR_THROUGH_SINGLE_STRAND_ANNEALING_SSA","REACTOME_HOMOLOGOUS_DNA_PAIRING_AND_STRAND_EXCHANGE","REACTOME_HOMOLOGY_DIRECTED_REPAIR","REACTOME_RESOLUTION_OF_D_LOOP_STRUCTURES","REACTOME_RESOLUTION_OF_D_LOOP_STRUCTURES_THROUGH_SYNTHESIS_DEPENDENT_STRAND_ANNEALING_SDSA")

repairIV = c("KEGG_HOMOLOGOUS_RECOMBINATION","REACTOME_DISEASES_OF_DNA_REPAIR","REACTOME_DISEASES_OF_MISMATCH_REPAIR_MMR","REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR","REACTOME_DNA_REPAIR","REACTOME_FANCONI_ANEMIA_PATHWAY","REACTOME_HDR_THROUGH_HOMOLOGOUS_RECOMBINATION_HRR","REACTOME_HDR_THROUGH_SINGLE_STRAND_ANNEALING_SSA","REACTOME_HOMOLOGOUS_DNA_PAIRING_AND_STRAND_EXCHANGE","REACTOME_HOMOLOGY_DIRECTED_REPAIR","REACTOME_RESOLUTION_OF_D_LOOP_STRUCTURES","REACTOME_RESOLUTION_OF_D_LOOP_STRUCTURES_THROUGH_SYNTHESIS_DEPENDENT_STRAND_ANNEALING_SDSA")

repairVsH %in% repairIV
repairIV %in% repairVsH


RnavVsH = c(
  "REACTOME_TRANSPORT_OF_THE_SLBP_DEPENDANT_MATURE_MRNA",
  "REACTOME_RESPONSE_OF_MTB_TO_PHAGOCYTOSIS",
  "REACTOME_LEISHMANIA_INFECTION",
  "REACTOME_NUCLEAR_ENVELOPE_NE_REASSEMBLY",
  "REACTOME_MRNA_SPLICING_MINOR_PATHWAY",
  "REACTOME_FORMATION_OF_THE_EARLY_ELONGATION_COMPLEX",
  "REACTOME_NEUROTOXICITY_OF_CLOSTRIDIUM_TOXINS",
  "REACTOME_HCMV_INFECTION",
  "REACTOME_FORMATION_OF_TC_NER_PRE_INCISION_COMPLEX",
  "REACTOME_FCERI_MEDIATED_NF_KB_ACTIVATION",
  "REACTOME_HIV_TRANSCRIPTION_ELONGATION",
  "REACTOME_METABOLISM_OF_RNA",
  "REACTOME_PROCESSING_OF_CAPPED_INTRON_CONTAINING_PRE_MRNA",
  "REACTOME_HIV_LIFE_CYCLE",
  "REACTOME_INFLUENZA_INFECTION",
  "REACTOME_HCMV_LATE_EVENTS",
  "REACTOME_HOST_INTERACTIONS_OF_HIV_FACTORS",
  "REACTOME_SNRNP_ASSEMBLY",
  "REACTOME_FORMATION_OF_RNA_POL_II_ELONGATION_COMPLEX",
  "REACTOME_TRANSPORT_OF_MATURE_MRNAS_DERIVED_FROM_INTRONLESS_TRANSCRIPTS",
  "REACTOME_SUPPRESSION_OF_PHAGOSOMAL_MATURATION",
  "REACTOME_RNA_POLYMERASE_II_PRE_TRANSCRIPTION_EVENTS",
  "REACTOME_TOXICITY_OF_BOTULINUM_TOXIN_TYPE_D_BOTD",
  "REACTOME_MRNA_SPLICING",
  "REACTOME_INFECTION_WITH_MYCOBACTERIUM_TUBERCULOSIS",
  "REACTOME_SUMOYLATION_OF_CHROMATIN_ORGANIZATION_PROTEINS",
  "REACTOME_HCMV_EARLY_EVENTS",
  "REACTOME_PREVENTION_OF_PHAGOSOMAL_LYSOSOMAL_FUSION",
  "REACTOME_HIV_INFECTION",
  "REACTOME_RNA_POLYMERASE_II_TRANSCRIPTION_TERMINATION",
  "REACTOME_SLBP_DEPENDENT_PROCESSING_OF_REPLICATION_DEPENDENT_HISTONE_PRE_MRNAS",
  "REACTOME_PROCESSING_OF_CAPPED_INTRONLESS_PRE_MRNA",
  "REACTOME_INFECTIOUS_DISEASE",
  "REACTOME_GLYCOLYSIS",
  "REACTOME_PROCESSING_OF_INTRONLESS_PRE_MRNAS",
  "KEGG_SPLICEOSOME",
  "REACTOME_ABORTIVE_ELONGATION_OF_HIV_1_TRANSCRIPT_IN_THE_ABSENCE_OF_TAT",
  "REACTOME_TRANSPORT_OF_MATURE_TRANSCRIPT_TO_CYTOPLASM"
)

RnavIV = c(
  "REACTOME_FORMATION_OF_RNA_POL_II_ELONGATION_COMPLEX",
  "REACTOME_HIV_INFECTION",
  "REACTOME_INFECTIOUS_DISEASE",
  "REACTOME_POTENTIAL_THERAPEUTICS_FOR_SARS",
  "REACTOME_RNA_POLYMERASE_II_PRE_TRANSCRIPTION_EVENTS",
  "REACTOME_RNA_POLYMERASE_II_TRANSCRIPTION",
  "REACTOME_SARS_COV_INFECTIONS",
  "REACTOME_TRANSCRIPTION_OF_THE_HIV_GENOME"
)


RnavIV[!RnavIV %in% RnavVsH]

# ------------------------------------------------------------------
# --- Update Enrichments ----- 
# ------------------------------------------------------------------

source("src/reu/RERConvergeFunctions.R")


getStat = function(res){
  stat=sign(res$Rho)*(-log10(res$P))
  names(stat)=rownames(res)
  #deal with duplicated genes
  genenames=sub("\\..*", "",names(stat))
  multname=names(which(table(genenames)>1))
  for(n in multname){
    ii=which(genenames==n)
    iimax=which(max(stat[ii])==max(abs(stat[ii])))
    stat[ii[-iimax]]=NA
  }
  sum(is.na(stat))
  stat=stat[!is.na(stat)]
  
  stat
}


fastwilcoxGMTall = function (vals, annotList, alternative = "two.sided", ...) 
{
  reslist = list()
  for (n in names(annotList)) {
    reslist[[n]] = fastwilcoxGMT(vals, annotList[[n]], alternative = alternative, 
                                 ...)
    message(paste0(nrow(reslist[[n]]), " results for annotation set ", 
                   n))
  }
  reslist
}


vals = rerStats
gmt = gmtAnnotations

fastwilcoxGMT=function(vals, gmt, simple=T, use.all=F, num.g=10,genes=NULL, outputGeneVals=F, order=F,
                       alternative = "two.sided"){
  vals=vals[!is.na(vals)]
  if(is.null(genes)){
    genes=unique(unlist(gmt$genesets))
  }
  out=matrix(nrow=length(gmt$genesets), ncol=5)
  rownames(out)=gmt$geneset.names
  colnames(out)=c("stat", "pval", "p.adj","num.genes", "gene.vals")
  out=as.data.frame(out)
  genes=intersect(genes, names(vals))
  
  valsr=rank(vals[genes])
  numg=length(vals)+1
  valsallr=rank(vals)
  for( i in 1:nrow(out)){
    
    curgenes=intersect(genes,gmt$genesets[[i]])
    
    bkgenes=setdiff(genes, curgenes)
    
    if (length(bkgenes)==0 || use.all){
      bkgenes=setdiff(names(vals), curgenes)
    }
    if(length(curgenes)>=num.g & length(bkgenes)>2){
      if(!simple){
        # change alternative = "greater" for the one-sided test
        res=wilcox.test(x = vals[curgenes], y=vals[bkgenes], exact=F, alternative = alternative)
        
        out[i, 1:2]=c(res$statistic/(as.numeric(length(bkgenes))*as.numeric(length(curgenes))), res$p.value)
      }
      else{
        # add an alternative parameter (can be "greater" or "two.sided")
        out[i, 1:2]=simpleAUCgenesRanks(valsr[curgenes],valsr[bkgenes], alt = alternative)
        
      }
      
      out[i,"num.genes"]=length(curgenes)
      if(outputGeneVals){
        if (out[i,1]>0.5){
          oo=order(vals[curgenes], decreasing = T)
          granks=numg-valsallr[curgenes]
        }
        else{
          oo=order(vals[curgenes], decreasing = F)
          granks=valsallr[curgenes]
        }
        
        
        nn=paste(curgenes[oo],round((granks[curgenes])[oo],2),sep=':' )
        out[i,"gene.vals"]=paste(nn, collapse = ", ")
      }
    }
    
  }
  # hist(out[,2])
  out[,1]=out[,1]-0.5
  out[, "p.adj"]=p.adjust(out[,2], method="BH")
  
  out=out[!is.na(out[,2]),]
  if(order){
    out=out[order(-abs(out[,1])),]
  }
  out
}


pos = valsr[curgenes]
neg = valsr[bkgenes]

simpleAUCgenesRanks=function(pos, neg, alt = "two.sided"){
  
  posn=length(pos)
  negn=length(neg)
  posn=as.numeric(posn)
  negn=as.numeric(negn)
  stat=sum(pos)-posn*(posn+1)/2 #Average pos
  auc=stat/(posn*negn)
  mu=posn*negn/2
  sd=sqrt((posn*negn*(posn+negn+1))/12)
  
  if(alt == "two.sided") {
    stattest=apply(cbind(stat, posn*negn-stat),1,max)
    pp=(2*pnorm(stattest, mu, sd, lower.tail = F))
  }
  
  else if(alt == "greater"){
    pp=(pnorm(stat,mu,sd,lower.tail=FALSE)) 
  }
  return(c(auc,pp))
}


# ------------------------------------------------------------------
# --- Make liam-non-liam comparison plots ----- 
# ------------------------------------------------------------------


# -- argument setup  -- 
significanceCutoff = 0.05
prefix = "CategoricalInsvertivoreTreeLiamInference"
pairwiseSets = c("Herbivore-Insectivore", "Herbivore-Vertivore", "Carnivore-Herbivore", "Herbivore-Omnivore", "Insectivore-Vertivore", "Omnivore-Vertivore", "Invertivore-Omnivore")
geneSet = "KeggReactome"
vennDiagramSet = c("Herbivore-Invertivore", "Herbivore-Vertivore", "Carnivore-Herbivore")  
vennColorset = c("darkblue", "red", "orange")
usingGo = !is.null(geneSet)
saveCombinedData = T
saveCombinedData = F
bothAxis = T
saveData = T
saveData = F

args = c("r=CategoricalInsvertivoreTree")
args = c("r=CategoricalInsvertivoreTreeLiamInference")


nonliamResults = combinedResults
nonliamGoResults = GoCombinedResults

liamResults = combinedResults
liamGoResults = GoCombinedResults

# -- Read Data 
combinedGeneDataFilename = paste0(outputFolderName, filePrefix, "combinedGeneResults.rds")
combinedResults = readRDS(combinedGeneDataFilename)
significanceColumns = names(combinedResults)[grep("significant", names(combinedResults))]
geneSignificanceResults = combinedResults[, names(combinedResults) %in% significanceColumns]


combinedGODataFilename = paste0(outputFolderName, filePrefix, "combinedGOResults-", geneSet, ".rds")
GoCombinedResults = readRDS(combinedGODataFilename)
GoSignificanceColumns = names(GoCombinedResults)[grep("significant", names(GoCombinedResults))]
GoSignificanceResults = GoCombinedResults[, names(GoCombinedResults) %in% GoSignificanceColumns]







# -- make resources to prefix-phenotype conversion 
prefixSet = NULL
prefixList = NULL
for(i in 1:length(pairwiseSets)){
  currentSet = pairwiseSets[i]
  correlationSubsetName = gsub("-", " - ", currentSet)
  correlationPrefix = paste(substr(strsplit(currentSet, split = "-")[[1]],1,1), collapse = '')
  prefixEntry = correlationPrefix; names(prefixEntry) = currentSet; prefixSet = append(prefixSet, prefixEntry)
  for(j in 1:2){
    prefixSingle = strsplit(correlationPrefix, split = "")[[1]][j]; names(prefixSingle) = strsplit(currentSet, split="-")[[1]][j]; prefixList = append(prefixList, prefixSingle)
  }
}


prefixList = unlist(prefixList); prefixList = prefixList[!duplicated(prefixList)]

# Debug code line for personal use 
names(prefixList)[which(names(prefixList) == "Insectivore")] = "Invertivore"

replacePrefixWithName = function(x) {
  inversePrefixList = setNames(names(prefixList), prefixList); parts = unlist(strsplit(x, "-")); fullNames = inversePrefixList[parts]; paste(fullNames, collapse = "-")
}

addDashes = function(vector) {
  sapply(vector, function(s) paste(strsplit(s, "")[[1]], collapse = "-"))
}


#-------------------------------------------------------------------



#gene
liamResults = liamGeneResults
nonliamResults = nonliamGeneResults

editedLiamResults = liamResults
colnames(editedLiamResults) = paste0(colnames(liamResults), "-Liam")
metaCombinedResults = cbind(nonliamResults, editedLiamResults)

combinedResults = metaCombinedResults
geneMetaCombined = metaCombinedResults


#GO
#nonliamGeneResults = nonliamResults
#liamGeneResults = liamResults

liamResults = liamGoResults
nonliamResults = nonliamGoResults

editedLiamResults = liamResults
colnames(editedLiamResults) = paste0(colnames(liamResults), "-Liam")
metaCombinedResults = cbind(nonliamResults, editedLiamResults)

names(combinedResults) = gsub("-stat", corrleationColumnType, names(combinedResults))

combinedResults = metaCombinedResults
GoMetaCombined = metaCombinedResults



#Prefix
prefixList = append(prefixList, c("L", "i", "a", "m", "-", ""))
names(prefixList) = append(names(prefixList)[1:5], c("Liam", "", "", "", "", ""))

bothAxis = F

corrleationColumnType = "-p.adj"
#-------------------------------------------------------------------

{
  
  
  
  grep(corrleationColumnType, names(combinedResults))
  rhoValues = combinedResults[,grep(corrleationColumnType, names(combinedResults))]
  
  rhoComparisions = names(rhoValues)
  
  #rhoComparisions = rhoComparisions[-c(5,6,7,9,13,14,15)]
  
  rhoPhenotypes = strsplit(gsub(corrleationColumnType, "", rhoComparisions), split = "")
  commonBackground = Reduce(intersect, rhoPhenotypes)
  if(length(commonBackground) == 1){
    for(i in 1:length(rhoPhenotypes)){ #invert tho if background in second postion so rho has consistent meaning relative to background
      if(rhoPhenotypes[[i]][1] != commonBackground){
        cat("Inverting rho of ", rhoPhenotypes[[i]] , "becuase background is in first position.")
        rhoValues[i] = -1*rhoValues[i]
      }
    }
  }
  
  densityScaleSet = NULL
  
  #generate the plots slim-ly to get the desired density scale 
  for(i in 1:length(rhoValues)){
    xName = names(rhoValues)[i]
    if(i+1 <= length(rhoValues)){
      for(j in (i+1):length(rhoValues)){
        yName = names(rhoValues)[j]
        
        rhoCorrellPlot = ggplot(rhoValues, aes(x = .data[[xName]], y = .data[[yName]])) + 
          geom_point() + geom_pointdensity() + scale_color_viridis()
        
        
        denstiyScaleValue = ggplot_build(rhoCorrellPlot)$plot$scales$scales[[1]]$get_limits()[2]
        densityScaleSet = append(densityScaleSet, denstiyScaleValue)
        rm(rhoCorrellPlot)
      }
    }
  }
  densityScale = c(1, max(densityScaleSet))
  
  
  rhoPlotSet = list()
  netIndex= 0
  for(i in 1:length(rhoValues)){
    xName = names(rhoValues)[i]
    if(i <= length(rhoValues)){
      if(bothAxis){jStart = 1}else{jStart = i+1}
      for(j in (jStart):length(rhoValues)){
        yName = names(rhoValues)[j]
        yLabel =  paste0(replacePrefixWithName(addDashes(gsub("-", "", gsub(corrleationColumnType, "", yName)))), " Dunn Z Statistic")
        xLabel =  paste0(replacePrefixWithName(addDashes(gsub("-", "", gsub(corrleationColumnType, "", xName)))), " Dunn Z Statistic")
        
        rhoCorrellPlot = ggplot(rhoValues, aes(x = .data[[xName]], y = .data[[yName]])) + 
          geom_point() + geom_pointdensity() + scale_color_viridis(name = "Density of genes", limits = densityScale) + 
          stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),formula = y ~ x,parse = TRUE, size = 6) +
          theme_classic()+
          xlab(xLabel) + ylab(yLabel)+
          theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))
        
        netIndex = netIndex +1
        rhoPlotSet[[netIndex]] = rhoCorrellPlot
        names(rhoPlotSet)[netIndex] = paste(xName, yName, sep="-")
        rm(rhoCorrellPlot)
      }
    }
  }
  # add a null plot by using two using permulations as a comparision 
  
  #pdf()
  #print(rhoPlotSet)
  #dev.off()
}

names(rhoPlotSet)

firstPhen = sapply(strsplit(names(rhoPlotSet), "-"), `[`, 1)
secondPhen = sapply(strsplit(names(rhoPlotSet), "-"), `[`, 3)

matchingPlots = which(firstPhen == secondPhen)

liamComaprePlots = rhoPlotSet[matchingPlots]
names(liamComaprePlots)

#Gene
CHPlot = liamComaprePlots[1]
CHPlot = CHPlot[[1]]

HIPlot = liamComaprePlots[2]
HIPlot = HIPlot[[1]]

HVPlot = liamComaprePlots[4]
HVPlot = HVPlot[[1]]


herbivoreCompare = grid.arrange(CHPlot, HIPlot, HVPlot, nrow = 2)

#geneComparePlot = herbivoreCompare
genePvalComparePlot = herbivoreCompare



#GO
CHPlot = liamComaprePlots[3]
CHPlot = CHPlot[[1]]

HIPlot = liamComaprePlots[1]
HIPlot = HIPlot[[1]]

HVPlot = liamComaprePlots[2]
HVPlot = HVPlot[[1]]


herbivoreCompare = grid.arrange(CHPlot, HIPlot, HVPlot, nrow = 2)

#GoComparePlot = herbivoreCompare
GoPvalComparePlot = herbivoreCompare


plot(geneComparePlot)
plot(GoComparePlot)


# ------------------------------------------------------------------
# ---------- Looking into the specific lost GO categories 

GoMetaCombined
names(GoMetaCombined)
which(GoMetaCombined$`HI-HV-Overlap` & !GoMetaCombined$`HI-HV-CH-Overlap`)

# --- Finding changed pathways 

HiDifference = which(GoMetaCombined$`HI-significant` != GoMetaCombined$`HI-significant-Liam`)
HvDifference = which(GoMetaCombined$`HV-significant` != GoMetaCombined$`HV-significant-Liam`)
ChDifference = which(GoMetaCombined$`CH-significant` != GoMetaCombined$`CH-significant-Liam`)






nonliamOnlyHI = length(which(GoMetaCombined[HiDifference, c(which(colnames(GoMetaCombined) == "HI-significant"),which(colnames(GoMetaCombined) == "HI-significant-Liam"))][,1]))
liamOnlyHI = length(which(GoMetaCombined[HiDifference, c(which(colnames(GoMetaCombined) == "HI-significant"),which(colnames(GoMetaCombined) == "HI-significant-Liam"))][,2]))
HiLoss = nonliamOnlyHI-liamOnlyHI
nonliamOnlyHV = length(which(GoMetaCombined[HiDifference, c(which(colnames(GoMetaCombined) == "HV-significant"),which(colnames(GoMetaCombined) == "HV-significant-Liam"))][,1]))
liamOnlyHV = length(which(GoMetaCombined[HiDifference, c(which(colnames(GoMetaCombined) == "HV-significant"),which(colnames(GoMetaCombined) == "HV-significant-Liam"))][,2]))
HvLoss = nonliamOnlyHV-liamOnlyHV
nonliamOnlyCH = length(which(GoMetaCombined[HiDifference, c(which(colnames(GoMetaCombined) == "CH-significant"),which(colnames(GoMetaCombined) == "CH-significant-Liam"))][,1]))
liamOnlyCH = length(which(GoMetaCombined[HiDifference, c(which(colnames(GoMetaCombined) == "CH-significant"),which(colnames(GoMetaCombined) == "CH-significant-Liam"))][,2]))
ChLoss = nonliamOnlyCH-liamOnlyCH

{
cat(paste0("Number of nonliam-only HI: \n", length(which(GoMetaCombined[HiDifference, c(which(colnames(GoMetaCombined) == "HI-significant"),which(colnames(GoMetaCombined) == "HI-significant-Liam"))][,1])), "\n"))
cat(paste0("Number of liam-only HI: \n", length(which(GoMetaCombined[HiDifference, c(which(colnames(GoMetaCombined) == "HI-significant"),which(colnames(GoMetaCombined) == "HI-significant-Liam"))][,2])), "\n"))
cat(paste0("Number of nonliam-only HV: \n", length(which(GoMetaCombined[HiDifference, c(which(colnames(GoMetaCombined) == "HV-significant"),which(colnames(GoMetaCombined) == "HV-significant-Liam"))][,1])), "\n"))
cat(paste0("Number of liam-only HV: \n", length(which(GoMetaCombined[HiDifference, c(which(colnames(GoMetaCombined) == "HV-significant"),which(colnames(GoMetaCombined) == "HV-significant-Liam"))][,2])), "\n"))
cat(paste0("Number of nonliam-only CH: \n", length(which(GoMetaCombined[HiDifference, c(which(colnames(GoMetaCombined) == "CH-significant"),which(colnames(GoMetaCombined) == "CH-significant-Liam"))][,1])), "\n"))
cat(paste0("Number of liam-only CH: \n", length(which(GoMetaCombined[HiDifference, c(which(colnames(GoMetaCombined) == "CH-significant"),which(colnames(GoMetaCombined) == "CH-significant-Liam"))][,2])), "\n"))
}

nonliamSignificant = (which(GoMetaCombined$`HI-significant` | GoMetaCombined$`HV-significant` | GoMetaCombined$`CH-significant`))
liamSignificant = (which(GoMetaCombined$`HI-significant-Liam` | GoMetaCombined$`HV-significant-Liam` | GoMetaCombined$`CH-significant-Liam`))

length(which(liamSignificant %in% nonliamSignificant))
#260, so almost all of the ones in liam are in the nonliam 

liamMissing = nonliamSignificant[which(!nonliamSignificant %in% liamSignificant)]
length(liamMissing)


liamMissingGO = GoMetaCombined[liamMissing,]
liamMissingGO = liamMissingGO[,-grep("Overlap", colnames(liamMissingGO))]
liamMissingGO = liamMissingGO[,-grep("Driver", colnames(liamMissingGO))]
liamMissingGO = liamMissingGO[,-grep("O", colnames(liamMissingGO))]

write.csv(liamMissingGO, "Results/LiamMissingGO.csv")

# --- Opposite pathways 
oppsitePathways = GoMetaCombined[which(GoMetaCombined$`HI-HV-Overlap` & !GoMetaCombined$`HI-HV-CH-Overlap`),]

grep("Overlap", colnames(oppsitePathways))

oppsitePathways = oppsitePathways[,-grep("Overlap", colnames(oppsitePathways))]
oppsitePathways = oppsitePathways[,-grep("Driver", colnames(oppsitePathways))]
oppsitePathways = oppsitePathways[,-grep("O", colnames(oppsitePathways))]


# ------------------------------------------------------------------
# --- Figure out paths files ----- 
# ------------------------------------------------------------------

OgPaths = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalPathsFile.rds")

#Make thoeretical non-liam paths

modelType = "ER"
ancestralTrait = NULL

mainTrees = readRDS("Data/zoonomiaAllMammalsTrees.rds")

nonliamPhenotype = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalPhenotypeVector.rds")
nonliamTree = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalTree.rds")
nonliamSpeciesFilter = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeSpeciesFilter.rds")

nonliamCharpaths = char2PathsCategorical(nonliamPhenotype, mainTrees, nonliamSpeciesFilter, model = modelType, anctrait = ancestralTrait) #make a path based on the phenotype vector
nonliamTreePaths = tree2Paths(nonliamTree, mainTrees, useSpecies = nonliamSpeciesFilter, categorical = TRUE) #do not binarize; the categorical data is already contained in the phenotype tree.

all.equal(OgPaths, nonliamCharpaths)
all.equal(OgPaths, nonliamTreePaths)

LiamInference

ogliamPaths = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCategoricalPathsFile.rds")


liamPhenotype = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCategoricalPhenotypeVector.rds")
liamTree = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCategoricalTree.rds")
liamSpeciesFilter = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceSpeciesFilter.rds")

liamCharpaths = char2PathsCategorical(liamPhenotype, mainTrees, liamSpeciesFilter, model = modelType, anctrait = ancestralTrait) #make a path based on the phenotype vector
liamTreePaths = tree2Paths(liamTree, mainTrees, useSpecies = liamSpeciesFilter, categorical = TRUE) #do not binarize; the categorical data is already contained in the phenotype tree.


all.equal(ogliamPaths, liamCharpaths)
all.equal(ogliamPaths, liamTreePaths)


# ------------------------------------------------------------------
# --- Compare GO Significance thresholds ----- 
# ------------------------------------------------------------------
combinedVenn = grid.arrange(geneVenn, goVenn, nrow = 1)

pdf("Output/CategoricalInsVertivoreTreeLiamInference/Categproca;")

vennPlotName = paste0(outputFolderName, filePrefix, "VennPlots.pdf")
overlapPlotName = paste0(outputFolderName, filePrefix, "OverlapPlots.pdf")



pdf(vennPlotName, 30,15)
grid.arrange(combinedVenn)
dev.off()



saveRDS(GoSignificanceColumns, paste0(combinedGODataFilename, ".rds"))


for(i in 1:length(pairwiseSets)){
  currentSet = pairwiseSets[i]
  correlationSubsetName = gsub("-", " - ", currentSet)
  correlationPrefix = paste(substr(strsplit(currentSet, split = "-")[[1]],1,1), collapse = '')
  
  goFilename = paste0(outputFolderName, currentSet, "/", filePrefix, currentSet,"Enrichment-", geneSet, ".rds")
  currentGoData = readRDS(goFilename)[[1]]
  currentGoData$significant = currentGoData$p.adj < 0.05
  message(length(which(currentGoData$significant)))
  currentGoData$significant = currentGoData$p.adj < 0.1
  message(length(which(currentGoData$significant)))
  
  names(currentGoData) = paste0(correlationPrefix, "-", names(currentGoData))
  
  
  driverFilename = paste0(outputFolderName, currentSet, "/", filePrefix, currentSet, "GoDriverTable-", geneSet, ".rds")
  if(file.exists(driverFilename)){
    driverTable = readRDS(driverFilename)
    if(all(rownames(driverTable) == rownames(currentGoData))){
      currentGoData$Driver = driverTable[which(names(driverTable) == "Driver")][[1]]
      currentGoData$DriverNumeric = driverTable[which(names(driverTable) == "DriverNumeric")][[1]]
      
      names(currentGoData)[which(names(currentGoData) == "Driver")] = paste0(correlationPrefix, "-", "Driver")
      names(currentGoData)[which(names(currentGoData) == "DriverNumeric")] = paste0(correlationPrefix, "-", "DriverNumeric")    
      
    }
    rm(driverTable)
  }
  GOResults[[i]] = currentGoData
  names(GOResults)[i] = correlationPrefix
  
}



# ------------------------------------------------------------------
# --- Compare liam results and non-liam results ----- 
# ------------------------------------------------------------------

saveRDS(GoCombinedResults, "Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsvertivoreTreeLiamInferencecombinedGOResults-KeggReactome.rds")

noLiamGenes = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsvertivoreTreecombinedGeneResults.rds")
liamGenes = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsvertivoreTreeLiamInferencecombinedGeneResults.rds")

noLiamGO = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsvertivoreTreecombinedGOResults-KeggReactome.rds")
liamGO = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsvertivoreTreeLiamInferencecombinedGOResults-KeggReactome.rds")

liam = liamGenes
nonliam = noLiamGenes
testCol = "HI-p.adj"

compareOverlap = function(testCol, liam, nonliam){
  liamColValue = which(colnames(liam) == testCol)
  liamCol = liam[liamColValue]
  
  nonliamColValue = which(colnames(nonliam) == testCol)
  nonliamCol = nonliam[nonliamColValue]
  
  isPval = !is.logical(liamCol[[1]])
  
  if(isPval){
    plot(liamCol[[1]], nonliamCol[[1]])
    
    model = lm(nonliamCol[[1]] ~ liamCol[[1]])
    r2 = summary(model)$r.squared
    text(x = min(liamCol[[1]], na.rm = T), y = max(nonliamCol[[1]], na.rm = T), labels = paste("R² =", round(r2, 3)), pos = 4)
    
    liamCol[[1]] = liamCol[[1]] < 0.1
    nonliamCol[[1]] = nonliamCol[[1]] < 0.1
  }
  
  
  
  liamSignificant = length(which(liamCol[[1]]))
  nonliamSignificant = length(which(nonliamCol[[1]]))
  
  liamInNonliam = length(which(which(liamCol[[1]]) %in% which(nonliamCol[[1]])))
  liamOverlapPercent = liamInNonliam/liamSignificant
  
  nonliamInLiam = length(which(which(nonliamCol[[1]]) %in% which(liamCol[[1]])))
  nonliamOverlapPercent = nonliamInLiam/nonliamSignificant
  
  cat("##########\n")
  cat(testCol)
  cat("\n")
  if(isPval){
    cat("R Squared: ")
    cat(r2)
    cat("\n")
  }
  
  #cat("#\n")
  cat("Liam value: ")
  cat(liamSignificant)
  cat("\n")
  cat("Non-Liam value: ")
  cat(nonliamSignificant)
  cat("\n")
  cat("Liam non-Liam ratio: ")
  cat(nonliamSignificant/liamSignificant)
  cat("\n")
  
  #cat("# \n")
  cat("Liam in Non-Liam: ")
  cat(liamInNonliam)
  cat("     ")
  cat(liamOverlapPercent)
  cat("\n")
  
  cat("Non-Liam in Liam: ")
  cat(nonliamInLiam)
  cat("     ")
  cat(nonliamOverlapPercent)
  cat("\n")
  
}

compareOverlap("HI-significant", liamGenes, noLiamGenes)
compareOverlap("HV-significant", liamGenes, noLiamGenes)
compareOverlap("IV-significant", liamGenes, noLiamGenes)
compareOverlap("CH-significant", liamGenes, noLiamGenes)

compareOverlap("HI-p.adj", liamGenes, noLiamGenes)
compareOverlap("HV-p.adj", liamGenes, noLiamGenes)
compareOverlap("IV-p.adj", liamGenes, noLiamGenes)
compareOverlap("CH-p.adj", liamGenes, noLiamGenes)

# MInimum overlap is 74% in HV, or 69% in CH, liam is larger than non-liam. 

compareOverlap("HI-significant", liamGO, noLiamGO)
compareOverlap("HV-significant", liamGO, noLiamGO)
compareOverlap("IV-significant", liamGO, noLiamGO)
compareOverlap("CH-significant", liamGO, noLiamGO)


compareOverlap("HI-p.adj", liamGO, noLiamGO)
compareOverlap("HV-p.adj", liamGO, noLiamGO)
compareOverlap("IV-p.adj", liamGO, noLiamGO)
compareOverlap("CH-p.adj", liamGO, noLiamGO)


# ------------------------------------------------------------------
# --- I-V numbers ----- 
# ------------------------------------------------------------------

genesResults = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsvertivoreTreecombinedGeneResults.rds")

length(which(genesResults$`IV-significant`))

GoResults = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsvertivoreTreecombinedGOResults-KeggReactome.rds")
length(which(GoResults$`IV-significant`))


# ------------------------------------------------------------------
# --- Examine liam results ----- 
# ------------------------------------------------------------------

liamResultsCore = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCombinedCategoricalCorrelationFile.rds")
omnibus = liamResultsCore[[1]]
length(which(omnibus$p.adj < 0.05))

liamGoOverall = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/Overall/CategoricalInsVertivoreTreeLiamInferenceOverallEnrichment-KeggReactome.rds")

liamGoOverall = liamGoOverall[[1]]
length(which(liamGoOverall$p.adj < 0.05))

liamCarnResults = readRDS("Output/CategoricalInsVertivoreTreeCarnivoreLiamInference/CategoricalInsvertivoreTreeCarnivoreLiamInferencecombinedGeneResults.rds")

length(which(liamCarnResults$`CH-significant`))
length(which(liamCarnResults$`CO-significant`))
length(which(liamCarnResults$`HO-significant`))

liamResultsCoreCarn = readRDS("Output/CategoricalInsVertivoreTreeCarnivoreLiamInference/CategoricalInsVertivoreTreeCarnivoreLiamInferenceCombinedCategoricalCorrelationFile.rds")
omnibusCarn = liamResultsCoreCarn[[1]]
length(which(omnibusCarn$p.adj < 0.05))


liamCarnGoResults = readRDS("Output/CategoricalInsVertivoreTreeCarnivoreLiamInference/CategoricalInsvertivoreTreeCarnivoreLiamInferencecombinedGOResults-KeggReactome.rds")
length(which(liamCarnGoResults$`CH-significant`))
length(which(liamCarnGoResults$`CO-significant`))
length(which(liamCarnGoResults$`HO-significant`))


liamCarnGoOverall = readRDS("Output/CategoricalInsVertivoreTreeCarnivoreLiamInference/Overall/CategoricalInsVertivoreTreeCarnivoreLiamInferenceOverallEnrichment-KeggReactome.rds")

liamCarnGoOverall = liamCarnGoOverall[[1]]
length(which(liamCarnGoOverall$p.adj < 0.05))




carnGoResults = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsvertivoreTreecombinedGOResults-KeggReactome.rds")
length(which(carnGoResults$`CH-significant`))
length(which(carnGoResults$`CO-significant`))
length(which(carnGoResults$`HO-significant`))

noliamResultsCore = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsvertivoreTreecombinedGeneResults.rds")
length(which(noliamResultsCore$`IV-significant`))

# Add carnivory results to main liam infrence

`Carnivore - Herbivore` = liamResultsCoreCarn[[2]][1] 
names(`Carnivore - Herbivore`) = "Carnivore - Herbivore"
`Carnivore - Omnivore` = liamResultsCoreCarn[[2]][2] 
names(`Carnivore - Omnivore`) = "Carnivore - Omnivore"

correlationResults[7] = `Carnivore - Herbivore`
names(correlationResults)[7] = "Carnivore - Herbivore"
correlationResults[8] = `Carnivore - Omnivore`
names(correlationResults)[8] = "Carnivore - Omnivore"

saveRDS(correlationResults, "Output/CategoricalInsvertivoreTreeLiamInference/CategoricalInsvertivoreTreeLiamInferencePairwiseCorrelationFile.rds")
write.csv(correlationResults, "Output/CategoricalInsvertivoreTreeLiamInference/CategoricalInsvertivoreTreeLiamInferencePairwiseCorrelationFile.csv")
# ------------------------------------------------------------------
# --- Make RER PLots ----- 
# ------------------------------------------------------------------

mainTrees = readRDS("Data/zoonomiaAllMammalsTrees.rds")
pathObject = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalPathsFile.rds")
RERObject = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeRERFile.rds")


commonRERs = RERObject
colnames(commonRERs) = ZonomNameConvertVectorCommon(colnames(commonRERs), tipColumn = "ZoonomiaTip")
palette(c( "darkgreen", "darkblue","black", "red"))

postiveSelectionOverlap = c("CTRL", "CPA1", "SLC36A1", "SLC6A19", "CELA3B", "CLPS", "SLC7A9", "PNLIPRP2", "FABP1", "ABCG8", "LCT")

postiveSelectionNoOverlap = c("APOB", "APOA1", "LIPF","NPC1L1","MEP1B","SLC3A2","HK1","MGAM", 
                        "APOB", "PLPP2", "APOA1", "SLC27A4", "PLA2G1B", "SLC8A2", "SLC3A1", "DPP4", "KCNN4", "SLC3A2", "MGAM2", "HK3", "G6PC",
                        "APOB", "PLA2G5", "SCAB1", "CPA3", 
                        "APOB", "MOGAT2", "MTTP", "SLC1A5", "SLC7A8", "SLC15A1", "PRKCB", "ATP1B1", "ATP1B3",
                        "SLC3A2", "SLC1A5",  "DPP4", "APOB", "CD36", "PLPP2", 
                        "APOB", "PIK3CD", "CPB2", "KCNK5"
)
length(unique(postiveSelectionNoOverlap))
length(unique(postiveSelectionOverlap))

unique(rownames(commonRERs))

i=1
{
plotRers(commonRERs, postiveSelectionOverlap[i], pathObject, sortrers = T)
i=i+1
}

i=1
{
  plotRers(commonRERs, postiveSelectionNoOverlap[i], pathObject, sortrers = T)
  i=i+1
}

i=1
{
  plotRers(commonRERs, unique(rownames(commonRERs))[i], pathObject, sortrers = T)
  i=i+1
}
#


plotRers(commonRERs, "RAD50", pathObject, sortrers = T)

# ------------------------------------------------------------------
# --- Make carnivory results for liam  ----- 
# ------------------------------------------------------------------

liamTree = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCategoricalTree.rds")
liamPhenotype = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCategoricalPhenotypeVector.rds")
liamFilter = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceSpeciesFilter.rds")


mergedData = read.csv("Data/mergedData.csv")
mixInvVertSpecies = which(mergedData$DerekDietClassification90InsVertivoreSorting == "C-InsVertivore-Mixed")

mixSpecies = mergedData$ZoonomiaTip[mixInvVertSpecies]

names(liamCarnivoryPhenotype)[which(names(liamCarnivoryPhenotype) %in% mixSpecies)]

ZonomNameConvertVectorCommon(names(liamCarnivoryPhenotype)[which(names(liamCarnivoryPhenotype) %in% mixSpecies)], tipColumn ="ZoonomiaTip")

liamCarnivoryPhenotype = liamPhenotype
liamCarnivoryPhenotype[which(liamCarnivoryPhenotype == "Vertivore" | liamCarnivoryPhenotype == "Insectivore")] = "Carnivore"
liamCarnivoryPhenotype[which(names(liamCarnivoryPhenotype) %in% mixSpecies)] = "Carnivore"

liamTreeCarnivory = liamTree
liamTreeCarnivory$edge.length[which(liamTreeCarnivory$edge.length %in% c(1))] = 5
liamTreeCarnivory$edge.length[which(liamTreeCarnivory$edge.length %in% c(3))] = 6
liamTreeCarnivory$edge.length[which(liamTreeCarnivory$edge.length %in% c(2,4))] = 1
liamTreeCarnivory$edge.length[which(liamTreeCarnivory$edge.length %in% c(5))] = 2
liamTreeCarnivory$edge.length[which(liamTreeCarnivory$edge.length %in% c(6))] = 3


speciesFilter = liamFilter
phenotypeVector = liamCarnivoryPhenotype
categoricalTree = liamTreeCarnivory
outputFolderName = "Output/CategoricalInsVertivoreTreeCarnivoreLiamInference/"
filePrefix = "CategoricalInsVertivoreTreeCarnivoreLiamInference"
phenotypeVectorFilename = paste(outputFolderName, filePrefix, "CategoricalPhenotypeVector.rds",sep="") #make a filename based on the prefix
saveRDS(liamCarnivoryPhenotype, file = phenotypeVectorFilename)                        #save the phenotype vector

speciesFilterFilename = paste(outputFolderName, filePrefix, "SpeciesFilter.rds",sep="") #set a filename for the species filter based on the prefix 
saveRDS(speciesFilter, file = speciesFilterFilename)                          #save that as the species filter




spreadSheetLocation = "Data/mergedData.csv"
nameColumn = "ZoonomiaTip"
mainTrees = readRDS("data/zoonomiaAllMammalsTrees.rds")
commonMainTrees = mainTrees
commonMainTrees$masterTree = ZoonomTreeNameToCommon(commonMainTrees$masterTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)
commonPhenotypeVector = phenotypeVector
names(commonPhenotypeVector) = ZonomNameConvertVectorCommon(names(commonPhenotypeVector), annotationLocation = spreadSheetLocation, tipColumn = nameColumn)
commonSpeciesFilter = ZonomNameConvertVectorCommon(speciesFilter, annotationLocation = spreadSheetLocation, tipColumn = nameColumn)
commonCategoricalTree = ZoonomTreeNameToCommon(categoricalTree, tipCol = "ZoonomiaTip")



treeImageFilename = paste(outputFolderName, filePrefix, "CategoricalTree.pdf", sep="") #make a filename based on the prefix
palette(c("red", "darkgreen", "black"))

pdf(treeImageFilename, height = length(phenotypeVector)/18, width = 10)                     #make a pdf to store the plot, sized based on tree size

plotTreeCategorical(commonCategoricalTree, c("Carnivore", "Herbivore", "Omnivore"), master = commonMainTrees$masterTree)
plotTreeCategorical(categoricalTree, c("Carnivore", "Herbivore", "Omnivore" ), master = mainTrees$masterTree)

#plotTreeCategorical(commonCategoricalTreeExtraTip, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = commonMasterAdded)
#plotTreeCategorical(categoricalTreeExtraTip, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = masterTreeAdded)
dev.off()  

categoricalTreeFilename = paste(outputFolderName, filePrefix, "CategoricalTree.rds", sep="") #make a filename based on the prefix
saveRDS(categoricalTree, categoricalTreeFilename)                               #save the tree
categoricalCommonTreeFilename = paste(outputFolderName, filePrefix, "CategoricalCommonTree.rds", sep="") #make a filename based on the prefix
saveRDS(commonCategoricalTree, categoricalCommonTreeFilename)


pathsFilename = paste(outputFolderName, filePrefix, "CategoricalPathsFile.rds", sep= "") #make a filename based on the prefix
paths = char2PathsCategorical(phenotypeVector, mainTrees, speciesFilter, model = modelType, anctrait = ancestralTrait) #make a path based on the phenotype vector
saveRDS(paths, file = pathsFilename)                                            #save the path 



# ------------------------------------------------------------------
# --- Emily CC figure ----- 
# ------------------------------------------------------------------

getPermsBinary=function(numperms, fg_vec, sisters_list, root_sp, RERmat, trees, mastertree, permmode="cc", method="k", min.pos=2, trees_list=NULL, calculateenrich=F, annotlist=NULL){
  pathvec = foreground2Paths(fg_vec, trees, clade="all",plotTree=F)
  col_labels = colnames(trees$paths)
  names(pathvec) = col_labels
  
  message("As of RERConverge [X.xx], permulation functions have been updated. Old versions have been moved to ccLegacy and ssmLegacy.")
  if(permmode=="cc"){
    print("Running CC permulation. sisters_list is required only for enrichments, otherwise sisters_list = NA is sufficient.")
    
    print("Generating permulated trees")
    
    # --- new code since legacy method; switching to categorical function to infer phenotype tree --
    #covert fg_vec to a categorical phenotypeVector
    phenotypeVector = rep(0, length(trees$masterTree$tip.label))
    names(phenotypeVector) = trees$masterTree$tip.label
    phenotypeVector[names(phenotypeVector) %in% fg_vec] = 1
    
    
    permulationData = categoricalPermulations(trees, phenotypeVector, rm = "ER", rp = "auto", ntrees = numperms)
    
    
    permulatedTrees = lapply(permulationData$trees, function(x) {
      tr = trees$masterTree
      tr$edge.length = c(x$tips, x$nodes)[tr$edge[,2]]
      tr$edge.length = tr$edge.length-1
      names(tr$edge.length) = NULL
      #tree2Paths(tr, treesObj, categorical = TRUE, useSpecies = names(phenvals))
      tr
    })
    permulated.binphens = list(permulatedTrees)
    #----
    #permulated.binphens = generatePermulatedBinPhen(trees$masterTree, numperms, trees, root_sp, fg_vec, sisters_list, pathvec, permmode="cc")
    permulated.fg = mapply(getForegroundsFromBinaryTree, permulated.binphens[[1]])
    permulated.fg.list = as.list(data.frame(permulated.fg))
    phenvec.table = mapply(foreground2Paths,permulated.fg.list,MoreArgs=list(treesObj=trees,clade="all"))
    phenvec.list = lapply(seq_len(ncol(phenvec.table)), function(i) phenvec.table[,i])
    
    print("Calculating correlations")
    corMatList = lapply(phenvec.list, correlateWithBinaryPhenotype, RERmat=RERmat)
    
    #make enrich list/matrices to fill
    permPvals=data.frame(matrix(ncol=numperms, nrow=nrow(RERmat)))
    rownames(permPvals)=rownames(RERmat)
    permRhovals=data.frame(matrix(ncol=numperms, nrow=nrow(RERmat)))
    rownames(permRhovals)=rownames(RERmat)
    permStatvals=data.frame(matrix(ncol=numperms, nrow=nrow(RERmat)))
    rownames(permStatvals)=rownames(RERmat)
    
    for (i in 1:length(corMatList)){
      permPvals[,i] = corMatList[[i]]$P
      permRhovals[,i] = corMatList[[i]]$Rho
      permStatvals[,i] = sign(corMatList[[i]]$Rho)*-log10(corMatList[[i]]$P)
    }
    
  }
  else if (permmode=="ccLegacy"){
    print("Running CC Legacy permulation")
    
    print("Generating permulated trees")
    permulated.binphens = generatePermulatedBinPhen(trees$masterTree, numperms, trees, root_sp, fg_vec, sisters_list, pathvec, permmode="cc")
    permulated.fg = mapply(getForegroundsFromBinaryTree, permulated.binphens[[1]])
    permulated.fg.list = as.list(data.frame(permulated.fg))
    phenvec.table = mapply(foreground2Paths,permulated.fg.list,MoreArgs=list(treesObj=trees,clade="all"))
    phenvec.list = lapply(seq_len(ncol(phenvec.table)), function(i) phenvec.table[,i])
    
    print("Calculating correlations")
    corMatList = lapply(phenvec.list, correlateWithBinaryPhenotype, RERmat=RERmat)
    
    #make enrich list/matrices to fill
    permPvals=data.frame(matrix(ncol=numperms, nrow=nrow(RERmat)))
    rownames(permPvals)=rownames(RERmat)
    permRhovals=data.frame(matrix(ncol=numperms, nrow=nrow(RERmat)))
    rownames(permRhovals)=rownames(RERmat)
    permStatvals=data.frame(matrix(ncol=numperms, nrow=nrow(RERmat)))
    rownames(permStatvals)=rownames(RERmat)
    
    for (i in 1:length(corMatList)){
      permPvals[,i] = corMatList[[i]]$P
      permRhovals[,i] = corMatList[[i]]$Rho
      permStatvals[,i] = sign(corMatList[[i]]$Rho)*-log10(corMatList[[i]]$P)
    }
    
  } else if (permmode=="ssm"){
    print("Running SSM permulation. sisters_list is required only for enrichments, otherwise sisters_list = NA is sufficient.")
    
    if (is.null(trees_list)){
      trees_list = trees$trees
    }
    
    RERmat = RERmat[match(names(trees_list), rownames(RERmat)),]
    
    print("Generating permulated trees")
    permulated.binphens = generatePermulatedBinPhenSSMBatched(trees_list,numperms,trees,root_sp,fg_vec,sisters_list,pathvec,permmode="ssm")
    
    # Get species membership of the trees
    df.list = lapply(trees_list,getSpeciesMembershipStats,masterTree=mastertree,foregrounds=fg_vec)
    df.converted = data.frame(matrix(unlist(df.list), nrow=length(df.list), byrow=T),stringsAsFactors=FALSE)
    attr = attributes(df.list[[1]])
    col_names = attr$names
    attr2 = attributes(df.list)
    row_names = attr2$names
    
    colnames(df.converted) = col_names
    rownames(df.converted) = row_names
    
    df.converted$num.fg = as.integer(df.converted$num.fg)
    df.converted$num.spec = as.integer(df.converted$num.spec)
    
    spec.members = df.converted$spec.members
    
    # Group gene trees based on the similarity of their species membership
    grouped.trees = groupTrees(spec.members)
    ind.unique.trees = grouped.trees$ind.unique.trees
    ind.unique.trees = unlist(ind.unique.trees)
    ind.tree.groups = grouped.trees$ind.tree.groups
    
    # For each unique tree, produce a permuted tree. We already have this function, but we need a list of trees to feed in.
    unique.trees = trees_list[ind.unique.trees]
    
    # precompute clade mapping for each unique tree
    unique.map.list = mapply(matchAllNodesClades,unique.trees,MoreArgs=list(treesObj=trees))
    
    # calculate paths for each permulation
    unique.permulated.binphens = permulated.binphens[ind.unique.trees]
    unique.permulated.paths = calculatePermulatedPaths_apply(unique.permulated.binphens,unique.map.list,trees)
    
    permulated.paths = vector("list", length = length(trees_list))
    for (j in 1:length(permulated.paths)){
      permulated.paths[[j]] = vector("list",length=numperms)
    }
    for (i in 1:length(unique.permulated.paths)){
      ind.unique.tree = ind.unique.trees[i]
      ind.tree.group = ind.tree.groups[[i]]
      unique.path = unique.permulated.paths[[i]]
      for (k in 1:length(ind.tree.group)){
        permulated.paths[[ind.tree.group[k]]] = unique.path
      }
    }
    attributes(permulated.paths)$names = row_names
    
    print("Calculating correlations")
    RERmat.list = lapply(seq_len(nrow(RERmat[])), function(i) RERmat[i,])
    corMatList = mapply(calculateCorPermuted,permulated.paths,RERmat.list)
    permPvals = extractCorResults(corMatList,numperms,mode="P")
    rownames(permPvals) = names(trees_list)
    permRhovals = extractCorResults(corMatList,numperms,mode="Rho")
    rownames(permRhovals) = names(trees_list)
    permStatvals = sign(permRhovals)*-log10(permPvals)
    rownames(permStatvals) = names(trees_list)
    
  } else if (permmode=="ssmLegacy"){
    print("Running SSM Legacy permulation")
    
    if (is.null(trees_list)){
      trees_list = trees$trees
    }
    
    RERmat = RERmat[match(names(trees_list), rownames(RERmat)),]
    
    print("Generating permulated trees")
    permulated.binphens = generatePermulatedBinPhenSSMBatched(trees_list,numperms,trees,root_sp,fg_vec,sisters_list,pathvec,permmode="ssmLegacy")
    
    # Get species membership of the trees
    df.list = lapply(trees_list,getSpeciesMembershipStats,masterTree=mastertree,foregrounds=fg_vec)
    df.converted = data.frame(matrix(unlist(df.list), nrow=length(df.list), byrow=T),stringsAsFactors=FALSE)
    attr = attributes(df.list[[1]])
    col_names = attr$names
    attr2 = attributes(df.list)
    row_names = attr2$names
    
    colnames(df.converted) = col_names
    rownames(df.converted) = row_names
    
    df.converted$num.fg = as.integer(df.converted$num.fg)
    df.converted$num.spec = as.integer(df.converted$num.spec)
    
    spec.members = df.converted$spec.members
    
    # Group gene trees based on the similarity of their species membership
    grouped.trees = groupTrees(spec.members)
    ind.unique.trees = grouped.trees$ind.unique.trees
    ind.unique.trees = unlist(ind.unique.trees)
    ind.tree.groups = grouped.trees$ind.tree.groups
    
    # For each unique tree, produce a permuted tree. We already have this function, but we need a list of trees to feed in.
    unique.trees = trees_list[ind.unique.trees]
    
    # precompute clade mapping for each unique tree
    unique.map.list = mapply(matchAllNodesClades,unique.trees,MoreArgs=list(treesObj=trees))
    
    # calculate paths for each permulation
    unique.permulated.binphens = permulated.binphens[ind.unique.trees]
    unique.permulated.paths = calculatePermulatedPaths_apply(unique.permulated.binphens,unique.map.list,trees)
    
    permulated.paths = vector("list", length = length(trees_list))
    for (j in 1:length(permulated.paths)){
      permulated.paths[[j]] = vector("list",length=numperms)
    }
    for (i in 1:length(unique.permulated.paths)){
      ind.unique.tree = ind.unique.trees[i]
      ind.tree.group = ind.tree.groups[[i]]
      unique.path = unique.permulated.paths[[i]]
      for (k in 1:length(ind.tree.group)){
        permulated.paths[[ind.tree.group[k]]] = unique.path
      }
    }
    attributes(permulated.paths)$names = row_names
    
    print("Calculating correlations")
    RERmat.list = lapply(seq_len(nrow(RERmat[])), function(i) RERmat[i,])
    corMatList = mapply(calculateCorPermuted,permulated.paths,RERmat.list)
    permPvals = extractCorResults(corMatList,numperms,mode="P")
    rownames(permPvals) = names(trees_list)
    permRhovals = extractCorResults(corMatList,numperms,mode="Rho")
    rownames(permRhovals) = names(trees_list)
    permStatvals = sign(permRhovals)*-log10(permPvals)
    rownames(permStatvals) = names(trees_list)
    
  } else {
    stop("Invalid binary permulation mode.")
  }
  
  if (calculateenrich){
    realFgtree = foreground2TreeClades(fg_vec, sisters_list, trees, plotTree=F)
    realpaths = tree2PathsClades(realFgtree, trees)
    realresults = getAllCor(RERmat, realpaths, method=method, min.pos=min.pos)
    realstat =sign(realresults$Rho)*-log10(realresults$P)
    names(realstat) = rownames(RERmat)
    realenrich = fastwilcoxGMTall(na.omit(realstat), annotlist, outputGeneVals=F)
    
    #sort real enrichments
    groups=length(realenrich)
    c=1
    while(c<=groups){
      current=realenrich[[c]]
      realenrich[[c]]=current[order(rownames(current)),]
      c=c+1
    }
    #make matrices to fill
    permenrichP=vector("list", length(realenrich))
    permenrichStat=vector("list", length(realenrich))
    c=1
    while(c<=length(realenrich)){
      newdf=data.frame(matrix(ncol=numperms, nrow=nrow(realenrich[[c]])))
      rownames(newdf)=rownames(realenrich[[c]])
      permenrichP[[c]]=newdf
      permenrichStat[[c]]=newdf
      c=c+1
    }
    
    counter=1;
    while (counter <= numperms){
      stat = permStatvals[,counter]
      names(stat) = rownames(RERmat)
      enrich=fastwilcoxGMTall(na.omit(stat), annotlist, outputGeneVals=F)
      #sort and store enrichment results
      groups=length(enrich)
      c=1
      while(c<=groups){
        current=enrich[[c]]
        enrich[[c]]=current[order(rownames(current)),]
        enrich[[c]]=enrich[[c]][match(rownames(permenrichP[[c]]), rownames(enrich[[c]])),]
        permenrichP[[c]][,counter]=enrich[[c]]$pval
        permenrichStat[[c]][,counter]=enrich[[c]]$stat
        c=c+1
      }
      counter = counter+1
    }
  }
  
  if(calculateenrich){
    data=vector("list", 5)
    data[[1]]=permPvals
    data[[2]]=permRhovals
    data[[3]]=permStatvals
    data[[4]]=permenrichP
    data[[5]]=permenrichStat
    names(data)=c("corP", "corRho", "corStat", "enrichP", "enrichStat")
  } else {
    data=vector("list", 3)
    data[[1]]=permPvals
    data[[2]]=permRhovals
    data[[3]]=permStatvals
    names(data)=c("corP", "corRho", "corStat")
  }
  data
}


#Load RERconverge output from murine rodent analysis
load("../../MiscData/RERconverge_output.logRTM_binary.OUmodel.RTMspeciesOnly.rds")


fgspec<-c("Pseudomys_novaehollandiae_ABTC08140", "Pseudomys_delicatulus_U1509", "Notomys_alexis_U1308", "Notomys_fuscus_M22830", "Notomys_mitchellii_M21518", "Zyzomys_pedunculatus_Z34925", "Bandicota_indica_ABTC119185", "Nesokia_indica_ABTC117074", "Hyomys_goliath_ABTC42697", "Pseudomys_shortridgei_Z25113", "Paruromys_dominator_JAE4870", "Eropeplus_canus_NMVZ21733")


newPerms = getPermsBinary(100, fgspec, NA, myTrees$masterTree$tip.label[1], myRER, myTrees, myTrees$masterTree, permmode = "cc")
legacyPerms = getPermsBinary(100, fgspec, NA, myTrees$masterTree$tip.label[1], myRER, myTrees, myTrees$masterTree, permmode = "ccLegacy")

numPerms = 100
trees = myTrees
permTrees<-list()
permFG_list<-list()
inFG<-c()
Sys.time()
pdf("Results/Emily/emilyNewPermulationPlot.pdf", onefile=TRUE, height=8.5, width=11)
for(i in 1:numPerms){
   
  phenotypeVector = rep(0, length(trees$masterTree$tip.label))
  names(phenotypeVector) = trees$masterTree$tip.label
  phenotypeVector[names(phenotypeVector) %in% fgspec] = 1
  
  permulationData = categoricalPermulations(trees, phenotypeVector, rm = "ER", rp = "auto", ntrees = 1)
    
    
    permulatedTrees = lapply(permulationData$trees, function(x) {
      tr = trees$masterTree
      tr$edge.length = c(x$tips, x$nodes)[tr$edge[,2]]
      tr$edge.length = tr$edge.length-1
      names(tr$edge.length) = NULL
      #tree2Paths(tr, treesObj, categorical = TRUE, useSpecies = names(phenvals))
      tr
    })
  permTree = permulatedTrees[[1]]
  #permTree<-getPermsBinary(1, fgspec, NA, myTrees$masterTree$tip.label[1], myRER, myTrees, myTrees$masterTree, permmode = "cc")
  permTrees[[i]]<-permTree
  fgEdges<-permTree$edge[which(permTree$edge.length==1),2]
  permFgs<-permTree$tip.label[fgEdges]
  permFG_list[[i]]<-permFgs
  inFG<-cbind(inFG, unlist(sapply(myTrees$masterTree$tip.label, function(x) if(x %in% permFgs){1} else{0})))
}
dev.off()
Sys.time()
inFG_sums<-rowSums(inFG)
inFG_props<-inFG_sums/numPerms
shortnames<-unlist(sapply(names(inFG_sums), function(x) paste(strsplit(x, "_")[[1]][1:2], collapse="_")))
names(inFG_sums)<-shortnames
names(inFG_props)<-shortnames
pdf("Results/Emily/newCCBarplot.pdf", onefile=TRUE, height=8.5, width=11)
barplot(height=inFG_sums, main=paste("Number of times each species appeared in the foreground out of", numPerms, "permulations\nccNew Permulations"), las=2, cex.names=0.5)
barplot(height=inFG_props, main=paste("Proportion of times each species appeared in the foreground out of", numPerms, "permulations\nccNew Permulations"), las=2, cex.names=0.5, ylim=c(0,0.8))
dev.off()

# ----
permTrees<-list()
permFG_list<-list()
inFG<-c()
sisters_list<-list("clade1"=c("Pseudomys_novaehollandiae_ABTC08140","Pseudomys_delicatulus_U1509"), "clade2"=c("Notomys_alexis_U1308","Notomys_fuscus_M22830"), "clade3"=c("clade2","Notomys_mitchellii_M21518"), "clade4"=c("Bandicota_indica_ABTC119185","Nesokia_indica_ABTC117074"), "clade5"=c("Paruromys_dominator_JAE4870", "Eropeplus_canus_NMVZ21733"))
fg_vec = fgspec
root_sp = trees$masterTree$tip.label[[1]]
Sys.time()
pdf("Results/Emily/emilyLegacyPermulationPlot.pdf", onefile=TRUE, height=8.5, width=11)
for(i in 1:numPerms){
  
  permulated.binphens = generatePermulatedBinPhen(trees$masterTree, 1, trees, root_sp, fg_vec, sisters_list, pathvec, permmode="cc")
  permTree = permulated.binphens[[1]][[1]]
  #permTree<-getPermsBinary(1, fgspec, NA, myTrees$masterTree$tip.label[1], myRER, myTrees, myTrees$masterTree, permmode = "cc")
  permTrees[[i]]<-permTree
  fgEdges<-permTree$edge[which(permTree$edge.length==1),2]
  permFgs<-permTree$tip.label[fgEdges]
  permFG_list[[i]]<-permFgs
  inFG<-cbind(inFG, unlist(sapply(myTrees$masterTree$tip.label, function(x) if(x %in% permFgs){1} else{0})))
}
dev.off()
Sys.time()
#inFGNew <- apply(inFG, 2, as.numeric)
#rownames(inFGNew) = rownames(inFG)
inFG_sums<-rowSums(inFG)
inFG_props<-inFG_sums/numPerms
shortnames<-unlist(sapply(names(inFG_sums), function(x) paste(strsplit(x, "_")[[1]][1:2], collapse="_")))
names(inFG_sums)<-shortnames
names(inFG_props)<-shortnames
pdf("Results/Emily/legacyCCBarPlot.pdf", onefile=TRUE, height=8.5, width=11)
barplot(height=inFG_sums, main=paste("Number of times each species appeared in the foreground out of", numPerms, "permulations\nccLegacy Permulations"), las=2, cex.names=0.5)
barplot(height=inFG_props, main=paste("Proportion of times each species appeared in the foreground out of", numPerms, "permulations\nccLegacy Permulations"), las=2, cex.names=0.5, ylim=c(0,0.8))
dev.off()



# ------------------------------------------------------------------
# --- Looking into differecen between liam and non-liam results  ----- 
# ------------------------------------------------------------------





nonLiamResults = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsvertivoreTreecombinedGeneResults.rds")
liamResults = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsvertivoreTreeLiamInferencecombinedGeneResults.rds")

which(!rownames(liamResults) %in% rownames(nonLiamResults))

liamHI = rownames(liamResults)[which(liamResults$`HI-p.adj` <0.05)]
nonliamHI = rownames(nonLiamResults)[which(nonLiamResults$`HI-p.adj` <0.05)]

length(which(liamHI %in% nonliamHI))


nonLiamResultsGO = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsvertivoreTreecombinedGOResults-KeggReactome.rds")
liamResultsGO = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsvertivoreTreeLiamInferencecombinedGOResults-KeggReactome.rds")



liamHIGO = rownames(liamResultsGO)[which(liamResultsGO$`HI-p.adj` <0.1)]
nonliamHIGO = rownames(nonLiamResultsGO)[which(nonLiamResultsGO$`HI-p.adj` <0.1)]


length(which(liamHIGO %in% nonliamHIGO))


names(liamResults) = paste0("liam-", names(liamResults))
names(liamResultsGO) = paste0("liam-", names(liamResultsGO))



which(liamResults$`HI-significant`)


# Get significance values for genes
significanceColumns = names(liamResults)[grep("significant", names(liamResults))]
geneSignificanceResults = liamResults[, names(liamResults) %in% significanceColumns]
colSums(geneSignificanceResults, na.rm = T)

significanceColumns = names(nonLiamResults)[grep("significant", names(nonLiamResults))]
geneSignificanceResultsNL = nonLiamResults[, names(nonLiamResults) %in% significanceColumns]
colSums(geneSignificanceResultsNL, na.rm = T)
#Remove the ch column from the non-liam results
geneSignificanceResultsNL = geneSignificanceResultsNL[,-1]


matchGeneSignificant = geneSignificanceResults == geneSignificanceResultsNL
matchGeneSignificant[!geneSignificanceResults & !geneSignificanceResultsNL] = NA
totalSignificant = colSums(!is.na(matchGeneSignificant))
sharedSignificant = colSums(matchGeneSignificant, na.rm = T)

sharedSignificant/totalSignificant 



significanceColumnsGO = names(liamResultsGO)[grep("significant", names(liamResultsGO))]
geneSignificanceResultsGO = liamResultsGO[, names(liamResultsGO) %in% significanceColumnsGO]
colSums(geneSignificanceResultsGO, na.rm = T)

significanceColumnsGO = names(nonLiamResultsGO)[grep("significant", names(nonLiamResultsGO))]
geneSignificanceResultsGONL = nonLiamResultsGO[, names(nonLiamResultsGO) %in% significanceColumnsGO]
colSums(geneSignificanceResultsGONL, na.rm = T)

geneSignificanceResultsGONL = geneSignificanceResultsGONL[,c(1,4,2,3,5,6)]

matchGeneSignificantGO = geneSignificanceResultsGO == geneSignificanceResultsGONL
matchGeneSignificantGO[!geneSignificanceResultsGO & !geneSignificanceResultsGONL] = NA
totalSignificantGO = colSums(!is.na(matchGeneSignificantGO))
sharedSignificantGO = colSums(matchGeneSignificantGO, na.rm = T)

sharedSignificantGO/totalSignificantGO 



# ------------------------------------------------------------------
# --- Making liam vs non-liam tree plot  ----- 
# ------------------------------------------------------------------

liamTree = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCategoricalTree.rds")
nonliamTree = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalTree.rds")

all.equal(liamTree$edge, nonliamTree$edge)
which(! liamTree$edge == nonliamTree$edge)



which(liamTree$edge.length != nonliamTree$edge.length)
length(which(liamTree$edge.length != nonliamTree$edge.length))

mainTrees = readRDS("Data/zoonomiaAllMammalsTrees.rds")
commonMainTrees = mainTrees
source("Src/Reu/ZoonomTreeNameToCommon.R")
nameColumn = "ZoonomiaTip"
spreadSheetLocation = "Data/MergedData.csv" 
commonMainTrees$masterTree = ZoonomTreeNameToCommon(commonMainTrees$masterTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)
commonLiamTree = ZoonomTreeNameToCommon(liamTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)
commonNonliamTree = ZoonomTreeNameToCommon(nonliamTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)

equalBranchMaster = commonMainTrees$masterTree
equalBranchMaster$edge.length = rep(1, length(equalBranchMaster$edge.length ))

edgelabelcolor = rep("black", length(liamTree$edge.length))
edgelabelcolor[which(liamTree$edge.length != nonliamTree$edge.length)] = "purple"

palette(c( "darkgreen", "darkblue","black", "red"))

pdf("results/CompareLiamTrees.pdf", height = length(liamTree$tip.label)/10, width = 20)     
par(mfrow = c(1,2))
plotTreeCategorical(commonNonliamTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = equalBranchMaster)
edgelabels(bg=NULL, frame = "none", col = edgelabelcolor, cex = 0.3)

plotTreeCategorical(commonLiamTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = equalBranchMaster)
edgelabels(bg=NULL, frame = "none", col = edgelabelcolor, cex = 0.3)
 


plotTreeCategorical(commonNonliamTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = commonMainTrees$masterTree)
edgelabels(bg=NULL, frame = "none", col = edgelabelcolor, cex = 0.3)

plotTreeCategorical(commonLiamTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = commonMainTrees$masterTree)
edgelabels(bg=NULL, frame = "none", col = edgelabelcolor, cex = 0.3)


dev.off()  



# ------------------------------------------------------------------
# ---  Make New SUpplementary File 3----- 
# ------------------------------------------------------------------

MGI = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsvertivoreTreeLiamInferencecombinedGOResults-MGI_Mammalian_Phenotype_Level_4.rds")
GO = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsvertivoreTreeLiamInferencecombinedGOResults-GO_Biological_Process_2023.rds")
DisGen = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsvertivoreTreeLiamInferencecombinedGOResults-DisGeNet.rds")
tissue = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsvertivoreTreeLiamInferencecombinedGOResults-tissue_specific.rds")
Enrichment = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsvertivoreTreeLiamInferencecombinedGOResults-EnrichmentHsSymbolsFile2.rds")


MGI$GeneSet = rep("MGI Mammalian Phenotype", nrow(MGI))
GO$GeneSet = rep("GO Biological Process", nrow(GO))
DisGen$GeneSet = rep("DisGeNET", nrow(DisGen))
tissue$GeneSet = rep("Tissue Specific", nrow(tissue))
Enrichment$GeneSet = rep("EnrichmentHsSymbols", nrow(Enrichment))


combinedGenesets = rbind(MGI, GO, DisGen, tissue, Enrichment)

combinedGenesets = combinedGenesets[,c(1,296, 2:295)]
write.csv(combinedGenesets, "Output/CategoricalInsVertivoreTreeLiamInference/CombinedGenesetOutput.csv")

# ------------------------------------------------------------------
# ---  Make newick files for trees ----- 
# ------------------------------------------------------------------

phenotypeTree3Diet = readRDS("Output/CategoricalInsVertivoreTreeCarnivoreLiamInference/CategoricalInsVertivoreTreeCarnivoreLiamInferenceCategoricalCommonTree.rds")
phenotypeTree4Diet = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCategoricalCommonTree.rds")
mainTrees = readRDS("Data/zoonomiaAllMammalsTrees.rds")
masterTree = mainTrees$masterTree


write.tree(phenotypeTree3Diet,"Results/PhenTree3.txt")
write.tree(phenotypeTree4Diet,"Results/PhenTree4.txt")
write.tree(masterTree,"Results/masterTree.txt")



# ------------------------------------------------------------------
# ---  Make Liam Radial Tree ----- 
# ------------------------------------------------------------------

liamTree = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCategoricalTree.rds")
nonliamTree = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalTree.rds")


match(liamTree$tip.label, nonliamTree$tip.label)

liamNonliamTipConversionIndex = match(nonliamTree$tip.label, liamTree$tip.label)


{
collapsedClades = data.frame()
collapsedClades[1,] = NA

collapsedClades$Platypus = MRCA(commonCategoricalTree, 195)
collapsedClades$Opossums = MRCA(commonCategoricalTree, c(194,192,193))
collapsedClades$Koala = MRCA(commonCategoricalTree, c(190,191))
collapsedClades$Kangaroos = MRCA(commonCategoricalTree, c(186,187,189,188))
collapsedClades$Anteaters = MRCA(commonCategoricalTree, c(183,182))
collapsedClades$Sloths = MRCA(commonCategoricalTree, c(181,180))
collapsedClades$Elephant= MRCA(commonCategoricalTree, c(178,177,176))
collapsedClades$Aardvark = MRCA(commonCategoricalTree, c(175,174,171,172,173))
collapsedClades$Strepsirrhini = MRCA(commonCategoricalTree, c(28,29,27,26,24,25,23,21,22,20,19))
collapsedClades$Atelidae = MRCA(commonCategoricalTree, c(6,5,4,3,2,1))
collapsedClades$Chimpanze = MRCA(commonCategoricalTree, c(18,17,16,15,14,12,13,10,9,8,7,11))
collapsedClades$Hares = MRCA(commonCategoricalTree, c(71,70))
collapsedClades$Squirrels = MRCA(commonCategoricalTree, c(69,67,68,61,66,65,63,62,64))
collapsedClades$Capybara = MRCA(commonCategoricalTree, c(33,32,31,30))
collapsedClades$Beaver = MRCA(commonCategoricalTree, c(36,35,34))
collapsedClades$Jerboa = MRCA(commonCategoricalTree, c(59,58,57))
collapsedClades$Deomyinae = MRCA(commonCategoricalTree, c(41,39,40,38,37,44,43,42))
collapsedClades$Vole = MRCA(commonCategoricalTree, c(48,47,46,45))
collapsedClades$Neotominae = MRCA(commonCategoricalTree, c(54,49,50,52,51,53))
collapsedClades$`African Hedgehogs` = MRCA(commonCategoricalTree, c(168,169))
collapsedClades$`Talpa europaea` = MRCA(commonCategoricalTree, c(167,165,164,166))
collapsedClades$`Flying Fox` = MRCA(commonCategoricalTree, c(163,161,162))
collapsedClades$Rhinolophidae = MRCA(commonCategoricalTree, c(155,156,158,157,159,160))
collapsedClades$`Big Brown Bat` = MRCA(commonCategoricalTree, c(143,142,139,140,141))
collapsedClades$Phyllostomidae = MRCA(commonCategoricalTree, c(144,145,147,146))
collapsedClades$Noctilio = MRCA(commonCategoricalTree, c(154))
collapsedClades$Horse = MRCA(commonCategoricalTree, c(100,99,98))
collapsedClades$Pig = MRCA(commonCategoricalTree, c(95,96))
collapsedClades$`bos bison` = MRCA(commonCategoricalTree, c(87,89,88))
collapsedClades$`Humpback Whale` = MRCA(commonCategoricalTree, c(75,74,72,73))
collapsedClades$`Dolphins` = MRCA(commonCategoricalTree, c(82,81,80,77,76))
collapsedClades$`Pangolin` = MRCA(commonCategoricalTree, c(138,137))
collapsedClades$`Lion` = MRCA(commonCategoricalTree, c(131,130,129))
collapsedClades$`Meerkat` = MRCA(commonCategoricalTree, c(135,134))
collapsedClades$`Dog` = MRCA(commonCategoricalTree, c(101,102))
collapsedClades$`Brown Bear` = MRCA(commonCategoricalTree, c(125,128,127))
collapsedClades$`Odobenus rosmarus` = MRCA(commonCategoricalTree, c(120,119,118,117))
collapsedClades$`Phocidae` = MRCA(commonCategoricalTree, c(124,123,121,122))
collapsedClades$`Procyon lotor` = MRCA(commonCategoricalTree, c(114,113,112))
collapsedClades$`Lontra provocax` = MRCA(commonCategoricalTree, c(108,107,106,105,104))
collapsedClades$`Tasmanian Devil` = MRCA(commonCategoricalTree, c(184,185))
}




{
  collapsedClades = data.frame()
  collapsedClades[1,] = NA
  
  collapsedClades$Platypus = MRCA(commonCategoricalTree, c(1))
  collapsedClades$Opossums = MRCA(commonCategoricalTree, c(3,4,5))
  collapsedClades$Koala = MRCA(commonCategoricalTree, c(8,9))
  collapsedClades$Kangaroos = MRCA(commonCategoricalTree, c(10,11,12,13))
  collapsedClades$Anteaters = MRCA(commonCategoricalTree, c(23,24))
  collapsedClades$Sloths = MRCA(commonCategoricalTree, c(25,26))
  collapsedClades$Elephant= MRCA(commonCategoricalTree, c(19,20,21))
  collapsedClades$Aardvark = MRCA(commonCategoricalTree, c(14,15,16,17,18))
  collapsedClades$Strepsirrhini = MRCA(commonCategoricalTree, c(27,28,29,30,31,32,33,34,35,36,37))
  collapsedClades$Atelidae = MRCA(commonCategoricalTree, c(38,39,40,41,42,43))
  collapsedClades$Chimpanze = MRCA(commonCategoricalTree, c(44:55))
  collapsedClades$Hares = MRCA(commonCategoricalTree, c(56,57))
  collapsedClades$Squirrels = MRCA(commonCategoricalTree, c(58:66))
  collapsedClades$Capybara = MRCA(commonCategoricalTree, c(67:70))
  collapsedClades$Beaver = MRCA(commonCategoricalTree, c(72:74))
  collapsedClades$Jerboa = MRCA(commonCategoricalTree, c(75:77))
  collapsedClades$Deomyinae = MRCA(commonCategoricalTree, c(90:97))
  collapsedClades$Vole = MRCA(commonCategoricalTree, c(86:89))
  collapsedClades$Neotominae = MRCA(commonCategoricalTree, c(80:85))
  collapsedClades$`African Hedgehogs` = MRCA(commonCategoricalTree, c(99,100))
  collapsedClades$`Talpa europaea` = MRCA(commonCategoricalTree, c(101:104))
  collapsedClades$`Flying Fox` = MRCA(commonCategoricalTree, c(105:107))
  collapsedClades$Rhinolophidae = MRCA(commonCategoricalTree, c(108:113))
  collapsedClades$`Big Brown Bat` = MRCA(commonCategoricalTree, c(125:129))
  collapsedClades$Phyllostomidae = MRCA(commonCategoricalTree, c(121:124))
  collapsedClades$Noctilio = MRCA(commonCategoricalTree, c(114))
  collapsedClades$Horse = MRCA(commonCategoricalTree, c(168:170))
  collapsedClades$Pig = MRCA(commonCategoricalTree, c(172:173))
  collapsedClades$`bos bison` = MRCA(commonCategoricalTree, c(194:196))
  collapsedClades$`Humpback Whale` = MRCA(commonCategoricalTree, c(175:178))
  collapsedClades$`Dolphins` = MRCA(commonCategoricalTree, c(182:186))
  collapsedClades$`Pangolin` = MRCA(commonCategoricalTree, c(130:131))
  collapsedClades$`Lion` = MRCA(commonCategoricalTree, c(137:139))
  collapsedClades$`Meerkat` = MRCA(commonCategoricalTree, c(135:136))
  collapsedClades$`Dog` = MRCA(commonCategoricalTree, c(140:141))
  collapsedClades$`Brown Bear` = MRCA(commonCategoricalTree, c(142:144))
  collapsedClades$`Odobenus rosmarus` = MRCA(commonCategoricalTree, c(146:149))
  collapsedClades$`Phocidae` = MRCA(commonCategoricalTree, c(150:153))
  collapsedClades$`Procyon lotor` = MRCA(commonCategoricalTree, c(156:158))
  collapsedClades$`Lontra provocax` = MRCA(commonCategoricalTree, c(162:166))
  collapsedClades$`Tasmanian Devil` = MRCA(commonCategoricalTree, c(6:7))
  
  #collapsedClades$Cats = MRCA(commonCategoricalTree, c("Jaguar", "Lion", "Cheetah"))
}


#collapsedClades$Cats = MRCA(commonCategoricalTree, c("Jaguar", "Lion", "Cheetah"))

# ------------------------------------------------------------------
# ---  check seize ----- 
# ------------------------------------------------------------------

phenoTree = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCategoricalTree.rds")

table(phenoTree$edge.length)


data = read.csv("Data/mergedData.csv")
table(data$DerekDietClassification90InsVertivoreSorting)
grep("Piscivore")

# ------------------------------------------------------------------
# ---  Examine cortisol results ----- 
# ------------------------------------------------------------------

which(rownames(GoMetaCombined) == "REACTOME_METABOLISM_OF_STEROID_HORMONES")

which(gmt$geneset.names == "REACTOME_METABOLISM_OF_STEROID_HORMONES")

steroidGenes = gmt$genesets[348][[1]]

steriodCHSignifiance = geneMetaCombined$`CH-significant-Liam`[match(steroidGenes, rownames(geneMetaCombined))]
steriodHISignifiance = geneMetaCombined$`HI-significant-Liam`[match(steroidGenes, rownames(geneMetaCombined))]
steriodHVSignifiance = geneMetaCombined$`HV-significant-Liam`[match(steroidGenes, rownames(geneMetaCombined))]



steriodSignificance = data.frame(steroidGenes, steriodCHSignifiance, steriodHISignifiance, steriodHVSignifiance)

genesAllwaysUnsignificant = which(apply(steriodSignificance[,2:4], 1, function(x) all(x == FALSE)))
genesNA = which(apply(steriodSignificance[,2:4], 1, function(x) all(is.na(x))))


steriodSignificant = steriodSignificance[-c(genesAllwaysUnsignificant, genesNA),]


CHOnlyGenes = steriodSignificant[c(which(apply(steriodSignificant[,3:4], 1, function(x) all(x == FALSE))),which(apply(steriodSignificant[,3:4], 1, function(x) all(x == TRUE)))),]
steriodDifferences = steriodSignificant[-c(which(apply(steriodSignificant[,3:4], 1, function(x) all(x == FALSE))),which(apply(steriodSignificant[,3:4], 1, function(x) all(x == TRUE)))),]




# ------------------------------------------------------------------
# ---  Make new venn diagram ----- 
# ------------------------------------------------------------------

vennInputDataframe = vennGoSignificanceResults

makeVennPlot = function(vennInputDataframe, mainTitle, plot = T){
  # Build logical vectors for each set
  set1 <- vennInputDataframe[1] == TRUE
  set2 <- vennInputDataframe[2] == TRUE
  set3 <- vennInputDataframe[3] == TRUE
  
  # Create the Venn counts for each region
  vennCounts = c(
    "Column One" = sum(set1 & !set2 & !set3, na.rm = T),
    "Column Two" = sum(!set1 & set2 & !set3, na.rm = T),
    "Column Three" = sum(!set1 & !set2 & set3, na.rm = T),
    "Column One&Column Two" = sum(set1 & set2 & !set3, na.rm = T),
    "Column One&Column Three" = sum(set1 & !set2 & set3, na.rm = T),
    "Column Two&Column Three" = sum(!set1 & set2 & set3, na.rm = T),
    "Column One&Column Two&Column Three" = sum(set1 & set2 & set3, na.rm = T)
  )
  
  comparisonPrefixes = gsub("-significant", "", names(vennInputDataframe))
  comparisonPrefixes = gsub(commonBackground, "", comparisonPrefixes)
  comparisonNames = sapply(comparisonPrefixes, replacePrefixWithName)
  names(vennCounts)=c(comparisonNames[1], comparisonNames[2], comparisonNames[3], paste(comparisonNames[1], comparisonNames[2], sep="&"),paste(comparisonNames[1], comparisonNames[3], sep="&"),paste(comparisonNames[2], comparisonNames[3], sep="&"), paste(comparisonNames[1], comparisonNames[2], comparisonNames[3], sep="&"))
  
  # - make venn diagram labels, with the combination section on two lines and the solo sections on one 
  totalVennValues = sum(vennCounts)
  vennLabels = paste0(
    vennCounts, "\n (", round(vennCounts / totalVennValues * 100, 1), "%)"
  )
  for(i in 1:3){
    vennLabels[i] = paste0(
      vennCounts[i], " (", round(vennCounts[i] / totalVennValues * 100, 1), "%)"
    )
  }
  
  
  vennCountsCustom = vennCounts
  
  vennCountsCustom[5] = 80
  vennCountsCustom[1] = 40 
  vennCountsCustom[2] = 15
  
  vennLabelsCustom = vennLabels
  vennLabelsCustom[4] = 3
  
  # Create Euler diagram
  fit = euler(vennCountsCustom)
  outPlot = plot(fit,
                 fills = list(fill = vennColorset, alpha = 0.5),
                 labels = list(font = 4),
                 quantities = list(labels = vennLabelsCustom, font = 3),
                 main = mainTitle, theme)
  if(plot){print(outPlot)}
  return(outPlot)
}

GoVennCustom = outPlot

png(width = 1120, height = 560, file = "Output/CategoricalInsvertivoreTreeLiamInference/VennDiagramFigure.png")
combinedPlot = grid.arrange(geneVenn, GoVennCustom, nrow = 1)
dev.off()

# ------------------------------------------------------------------
# --- Add cytoscape index column ----- 
# ------------------------------------------------------------------

cytoscapeNodes1 = read.csv("Output/CategoricalInsVertivoreTreeLiamInference/Herbivore-Insectivore/Cytoscape/CarnivoreConserveddefaultnode.csv")
cytoscapeNodes2 = read.csv("Output/CategoricalInsVertivoreTreeLiamInference/Herbivore-Insectivore/Cytoscape/Bidirectionalnodes.csv")
cytoscapeNodes3 = read.csv("Output/CategoricalInsVertivoreTreeLiamInference/Herbivore-Insectivore/Cytoscape/Herbivorenodes.csv")
cytoscapeNodes4 = read.csv("Output/CategoricalInsVertivoreTreeLiamInference/Herbivore-Insectivore/Cytoscape/BidectionalNewNodes.csv")
cytoscapeNodes5 = read.csv("Output/CategoricalInsVertivoreTreeLiamInference/Insectivore-Vertivore/Cytoscape/IVNodes.csv")
cytoscapeNodes6 = read.csv("Output/CategoricalInsVertivoreTreeLiamInference/Herbivore-Insectivore/Cytoscape/NewPredatorNodes.csv")
cytoscapeNodes7 = read.csv("Output/CategoricalInsVertivoreTreeLiamInference/Herbivore-Insectivore/Cytoscape/NewBidirectionalNodes.csv")


cytoscapeNodes = cytoscapeNodes7

GOOutput = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsvertivoreTreeLiamInferencecombinedGOResults-KeggReactome.rds")

rownames(GOOutput)

cytoscapeNodes$shared.name

test = match(cytoscapeNodes$shared.name, rownames(GOOutput))

which(rownames(GOOutput) == "KEGG_BETA_ALANINE_METABOLISM")
which(rownames(GOOutput) == "REACTOME_DNA_REPAIR")



rownames(GOOutput)[436]

nodeIdexes = data.frame(cytoscapeNodes$shared.name, (match(cytoscapeNodes$shared.name, rownames(GOOutput)))+1)

write.csv(nodeIdexes, "Output/CategoricalInsVertivoreTreeLiamInference/Herbivore-Insectivore/Cytoscape/fixedBidirectionalIndex.csv")


length(which(GOOutput$`IV-significant`))

length(which(geneMetaCombined$`IV-significant-Liam`))

repairVsH = c("PID_FANCONI_PATHWAY","REACTOME_DISEASES_OF_DNA_REPAIR","REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR","REACTOME_DNA_REPAIR","REACTOME_FANCONI_ANEMIA_PATHWAY","REACTOME_HDR_THROUGH_HOMOLOGOUS_RECOMBINATION_HRR","REACTOME_HDR_THROUGH_SINGLE_STRAND_ANNEALING_SSA","REACTOME_HOMOLOGOUS_DNA_PAIRING_AND_STRAND_EXCHANGE","REACTOME_HOMOLOGY_DIRECTED_REPAIR","REACTOME_RESOLUTION_OF_D_LOOP_STRUCTURES","REACTOME_RESOLUTION_OF_D_LOOP_STRUCTURES_THROUGH_SYNTHESIS_DEPENDENT_STRAND_ANNEALING_SDSA")

repairIV = c("KEGG_HOMOLOGOUS_RECOMBINATION","REACTOME_DISEASES_OF_DNA_REPAIR","REACTOME_DISEASES_OF_MISMATCH_REPAIR_MMR","REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR","REACTOME_DNA_REPAIR","REACTOME_FANCONI_ANEMIA_PATHWAY","REACTOME_HDR_THROUGH_HOMOLOGOUS_RECOMBINATION_HRR","REACTOME_HDR_THROUGH_SINGLE_STRAND_ANNEALING_SSA","REACTOME_HOMOLOGOUS_DNA_PAIRING_AND_STRAND_EXCHANGE","REACTOME_HOMOLOGY_DIRECTED_REPAIR","REACTOME_RESOLUTION_OF_D_LOOP_STRUCTURES","REACTOME_RESOLUTION_OF_D_LOOP_STRUCTURES_THROUGH_SYNTHESIS_DEPENDENT_STRAND_ANNEALING_SDSA")

repairVsH %in% repairIV
repairIV %in% repairVsH


RnavVsH = c(
  "REACTOME_TRANSPORT_OF_THE_SLBP_DEPENDANT_MATURE_MRNA",
  "REACTOME_RESPONSE_OF_MTB_TO_PHAGOCYTOSIS",
  "REACTOME_LEISHMANIA_INFECTION",
  "REACTOME_NUCLEAR_ENVELOPE_NE_REASSEMBLY",
  "REACTOME_MRNA_SPLICING_MINOR_PATHWAY",
  "REACTOME_FORMATION_OF_THE_EARLY_ELONGATION_COMPLEX",
  "REACTOME_NEUROTOXICITY_OF_CLOSTRIDIUM_TOXINS",
  "REACTOME_HCMV_INFECTION",
  "REACTOME_FORMATION_OF_TC_NER_PRE_INCISION_COMPLEX",
  "REACTOME_FCERI_MEDIATED_NF_KB_ACTIVATION",
  "REACTOME_HIV_TRANSCRIPTION_ELONGATION",
  "REACTOME_METABOLISM_OF_RNA",
  "REACTOME_PROCESSING_OF_CAPPED_INTRON_CONTAINING_PRE_MRNA",
  "REACTOME_HIV_LIFE_CYCLE",
  "REACTOME_INFLUENZA_INFECTION",
  "REACTOME_HCMV_LATE_EVENTS",
  "REACTOME_HOST_INTERACTIONS_OF_HIV_FACTORS",
  "REACTOME_SNRNP_ASSEMBLY",
  "REACTOME_FORMATION_OF_RNA_POL_II_ELONGATION_COMPLEX",
  "REACTOME_TRANSPORT_OF_MATURE_MRNAS_DERIVED_FROM_INTRONLESS_TRANSCRIPTS",
  "REACTOME_SUPPRESSION_OF_PHAGOSOMAL_MATURATION",
  "REACTOME_RNA_POLYMERASE_II_PRE_TRANSCRIPTION_EVENTS",
  "REACTOME_TOXICITY_OF_BOTULINUM_TOXIN_TYPE_D_BOTD",
  "REACTOME_MRNA_SPLICING",
  "REACTOME_INFECTION_WITH_MYCOBACTERIUM_TUBERCULOSIS",
  "REACTOME_SUMOYLATION_OF_CHROMATIN_ORGANIZATION_PROTEINS",
  "REACTOME_HCMV_EARLY_EVENTS",
  "REACTOME_PREVENTION_OF_PHAGOSOMAL_LYSOSOMAL_FUSION",
  "REACTOME_HIV_INFECTION",
  "REACTOME_RNA_POLYMERASE_II_TRANSCRIPTION_TERMINATION",
  "REACTOME_SLBP_DEPENDENT_PROCESSING_OF_REPLICATION_DEPENDENT_HISTONE_PRE_MRNAS",
  "REACTOME_PROCESSING_OF_CAPPED_INTRONLESS_PRE_MRNA",
  "REACTOME_INFECTIOUS_DISEASE",
  "REACTOME_GLYCOLYSIS",
  "REACTOME_PROCESSING_OF_INTRONLESS_PRE_MRNAS",
  "KEGG_SPLICEOSOME",
  "REACTOME_ABORTIVE_ELONGATION_OF_HIV_1_TRANSCRIPT_IN_THE_ABSENCE_OF_TAT",
  "REACTOME_TRANSPORT_OF_MATURE_TRANSCRIPT_TO_CYTOPLASM"
)

RnavIV = c(
  "REACTOME_FORMATION_OF_RNA_POL_II_ELONGATION_COMPLEX",
  "REACTOME_HIV_INFECTION",
  "REACTOME_INFECTIOUS_DISEASE",
  "REACTOME_POTENTIAL_THERAPEUTICS_FOR_SARS",
  "REACTOME_RNA_POLYMERASE_II_PRE_TRANSCRIPTION_EVENTS",
  "REACTOME_RNA_POLYMERASE_II_TRANSCRIPTION",
  "REACTOME_SARS_COV_INFECTIONS",
  "REACTOME_TRANSCRIPTION_OF_THE_HIV_GENOME"
)


RnavIV[!RnavIV %in% RnavVsH]

# ------------------------------------------------------------------
# --- Update Enrichments ----- 
# ------------------------------------------------------------------

source("src/reu/RERConvergeFunctions.R")


getStat = function(res){
  stat=sign(res$Rho)*(-log10(res$P))
  names(stat)=rownames(res)
  #deal with duplicated genes
  genenames=sub("\\..*", "",names(stat))
  multname=names(which(table(genenames)>1))
  for(n in multname){
    ii=which(genenames==n)
    iimax=which(max(stat[ii])==max(abs(stat[ii])))
    stat[ii[-iimax]]=NA
  }
  sum(is.na(stat))
  stat=stat[!is.na(stat)]
  
  stat
}


fastwilcoxGMTall = function (vals, annotList, alternative = "two.sided", ...) 
{
  reslist = list()
  for (n in names(annotList)) {
    reslist[[n]] = fastwilcoxGMT(vals, annotList[[n]], alternative = alternative, 
                                 ...)
    message(paste0(nrow(reslist[[n]]), " results for annotation set ", 
                   n))
  }
  reslist
}


vals = rerStats
gmt = gmtAnnotations

fastwilcoxGMT=function(vals, gmt, simple=T, use.all=F, num.g=10,genes=NULL, outputGeneVals=F, order=F,
                       alternative = "two.sided"){
  vals=vals[!is.na(vals)]
  if(is.null(genes)){
    genes=unique(unlist(gmt$genesets))
  }
  out=matrix(nrow=length(gmt$genesets), ncol=5)
  rownames(out)=gmt$geneset.names
  colnames(out)=c("stat", "pval", "p.adj","num.genes", "gene.vals")
  out=as.data.frame(out)
  genes=intersect(genes, names(vals))
  
  valsr=rank(vals[genes])
  numg=length(vals)+1
  valsallr=rank(vals)
  for( i in 1:nrow(out)){
    
    curgenes=intersect(genes,gmt$genesets[[i]])
    
    bkgenes=setdiff(genes, curgenes)
    
    if (length(bkgenes)==0 || use.all){
      bkgenes=setdiff(names(vals), curgenes)
    }
    if(length(curgenes)>=num.g & length(bkgenes)>2){
      if(!simple){
        # change alternative = "greater" for the one-sided test
        res=wilcox.test(x = vals[curgenes], y=vals[bkgenes], exact=F, alternative = alternative)
        
        out[i, 1:2]=c(res$statistic/(as.numeric(length(bkgenes))*as.numeric(length(curgenes))), res$p.value)
      }
      else{
        # add an alternative parameter (can be "greater" or "two.sided")
        out[i, 1:2]=simpleAUCgenesRanks(valsr[curgenes],valsr[bkgenes], alt = alternative)
        
      }
      
      out[i,"num.genes"]=length(curgenes)
      if(outputGeneVals){
        if (out[i,1]>0.5){
          oo=order(vals[curgenes], decreasing = T)
          granks=numg-valsallr[curgenes]
        }
        else{
          oo=order(vals[curgenes], decreasing = F)
          granks=valsallr[curgenes]
        }
        
        
        nn=paste(curgenes[oo],round((granks[curgenes])[oo],2),sep=':' )
        out[i,"gene.vals"]=paste(nn, collapse = ", ")
      }
    }
    
  }
  # hist(out[,2])
  out[,1]=out[,1]-0.5
  out[, "p.adj"]=p.adjust(out[,2], method="BH")
  
  out=out[!is.na(out[,2]),]
  if(order){
    out=out[order(-abs(out[,1])),]
  }
  out
}


pos = valsr[curgenes]
neg = valsr[bkgenes]

simpleAUCgenesRanks=function(pos, neg, alt = "two.sided"){
  
  posn=length(pos)
  negn=length(neg)
  posn=as.numeric(posn)
  negn=as.numeric(negn)
  stat=sum(pos)-posn*(posn+1)/2 #Average pos
  auc=stat/(posn*negn)
  mu=posn*negn/2
  sd=sqrt((posn*negn*(posn+negn+1))/12)
  
  if(alt == "two.sided") {
    stattest=apply(cbind(stat, posn*negn-stat),1,max)
    pp=(2*pnorm(stattest, mu, sd, lower.tail = F))
  }
  
  else if(alt == "greater"){
    pp=(pnorm(stat,mu,sd,lower.tail=FALSE)) 
  }
  return(c(auc,pp))
}


# ------------------------------------------------------------------
# --- Make liam-non-liam comparison plots ----- 
# ------------------------------------------------------------------


# -- argument setup  -- 
significanceCutoff = 0.05
prefix = "CategoricalInsvertivoreTreeLiamInference"
pairwiseSets = c("Herbivore-Insectivore", "Herbivore-Vertivore", "Carnivore-Herbivore", "Herbivore-Omnivore", "Insectivore-Vertivore", "Omnivore-Vertivore", "Invertivore-Omnivore")
geneSet = "KeggReactome"
vennDiagramSet = c("Herbivore-Invertivore", "Herbivore-Vertivore", "Carnivore-Herbivore")  
vennColorset = c("darkblue", "red", "orange")
usingGo = !is.null(geneSet)
saveCombinedData = T
saveCombinedData = F
bothAxis = T
saveData = T
saveData = F

args = c("r=CategoricalInsvertivoreTree")
args = c("r=CategoricalInsvertivoreTreeLiamInference")


nonliamResults = combinedResults
nonliamGoResults = GoCombinedResults

liamResults = combinedResults
liamGoResults = GoCombinedResults

# -- Read Data 
combinedGeneDataFilename = paste0(outputFolderName, filePrefix, "combinedGeneResults.rds")
combinedResults = readRDS(combinedGeneDataFilename)
significanceColumns = names(combinedResults)[grep("significant", names(combinedResults))]
geneSignificanceResults = combinedResults[, names(combinedResults) %in% significanceColumns]


combinedGODataFilename = paste0(outputFolderName, filePrefix, "combinedGOResults-", geneSet, ".rds")
GoCombinedResults = readRDS(combinedGODataFilename)
GoSignificanceColumns = names(GoCombinedResults)[grep("significant", names(GoCombinedResults))]
GoSignificanceResults = GoCombinedResults[, names(GoCombinedResults) %in% GoSignificanceColumns]







# -- make resources to prefix-phenotype conversion 
prefixSet = NULL
prefixList = NULL
for(i in 1:length(pairwiseSets)){
  currentSet = pairwiseSets[i]
  correlationSubsetName = gsub("-", " - ", currentSet)
  correlationPrefix = paste(substr(strsplit(currentSet, split = "-")[[1]],1,1), collapse = '')
  prefixEntry = correlationPrefix; names(prefixEntry) = currentSet; prefixSet = append(prefixSet, prefixEntry)
  for(j in 1:2){
    prefixSingle = strsplit(correlationPrefix, split = "")[[1]][j]; names(prefixSingle) = strsplit(currentSet, split="-")[[1]][j]; prefixList = append(prefixList, prefixSingle)
  }
}


prefixList = unlist(prefixList); prefixList = prefixList[!duplicated(prefixList)]

# Debug code line for personal use 
names(prefixList)[which(names(prefixList) == "Insectivore")] = "Invertivore"

replacePrefixWithName = function(x) {
  inversePrefixList = setNames(names(prefixList), prefixList); parts = unlist(strsplit(x, "-")); fullNames = inversePrefixList[parts]; paste(fullNames, collapse = "-")
}

addDashes = function(vector) {
  sapply(vector, function(s) paste(strsplit(s, "")[[1]], collapse = "-"))
}


#-------------------------------------------------------------------



#gene
liamResults = liamGeneResults
nonliamResults = nonliamGeneResults

editedLiamResults = liamResults
colnames(editedLiamResults) = paste0(colnames(liamResults), "-Liam")
metaCombinedResults = cbind(nonliamResults, editedLiamResults)

combinedResults = metaCombinedResults
geneMetaCombined = metaCombinedResults


#GO
#nonliamGeneResults = nonliamResults
#liamGeneResults = liamResults

liamResults = liamGoResults
nonliamResults = nonliamGoResults

editedLiamResults = liamResults
colnames(editedLiamResults) = paste0(colnames(liamResults), "-Liam")
metaCombinedResults = cbind(nonliamResults, editedLiamResults)

names(combinedResults) = gsub("-stat", corrleationColumnType, names(combinedResults))

combinedResults = metaCombinedResults
GoMetaCombined = metaCombinedResults



#Prefix
prefixList = append(prefixList, c("L", "i", "a", "m", "-", ""))
names(prefixList) = append(names(prefixList)[1:5], c("Liam", "", "", "", "", ""))

bothAxis = F

corrleationColumnType = "-p.adj"
#-------------------------------------------------------------------

{
  
  
  
  grep(corrleationColumnType, names(combinedResults))
  rhoValues = combinedResults[,grep(corrleationColumnType, names(combinedResults))]
  
  rhoComparisions = names(rhoValues)
  
  #rhoComparisions = rhoComparisions[-c(5,6,7,9,13,14,15)]
  
  rhoPhenotypes = strsplit(gsub(corrleationColumnType, "", rhoComparisions), split = "")
  commonBackground = Reduce(intersect, rhoPhenotypes)
  if(length(commonBackground) == 1){
    for(i in 1:length(rhoPhenotypes)){ #invert tho if background in second postion so rho has consistent meaning relative to background
      if(rhoPhenotypes[[i]][1] != commonBackground){
        cat("Inverting rho of ", rhoPhenotypes[[i]] , "becuase background is in first position.")
        rhoValues[i] = -1*rhoValues[i]
      }
    }
  }
  
  densityScaleSet = NULL
  
  #generate the plots slim-ly to get the desired density scale 
  for(i in 1:length(rhoValues)){
    xName = names(rhoValues)[i]
    if(i+1 <= length(rhoValues)){
      for(j in (i+1):length(rhoValues)){
        yName = names(rhoValues)[j]
        
        rhoCorrellPlot = ggplot(rhoValues, aes(x = .data[[xName]], y = .data[[yName]])) + 
          geom_point() + geom_pointdensity() + scale_color_viridis()
        
        
        denstiyScaleValue = ggplot_build(rhoCorrellPlot)$plot$scales$scales[[1]]$get_limits()[2]
        densityScaleSet = append(densityScaleSet, denstiyScaleValue)
        rm(rhoCorrellPlot)
      }
    }
  }
  densityScale = c(1, max(densityScaleSet))
  
  
  rhoPlotSet = list()
  netIndex= 0
  for(i in 1:length(rhoValues)){
    xName = names(rhoValues)[i]
    if(i <= length(rhoValues)){
      if(bothAxis){jStart = 1}else{jStart = i+1}
      for(j in (jStart):length(rhoValues)){
        yName = names(rhoValues)[j]
        yLabel =  paste0(replacePrefixWithName(addDashes(gsub("-", "", gsub(corrleationColumnType, "", yName)))), " Dunn Z Statistic")
        xLabel =  paste0(replacePrefixWithName(addDashes(gsub("-", "", gsub(corrleationColumnType, "", xName)))), " Dunn Z Statistic")
        
        rhoCorrellPlot = ggplot(rhoValues, aes(x = .data[[xName]], y = .data[[yName]])) + 
          geom_point() + geom_pointdensity() + scale_color_viridis(name = "Density of genes", limits = densityScale) + 
          stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),formula = y ~ x,parse = TRUE, size = 6) +
          theme_classic()+
          xlab(xLabel) + ylab(yLabel)+
          theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))
        
        netIndex = netIndex +1
        rhoPlotSet[[netIndex]] = rhoCorrellPlot
        names(rhoPlotSet)[netIndex] = paste(xName, yName, sep="-")
        rm(rhoCorrellPlot)
      }
    }
  }
  # add a null plot by using two using permulations as a comparision 
  
  #pdf()
  #print(rhoPlotSet)
  #dev.off()
}

names(rhoPlotSet)

firstPhen = sapply(strsplit(names(rhoPlotSet), "-"), `[`, 1)
secondPhen = sapply(strsplit(names(rhoPlotSet), "-"), `[`, 3)

matchingPlots = which(firstPhen == secondPhen)

liamComaprePlots = rhoPlotSet[matchingPlots]
names(liamComaprePlots)

#Gene
CHPlot = liamComaprePlots[1]
CHPlot = CHPlot[[1]]

HIPlot = liamComaprePlots[2]
HIPlot = HIPlot[[1]]

HVPlot = liamComaprePlots[4]
HVPlot = HVPlot[[1]]


herbivoreCompare = grid.arrange(CHPlot, HIPlot, HVPlot, nrow = 2)

#geneComparePlot = herbivoreCompare
genePvalComparePlot = herbivoreCompare



#GO
CHPlot = liamComaprePlots[3]
CHPlot = CHPlot[[1]]

HIPlot = liamComaprePlots[1]
HIPlot = HIPlot[[1]]

HVPlot = liamComaprePlots[2]
HVPlot = HVPlot[[1]]


herbivoreCompare = grid.arrange(CHPlot, HIPlot, HVPlot, nrow = 2)

#GoComparePlot = herbivoreCompare
GoPvalComparePlot = herbivoreCompare


plot(geneComparePlot)
plot(GoComparePlot)


# ------------------------------------------------------------------
# ---------- Looking into the specific lost GO categories 

GoMetaCombined
names(GoMetaCombined)
which(GoMetaCombined$`HI-HV-Overlap` & !GoMetaCombined$`HI-HV-CH-Overlap`)

# --- Finding changed pathways 

HiDifference = which(GoMetaCombined$`HI-significant` != GoMetaCombined$`HI-significant-Liam`)
HvDifference = which(GoMetaCombined$`HV-significant` != GoMetaCombined$`HV-significant-Liam`)
ChDifference = which(GoMetaCombined$`CH-significant` != GoMetaCombined$`CH-significant-Liam`)






nonliamOnlyHI = length(which(GoMetaCombined[HiDifference, c(which(colnames(GoMetaCombined) == "HI-significant"),which(colnames(GoMetaCombined) == "HI-significant-Liam"))][,1]))
liamOnlyHI = length(which(GoMetaCombined[HiDifference, c(which(colnames(GoMetaCombined) == "HI-significant"),which(colnames(GoMetaCombined) == "HI-significant-Liam"))][,2]))
HiLoss = nonliamOnlyHI-liamOnlyHI
nonliamOnlyHV = length(which(GoMetaCombined[HiDifference, c(which(colnames(GoMetaCombined) == "HV-significant"),which(colnames(GoMetaCombined) == "HV-significant-Liam"))][,1]))
liamOnlyHV = length(which(GoMetaCombined[HiDifference, c(which(colnames(GoMetaCombined) == "HV-significant"),which(colnames(GoMetaCombined) == "HV-significant-Liam"))][,2]))
HvLoss = nonliamOnlyHV-liamOnlyHV
nonliamOnlyCH = length(which(GoMetaCombined[HiDifference, c(which(colnames(GoMetaCombined) == "CH-significant"),which(colnames(GoMetaCombined) == "CH-significant-Liam"))][,1]))
liamOnlyCH = length(which(GoMetaCombined[HiDifference, c(which(colnames(GoMetaCombined) == "CH-significant"),which(colnames(GoMetaCombined) == "CH-significant-Liam"))][,2]))
ChLoss = nonliamOnlyCH-liamOnlyCH

{
cat(paste0("Number of nonliam-only HI: \n", length(which(GoMetaCombined[HiDifference, c(which(colnames(GoMetaCombined) == "HI-significant"),which(colnames(GoMetaCombined) == "HI-significant-Liam"))][,1])), "\n"))
cat(paste0("Number of liam-only HI: \n", length(which(GoMetaCombined[HiDifference, c(which(colnames(GoMetaCombined) == "HI-significant"),which(colnames(GoMetaCombined) == "HI-significant-Liam"))][,2])), "\n"))
cat(paste0("Number of nonliam-only HV: \n", length(which(GoMetaCombined[HiDifference, c(which(colnames(GoMetaCombined) == "HV-significant"),which(colnames(GoMetaCombined) == "HV-significant-Liam"))][,1])), "\n"))
cat(paste0("Number of liam-only HV: \n", length(which(GoMetaCombined[HiDifference, c(which(colnames(GoMetaCombined) == "HV-significant"),which(colnames(GoMetaCombined) == "HV-significant-Liam"))][,2])), "\n"))
cat(paste0("Number of nonliam-only CH: \n", length(which(GoMetaCombined[HiDifference, c(which(colnames(GoMetaCombined) == "CH-significant"),which(colnames(GoMetaCombined) == "CH-significant-Liam"))][,1])), "\n"))
cat(paste0("Number of liam-only CH: \n", length(which(GoMetaCombined[HiDifference, c(which(colnames(GoMetaCombined) == "CH-significant"),which(colnames(GoMetaCombined) == "CH-significant-Liam"))][,2])), "\n"))
}

nonliamSignificant = (which(GoMetaCombined$`HI-significant` | GoMetaCombined$`HV-significant` | GoMetaCombined$`CH-significant`))
liamSignificant = (which(GoMetaCombined$`HI-significant-Liam` | GoMetaCombined$`HV-significant-Liam` | GoMetaCombined$`CH-significant-Liam`))

length(which(liamSignificant %in% nonliamSignificant))
#260, so almost all of the ones in liam are in the nonliam 

liamMissing = nonliamSignificant[which(!nonliamSignificant %in% liamSignificant)]
length(liamMissing)


liamMissingGO = GoMetaCombined[liamMissing,]
liamMissingGO = liamMissingGO[,-grep("Overlap", colnames(liamMissingGO))]
liamMissingGO = liamMissingGO[,-grep("Driver", colnames(liamMissingGO))]
liamMissingGO = liamMissingGO[,-grep("O", colnames(liamMissingGO))]

write.csv(liamMissingGO, "Results/LiamMissingGO.csv")

# --- Opposite pathways 
oppsitePathways = GoMetaCombined[which(GoMetaCombined$`HI-HV-Overlap` & !GoMetaCombined$`HI-HV-CH-Overlap`),]

grep("Overlap", colnames(oppsitePathways))

oppsitePathways = oppsitePathways[,-grep("Overlap", colnames(oppsitePathways))]
oppsitePathways = oppsitePathways[,-grep("Driver", colnames(oppsitePathways))]
oppsitePathways = oppsitePathways[,-grep("O", colnames(oppsitePathways))]


# ------------------------------------------------------------------
# --- Figure out paths files ----- 
# ------------------------------------------------------------------

OgPaths = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalPathsFile.rds")

#Make thoeretical non-liam paths

modelType = "ER"
ancestralTrait = NULL

mainTrees = readRDS("Data/zoonomiaAllMammalsTrees.rds")

nonliamPhenotype = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalPhenotypeVector.rds")
nonliamTree = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalTree.rds")
nonliamSpeciesFilter = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeSpeciesFilter.rds")

nonliamCharpaths = char2PathsCategorical(nonliamPhenotype, mainTrees, nonliamSpeciesFilter, model = modelType, anctrait = ancestralTrait) #make a path based on the phenotype vector
nonliamTreePaths = tree2Paths(nonliamTree, mainTrees, useSpecies = nonliamSpeciesFilter, categorical = TRUE) #do not binarize; the categorical data is already contained in the phenotype tree.

all.equal(OgPaths, nonliamCharpaths)
all.equal(OgPaths, nonliamTreePaths)

LiamInference

ogliamPaths = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCategoricalPathsFile.rds")


liamPhenotype = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCategoricalPhenotypeVector.rds")
liamTree = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCategoricalTree.rds")
liamSpeciesFilter = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceSpeciesFilter.rds")

liamCharpaths = char2PathsCategorical(liamPhenotype, mainTrees, liamSpeciesFilter, model = modelType, anctrait = ancestralTrait) #make a path based on the phenotype vector
liamTreePaths = tree2Paths(liamTree, mainTrees, useSpecies = liamSpeciesFilter, categorical = TRUE) #do not binarize; the categorical data is already contained in the phenotype tree.


all.equal(ogliamPaths, liamCharpaths)
all.equal(ogliamPaths, liamTreePaths)


# ------------------------------------------------------------------
# --- Compare GO Significance thresholds ----- 
# ------------------------------------------------------------------
combinedVenn = grid.arrange(geneVenn, goVenn, nrow = 1)

pdf("Output/CategoricalInsVertivoreTreeLiamInference/Categproca;")

vennPlotName = paste0(outputFolderName, filePrefix, "VennPlots.pdf")
overlapPlotName = paste0(outputFolderName, filePrefix, "OverlapPlots.pdf")



pdf(vennPlotName, 30,15)
grid.arrange(combinedVenn)
dev.off()



saveRDS(GoSignificanceColumns, paste0(combinedGODataFilename, ".rds"))


for(i in 1:length(pairwiseSets)){
  currentSet = pairwiseSets[i]
  correlationSubsetName = gsub("-", " - ", currentSet)
  correlationPrefix = paste(substr(strsplit(currentSet, split = "-")[[1]],1,1), collapse = '')
  
  goFilename = paste0(outputFolderName, currentSet, "/", filePrefix, currentSet,"Enrichment-", geneSet, ".rds")
  currentGoData = readRDS(goFilename)[[1]]
  currentGoData$significant = currentGoData$p.adj < 0.05
  message(length(which(currentGoData$significant)))
  currentGoData$significant = currentGoData$p.adj < 0.1
  message(length(which(currentGoData$significant)))
  
  names(currentGoData) = paste0(correlationPrefix, "-", names(currentGoData))
  
  
  driverFilename = paste0(outputFolderName, currentSet, "/", filePrefix, currentSet, "GoDriverTable-", geneSet, ".rds")
  if(file.exists(driverFilename)){
    driverTable = readRDS(driverFilename)
    if(all(rownames(driverTable) == rownames(currentGoData))){
      currentGoData$Driver = driverTable[which(names(driverTable) == "Driver")][[1]]
      currentGoData$DriverNumeric = driverTable[which(names(driverTable) == "DriverNumeric")][[1]]
      
      names(currentGoData)[which(names(currentGoData) == "Driver")] = paste0(correlationPrefix, "-", "Driver")
      names(currentGoData)[which(names(currentGoData) == "DriverNumeric")] = paste0(correlationPrefix, "-", "DriverNumeric")    
      
    }
    rm(driverTable)
  }
  GOResults[[i]] = currentGoData
  names(GOResults)[i] = correlationPrefix
  
}



# ------------------------------------------------------------------
# --- Compare liam results and non-liam results ----- 
# ------------------------------------------------------------------

saveRDS(GoCombinedResults, "Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsvertivoreTreeLiamInferencecombinedGOResults-KeggReactome.rds")

noLiamGenes = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsvertivoreTreecombinedGeneResults.rds")
liamGenes = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsvertivoreTreeLiamInferencecombinedGeneResults.rds")

noLiamGO = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsvertivoreTreecombinedGOResults-KeggReactome.rds")
liamGO = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsvertivoreTreeLiamInferencecombinedGOResults-KeggReactome.rds")

liam = liamGenes
nonliam = noLiamGenes
testCol = "HI-p.adj"

compareOverlap = function(testCol, liam, nonliam){
  liamColValue = which(colnames(liam) == testCol)
  liamCol = liam[liamColValue]
  
  nonliamColValue = which(colnames(nonliam) == testCol)
  nonliamCol = nonliam[nonliamColValue]
  
  isPval = !is.logical(liamCol[[1]])
  
  if(isPval){
    plot(liamCol[[1]], nonliamCol[[1]])
    
    model = lm(nonliamCol[[1]] ~ liamCol[[1]])
    r2 = summary(model)$r.squared
    text(x = min(liamCol[[1]], na.rm = T), y = max(nonliamCol[[1]], na.rm = T), labels = paste("R² =", round(r2, 3)), pos = 4)
    
    liamCol[[1]] = liamCol[[1]] < 0.1
    nonliamCol[[1]] = nonliamCol[[1]] < 0.1
  }
  
  
  
  liamSignificant = length(which(liamCol[[1]]))
  nonliamSignificant = length(which(nonliamCol[[1]]))
  
  liamInNonliam = length(which(which(liamCol[[1]]) %in% which(nonliamCol[[1]])))
  liamOverlapPercent = liamInNonliam/liamSignificant
  
  nonliamInLiam = length(which(which(nonliamCol[[1]]) %in% which(liamCol[[1]])))
  nonliamOverlapPercent = nonliamInLiam/nonliamSignificant
  
  cat("##########\n")
  cat(testCol)
  cat("\n")
  if(isPval){
    cat("R Squared: ")
    cat(r2)
    cat("\n")
  }
  
  #cat("#\n")
  cat("Liam value: ")
  cat(liamSignificant)
  cat("\n")
  cat("Non-Liam value: ")
  cat(nonliamSignificant)
  cat("\n")
  cat("Liam non-Liam ratio: ")
  cat(nonliamSignificant/liamSignificant)
  cat("\n")
  
  #cat("# \n")
  cat("Liam in Non-Liam: ")
  cat(liamInNonliam)
  cat("     ")
  cat(liamOverlapPercent)
  cat("\n")
  
  cat("Non-Liam in Liam: ")
  cat(nonliamInLiam)
  cat("     ")
  cat(nonliamOverlapPercent)
  cat("\n")
  
}

compareOverlap("HI-significant", liamGenes, noLiamGenes)
compareOverlap("HV-significant", liamGenes, noLiamGenes)
compareOverlap("IV-significant", liamGenes, noLiamGenes)
compareOverlap("CH-significant", liamGenes, noLiamGenes)

compareOverlap("HI-p.adj", liamGenes, noLiamGenes)
compareOverlap("HV-p.adj", liamGenes, noLiamGenes)
compareOverlap("IV-p.adj", liamGenes, noLiamGenes)
compareOverlap("CH-p.adj", liamGenes, noLiamGenes)

# MInimum overlap is 74% in HV, or 69% in CH, liam is larger than non-liam. 

compareOverlap("HI-significant", liamGO, noLiamGO)
compareOverlap("HV-significant", liamGO, noLiamGO)
compareOverlap("IV-significant", liamGO, noLiamGO)
compareOverlap("CH-significant", liamGO, noLiamGO)


compareOverlap("HI-p.adj", liamGO, noLiamGO)
compareOverlap("HV-p.adj", liamGO, noLiamGO)
compareOverlap("IV-p.adj", liamGO, noLiamGO)
compareOverlap("CH-p.adj", liamGO, noLiamGO)


# ------------------------------------------------------------------
# --- I-V numbers ----- 
# ------------------------------------------------------------------

genesResults = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsvertivoreTreecombinedGeneResults.rds")

length(which(genesResults$`IV-significant`))

GoResults = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsvertivoreTreecombinedGOResults-KeggReactome.rds")
length(which(GoResults$`IV-significant`))


# ------------------------------------------------------------------
# --- Examine liam results ----- 
# ------------------------------------------------------------------

liamResultsCore = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCombinedCategoricalCorrelationFile.rds")
omnibus = liamResultsCore[[1]]
length(which(omnibus$p.adj < 0.05))

liamGoOverall = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/Overall/CategoricalInsVertivoreTreeLiamInferenceOverallEnrichment-KeggReactome.rds")

liamGoOverall = liamGoOverall[[1]]
length(which(liamGoOverall$p.adj < 0.05))

liamCarnResults = readRDS("Output/CategoricalInsVertivoreTreeCarnivoreLiamInference/CategoricalInsvertivoreTreeCarnivoreLiamInferencecombinedGeneResults.rds")

length(which(liamCarnResults$`CH-significant`))
length(which(liamCarnResults$`CO-significant`))
length(which(liamCarnResults$`HO-significant`))

liamResultsCoreCarn = readRDS("Output/CategoricalInsVertivoreTreeCarnivoreLiamInference/CategoricalInsVertivoreTreeCarnivoreLiamInferenceCombinedCategoricalCorrelationFile.rds")
omnibusCarn = liamResultsCoreCarn[[1]]
length(which(omnibusCarn$p.adj < 0.05))


liamCarnGoResults = readRDS("Output/CategoricalInsVertivoreTreeCarnivoreLiamInference/CategoricalInsvertivoreTreeCarnivoreLiamInferencecombinedGOResults-KeggReactome.rds")
length(which(liamCarnGoResults$`CH-significant`))
length(which(liamCarnGoResults$`CO-significant`))
length(which(liamCarnGoResults$`HO-significant`))


liamCarnGoOverall = readRDS("Output/CategoricalInsVertivoreTreeCarnivoreLiamInference/Overall/CategoricalInsVertivoreTreeCarnivoreLiamInferenceOverallEnrichment-KeggReactome.rds")

liamCarnGoOverall = liamCarnGoOverall[[1]]
length(which(liamCarnGoOverall$p.adj < 0.05))




carnGoResults = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsvertivoreTreecombinedGOResults-KeggReactome.rds")
length(which(carnGoResults$`CH-significant`))
length(which(carnGoResults$`CO-significant`))
length(which(carnGoResults$`HO-significant`))

noliamResultsCore = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsvertivoreTreecombinedGeneResults.rds")
length(which(noliamResultsCore$`IV-significant`))

# Add carnivory results to main liam infrence

`Carnivore - Herbivore` = liamResultsCoreCarn[[2]][1] 
names(`Carnivore - Herbivore`) = "Carnivore - Herbivore"
`Carnivore - Omnivore` = liamResultsCoreCarn[[2]][2] 
names(`Carnivore - Omnivore`) = "Carnivore - Omnivore"

correlationResults[7] = `Carnivore - Herbivore`
names(correlationResults)[7] = "Carnivore - Herbivore"
correlationResults[8] = `Carnivore - Omnivore`
names(correlationResults)[8] = "Carnivore - Omnivore"

saveRDS(correlationResults, "Output/CategoricalInsvertivoreTreeLiamInference/CategoricalInsvertivoreTreeLiamInferencePairwiseCorrelationFile.rds")
write.csv(correlationResults, "Output/CategoricalInsvertivoreTreeLiamInference/CategoricalInsvertivoreTreeLiamInferencePairwiseCorrelationFile.csv")
# ------------------------------------------------------------------
# --- Make RER PLots ----- 
# ------------------------------------------------------------------

mainTrees = readRDS("Data/zoonomiaAllMammalsTrees.rds")
pathObject = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalPathsFile.rds")
RERObject = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeRERFile.rds")


commonRERs = RERObject
colnames(commonRERs) = ZonomNameConvertVectorCommon(colnames(commonRERs), tipColumn = "ZoonomiaTip")
palette(c( "darkgreen", "darkblue","black", "red"))

postiveSelectionOverlap = c("CTRL", "CPA1", "SLC36A1", "SLC6A19", "CELA3B", "CLPS", "SLC7A9", "PNLIPRP2", "FABP1", "ABCG8", "LCT")

postiveSelectionNoOverlap = c("APOB", "APOA1", "LIPF","NPC1L1","MEP1B","SLC3A2","HK1","MGAM", 
                        "APOB", "PLPP2", "APOA1", "SLC27A4", "PLA2G1B", "SLC8A2", "SLC3A1", "DPP4", "KCNN4", "SLC3A2", "MGAM2", "HK3", "G6PC",
                        "APOB", "PLA2G5", "SCAB1", "CPA3", 
                        "APOB", "MOGAT2", "MTTP", "SLC1A5", "SLC7A8", "SLC15A1", "PRKCB", "ATP1B1", "ATP1B3",
                        "SLC3A2", "SLC1A5",  "DPP4", "APOB", "CD36", "PLPP2", 
                        "APOB", "PIK3CD", "CPB2", "KCNK5"
)
length(unique(postiveSelectionNoOverlap))
length(unique(postiveSelectionOverlap))

unique(rownames(commonRERs))

i=1
{
plotRers(commonRERs, postiveSelectionOverlap[i], pathObject, sortrers = T)
i=i+1
}

i=1
{
  plotRers(commonRERs, postiveSelectionNoOverlap[i], pathObject, sortrers = T)
  i=i+1
}

i=1
{
  plotRers(commonRERs, unique(rownames(commonRERs))[i], pathObject, sortrers = T)
  i=i+1
}
#


plotRers(commonRERs, "RAD50", pathObject, sortrers = T)

# ------------------------------------------------------------------
# --- Make carnivory results for liam  ----- 
# ------------------------------------------------------------------

liamTree = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCategoricalTree.rds")
liamPhenotype = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCategoricalPhenotypeVector.rds")
liamFilter = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceSpeciesFilter.rds")


mergedData = read.csv("Data/mergedData.csv")
mixInvVertSpecies = which(mergedData$DerekDietClassification90InsVertivoreSorting == "C-InsVertivore-Mixed")

mixSpecies = mergedData$ZoonomiaTip[mixInvVertSpecies]

names(liamCarnivoryPhenotype)[which(names(liamCarnivoryPhenotype) %in% mixSpecies)]

ZonomNameConvertVectorCommon(names(liamCarnivoryPhenotype)[which(names(liamCarnivoryPhenotype) %in% mixSpecies)], tipColumn ="ZoonomiaTip")

liamCarnivoryPhenotype = liamPhenotype
liamCarnivoryPhenotype[which(liamCarnivoryPhenotype == "Vertivore" | liamCarnivoryPhenotype == "Insectivore")] = "Carnivore"
liamCarnivoryPhenotype[which(names(liamCarnivoryPhenotype) %in% mixSpecies)] = "Carnivore"

liamTreeCarnivory = liamTree
liamTreeCarnivory$edge.length[which(liamTreeCarnivory$edge.length %in% c(1))] = 5
liamTreeCarnivory$edge.length[which(liamTreeCarnivory$edge.length %in% c(3))] = 6
liamTreeCarnivory$edge.length[which(liamTreeCarnivory$edge.length %in% c(2,4))] = 1
liamTreeCarnivory$edge.length[which(liamTreeCarnivory$edge.length %in% c(5))] = 2
liamTreeCarnivory$edge.length[which(liamTreeCarnivory$edge.length %in% c(6))] = 3


speciesFilter = liamFilter
phenotypeVector = liamCarnivoryPhenotype
categoricalTree = liamTreeCarnivory
outputFolderName = "Output/CategoricalInsVertivoreTreeCarnivoreLiamInference/"
filePrefix = "CategoricalInsVertivoreTreeCarnivoreLiamInference"
phenotypeVectorFilename = paste(outputFolderName, filePrefix, "CategoricalPhenotypeVector.rds",sep="") #make a filename based on the prefix
saveRDS(liamCarnivoryPhenotype, file = phenotypeVectorFilename)                        #save the phenotype vector

speciesFilterFilename = paste(outputFolderName, filePrefix, "SpeciesFilter.rds",sep="") #set a filename for the species filter based on the prefix 
saveRDS(speciesFilter, file = speciesFilterFilename)                          #save that as the species filter




spreadSheetLocation = "Data/mergedData.csv"
nameColumn = "ZoonomiaTip"
mainTrees = readRDS("data/zoonomiaAllMammalsTrees.rds")
commonMainTrees = mainTrees
commonMainTrees$masterTree = ZoonomTreeNameToCommon(commonMainTrees$masterTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)
commonPhenotypeVector = phenotypeVector
names(commonPhenotypeVector) = ZonomNameConvertVectorCommon(names(commonPhenotypeVector), annotationLocation = spreadSheetLocation, tipColumn = nameColumn)
commonSpeciesFilter = ZonomNameConvertVectorCommon(speciesFilter, annotationLocation = spreadSheetLocation, tipColumn = nameColumn)
commonCategoricalTree = ZoonomTreeNameToCommon(categoricalTree, tipCol = "ZoonomiaTip")



treeImageFilename = paste(outputFolderName, filePrefix, "CategoricalTree.pdf", sep="") #make a filename based on the prefix
palette(c("red", "darkgreen", "black"))

pdf(treeImageFilename, height = length(phenotypeVector)/18, width = 10)                     #make a pdf to store the plot, sized based on tree size

plotTreeCategorical(commonCategoricalTree, c("Carnivore", "Herbivore", "Omnivore"), master = commonMainTrees$masterTree)
plotTreeCategorical(categoricalTree, c("Carnivore", "Herbivore", "Omnivore" ), master = mainTrees$masterTree)

#plotTreeCategorical(commonCategoricalTreeExtraTip, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = commonMasterAdded)
#plotTreeCategorical(categoricalTreeExtraTip, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = masterTreeAdded)
dev.off()  

categoricalTreeFilename = paste(outputFolderName, filePrefix, "CategoricalTree.rds", sep="") #make a filename based on the prefix
saveRDS(categoricalTree, categoricalTreeFilename)                               #save the tree
categoricalCommonTreeFilename = paste(outputFolderName, filePrefix, "CategoricalCommonTree.rds", sep="") #make a filename based on the prefix
saveRDS(commonCategoricalTree, categoricalCommonTreeFilename)


pathsFilename = paste(outputFolderName, filePrefix, "CategoricalPathsFile.rds", sep= "") #make a filename based on the prefix
paths = char2PathsCategorical(phenotypeVector, mainTrees, speciesFilter, model = modelType, anctrait = ancestralTrait) #make a path based on the phenotype vector
saveRDS(paths, file = pathsFilename)                                            #save the path 



# ------------------------------------------------------------------
# --- Emily CC figure ----- 
# ------------------------------------------------------------------

getPermsBinary=function(numperms, fg_vec, sisters_list, root_sp, RERmat, trees, mastertree, permmode="cc", method="k", min.pos=2, trees_list=NULL, calculateenrich=F, annotlist=NULL){
  pathvec = foreground2Paths(fg_vec, trees, clade="all",plotTree=F)
  col_labels = colnames(trees$paths)
  names(pathvec) = col_labels
  
  message("As of RERConverge [X.xx], permulation functions have been updated. Old versions have been moved to ccLegacy and ssmLegacy.")
  if(permmode=="cc"){
    print("Running CC permulation. sisters_list is required only for enrichments, otherwise sisters_list = NA is sufficient.")
    
    print("Generating permulated trees")
    
    # --- new code since legacy method; switching to categorical function to infer phenotype tree --
    #covert fg_vec to a categorical phenotypeVector
    phenotypeVector = rep(0, length(trees$masterTree$tip.label))
    names(phenotypeVector) = trees$masterTree$tip.label
    phenotypeVector[names(phenotypeVector) %in% fg_vec] = 1
    
    
    permulationData = categoricalPermulations(trees, phenotypeVector, rm = "ER", rp = "auto", ntrees = numperms)
    
    
    permulatedTrees = lapply(permulationData$trees, function(x) {
      tr = trees$masterTree
      tr$edge.length = c(x$tips, x$nodes)[tr$edge[,2]]
      tr$edge.length = tr$edge.length-1
      names(tr$edge.length) = NULL
      #tree2Paths(tr, treesObj, categorical = TRUE, useSpecies = names(phenvals))
      tr
    })
    permulated.binphens = list(permulatedTrees)
    #----
    #permulated.binphens = generatePermulatedBinPhen(trees$masterTree, numperms, trees, root_sp, fg_vec, sisters_list, pathvec, permmode="cc")
    permulated.fg = mapply(getForegroundsFromBinaryTree, permulated.binphens[[1]])
    permulated.fg.list = as.list(data.frame(permulated.fg))
    phenvec.table = mapply(foreground2Paths,permulated.fg.list,MoreArgs=list(treesObj=trees,clade="all"))
    phenvec.list = lapply(seq_len(ncol(phenvec.table)), function(i) phenvec.table[,i])
    
    print("Calculating correlations")
    corMatList = lapply(phenvec.list, correlateWithBinaryPhenotype, RERmat=RERmat)
    
    #make enrich list/matrices to fill
    permPvals=data.frame(matrix(ncol=numperms, nrow=nrow(RERmat)))
    rownames(permPvals)=rownames(RERmat)
    permRhovals=data.frame(matrix(ncol=numperms, nrow=nrow(RERmat)))
    rownames(permRhovals)=rownames(RERmat)
    permStatvals=data.frame(matrix(ncol=numperms, nrow=nrow(RERmat)))
    rownames(permStatvals)=rownames(RERmat)
    
    for (i in 1:length(corMatList)){
      permPvals[,i] = corMatList[[i]]$P
      permRhovals[,i] = corMatList[[i]]$Rho
      permStatvals[,i] = sign(corMatList[[i]]$Rho)*-log10(corMatList[[i]]$P)
    }
    
  }
  else if (permmode=="ccLegacy"){
    print("Running CC Legacy permulation")
    
    print("Generating permulated trees")
    permulated.binphens = generatePermulatedBinPhen(trees$masterTree, numperms, trees, root_sp, fg_vec, sisters_list, pathvec, permmode="cc")
    permulated.fg = mapply(getForegroundsFromBinaryTree, permulated.binphens[[1]])
    permulated.fg.list = as.list(data.frame(permulated.fg))
    phenvec.table = mapply(foreground2Paths,permulated.fg.list,MoreArgs=list(treesObj=trees,clade="all"))
    phenvec.list = lapply(seq_len(ncol(phenvec.table)), function(i) phenvec.table[,i])
    
    print("Calculating correlations")
    corMatList = lapply(phenvec.list, correlateWithBinaryPhenotype, RERmat=RERmat)
    
    #make enrich list/matrices to fill
    permPvals=data.frame(matrix(ncol=numperms, nrow=nrow(RERmat)))
    rownames(permPvals)=rownames(RERmat)
    permRhovals=data.frame(matrix(ncol=numperms, nrow=nrow(RERmat)))
    rownames(permRhovals)=rownames(RERmat)
    permStatvals=data.frame(matrix(ncol=numperms, nrow=nrow(RERmat)))
    rownames(permStatvals)=rownames(RERmat)
    
    for (i in 1:length(corMatList)){
      permPvals[,i] = corMatList[[i]]$P
      permRhovals[,i] = corMatList[[i]]$Rho
      permStatvals[,i] = sign(corMatList[[i]]$Rho)*-log10(corMatList[[i]]$P)
    }
    
  } else if (permmode=="ssm"){
    print("Running SSM permulation. sisters_list is required only for enrichments, otherwise sisters_list = NA is sufficient.")
    
    if (is.null(trees_list)){
      trees_list = trees$trees
    }
    
    RERmat = RERmat[match(names(trees_list), rownames(RERmat)),]
    
    print("Generating permulated trees")
    permulated.binphens = generatePermulatedBinPhenSSMBatched(trees_list,numperms,trees,root_sp,fg_vec,sisters_list,pathvec,permmode="ssm")
    
    # Get species membership of the trees
    df.list = lapply(trees_list,getSpeciesMembershipStats,masterTree=mastertree,foregrounds=fg_vec)
    df.converted = data.frame(matrix(unlist(df.list), nrow=length(df.list), byrow=T),stringsAsFactors=FALSE)
    attr = attributes(df.list[[1]])
    col_names = attr$names
    attr2 = attributes(df.list)
    row_names = attr2$names
    
    colnames(df.converted) = col_names
    rownames(df.converted) = row_names
    
    df.converted$num.fg = as.integer(df.converted$num.fg)
    df.converted$num.spec = as.integer(df.converted$num.spec)
    
    spec.members = df.converted$spec.members
    
    # Group gene trees based on the similarity of their species membership
    grouped.trees = groupTrees(spec.members)
    ind.unique.trees = grouped.trees$ind.unique.trees
    ind.unique.trees = unlist(ind.unique.trees)
    ind.tree.groups = grouped.trees$ind.tree.groups
    
    # For each unique tree, produce a permuted tree. We already have this function, but we need a list of trees to feed in.
    unique.trees = trees_list[ind.unique.trees]
    
    # precompute clade mapping for each unique tree
    unique.map.list = mapply(matchAllNodesClades,unique.trees,MoreArgs=list(treesObj=trees))
    
    # calculate paths for each permulation
    unique.permulated.binphens = permulated.binphens[ind.unique.trees]
    unique.permulated.paths = calculatePermulatedPaths_apply(unique.permulated.binphens,unique.map.list,trees)
    
    permulated.paths = vector("list", length = length(trees_list))
    for (j in 1:length(permulated.paths)){
      permulated.paths[[j]] = vector("list",length=numperms)
    }
    for (i in 1:length(unique.permulated.paths)){
      ind.unique.tree = ind.unique.trees[i]
      ind.tree.group = ind.tree.groups[[i]]
      unique.path = unique.permulated.paths[[i]]
      for (k in 1:length(ind.tree.group)){
        permulated.paths[[ind.tree.group[k]]] = unique.path
      }
    }
    attributes(permulated.paths)$names = row_names
    
    print("Calculating correlations")
    RERmat.list = lapply(seq_len(nrow(RERmat[])), function(i) RERmat[i,])
    corMatList = mapply(calculateCorPermuted,permulated.paths,RERmat.list)
    permPvals = extractCorResults(corMatList,numperms,mode="P")
    rownames(permPvals) = names(trees_list)
    permRhovals = extractCorResults(corMatList,numperms,mode="Rho")
    rownames(permRhovals) = names(trees_list)
    permStatvals = sign(permRhovals)*-log10(permPvals)
    rownames(permStatvals) = names(trees_list)
    
  } else if (permmode=="ssmLegacy"){
    print("Running SSM Legacy permulation")
    
    if (is.null(trees_list)){
      trees_list = trees$trees
    }
    
    RERmat = RERmat[match(names(trees_list), rownames(RERmat)),]
    
    print("Generating permulated trees")
    permulated.binphens = generatePermulatedBinPhenSSMBatched(trees_list,numperms,trees,root_sp,fg_vec,sisters_list,pathvec,permmode="ssmLegacy")
    
    # Get species membership of the trees
    df.list = lapply(trees_list,getSpeciesMembershipStats,masterTree=mastertree,foregrounds=fg_vec)
    df.converted = data.frame(matrix(unlist(df.list), nrow=length(df.list), byrow=T),stringsAsFactors=FALSE)
    attr = attributes(df.list[[1]])
    col_names = attr$names
    attr2 = attributes(df.list)
    row_names = attr2$names
    
    colnames(df.converted) = col_names
    rownames(df.converted) = row_names
    
    df.converted$num.fg = as.integer(df.converted$num.fg)
    df.converted$num.spec = as.integer(df.converted$num.spec)
    
    spec.members = df.converted$spec.members
    
    # Group gene trees based on the similarity of their species membership
    grouped.trees = groupTrees(spec.members)
    ind.unique.trees = grouped.trees$ind.unique.trees
    ind.unique.trees = unlist(ind.unique.trees)
    ind.tree.groups = grouped.trees$ind.tree.groups
    
    # For each unique tree, produce a permuted tree. We already have this function, but we need a list of trees to feed in.
    unique.trees = trees_list[ind.unique.trees]
    
    # precompute clade mapping for each unique tree
    unique.map.list = mapply(matchAllNodesClades,unique.trees,MoreArgs=list(treesObj=trees))
    
    # calculate paths for each permulation
    unique.permulated.binphens = permulated.binphens[ind.unique.trees]
    unique.permulated.paths = calculatePermulatedPaths_apply(unique.permulated.binphens,unique.map.list,trees)
    
    permulated.paths = vector("list", length = length(trees_list))
    for (j in 1:length(permulated.paths)){
      permulated.paths[[j]] = vector("list",length=numperms)
    }
    for (i in 1:length(unique.permulated.paths)){
      ind.unique.tree = ind.unique.trees[i]
      ind.tree.group = ind.tree.groups[[i]]
      unique.path = unique.permulated.paths[[i]]
      for (k in 1:length(ind.tree.group)){
        permulated.paths[[ind.tree.group[k]]] = unique.path
      }
    }
    attributes(permulated.paths)$names = row_names
    
    print("Calculating correlations")
    RERmat.list = lapply(seq_len(nrow(RERmat[])), function(i) RERmat[i,])
    corMatList = mapply(calculateCorPermuted,permulated.paths,RERmat.list)
    permPvals = extractCorResults(corMatList,numperms,mode="P")
    rownames(permPvals) = names(trees_list)
    permRhovals = extractCorResults(corMatList,numperms,mode="Rho")
    rownames(permRhovals) = names(trees_list)
    permStatvals = sign(permRhovals)*-log10(permPvals)
    rownames(permStatvals) = names(trees_list)
    
  } else {
    stop("Invalid binary permulation mode.")
  }
  
  if (calculateenrich){
    realFgtree = foreground2TreeClades(fg_vec, sisters_list, trees, plotTree=F)
    realpaths = tree2PathsClades(realFgtree, trees)
    realresults = getAllCor(RERmat, realpaths, method=method, min.pos=min.pos)
    realstat =sign(realresults$Rho)*-log10(realresults$P)
    names(realstat) = rownames(RERmat)
    realenrich = fastwilcoxGMTall(na.omit(realstat), annotlist, outputGeneVals=F)
    
    #sort real enrichments
    groups=length(realenrich)
    c=1
    while(c<=groups){
      current=realenrich[[c]]
      realenrich[[c]]=current[order(rownames(current)),]
      c=c+1
    }
    #make matrices to fill
    permenrichP=vector("list", length(realenrich))
    permenrichStat=vector("list", length(realenrich))
    c=1
    while(c<=length(realenrich)){
      newdf=data.frame(matrix(ncol=numperms, nrow=nrow(realenrich[[c]])))
      rownames(newdf)=rownames(realenrich[[c]])
      permenrichP[[c]]=newdf
      permenrichStat[[c]]=newdf
      c=c+1
    }
    
    counter=1;
    while (counter <= numperms){
      stat = permStatvals[,counter]
      names(stat) = rownames(RERmat)
      enrich=fastwilcoxGMTall(na.omit(stat), annotlist, outputGeneVals=F)
      #sort and store enrichment results
      groups=length(enrich)
      c=1
      while(c<=groups){
        current=enrich[[c]]
        enrich[[c]]=current[order(rownames(current)),]
        enrich[[c]]=enrich[[c]][match(rownames(permenrichP[[c]]), rownames(enrich[[c]])),]
        permenrichP[[c]][,counter]=enrich[[c]]$pval
        permenrichStat[[c]][,counter]=enrich[[c]]$stat
        c=c+1
      }
      counter = counter+1
    }
  }
  
  if(calculateenrich){
    data=vector("list", 5)
    data[[1]]=permPvals
    data[[2]]=permRhovals
    data[[3]]=permStatvals
    data[[4]]=permenrichP
    data[[5]]=permenrichStat
    names(data)=c("corP", "corRho", "corStat", "enrichP", "enrichStat")
  } else {
    data=vector("list", 3)
    data[[1]]=permPvals
    data[[2]]=permRhovals
    data[[3]]=permStatvals
    names(data)=c("corP", "corRho", "corStat")
  }
  data
}


#Load RERconverge output from murine rodent analysis
load("../../MiscData/RERconverge_output.logRTM_binary.OUmodel.RTMspeciesOnly.rds")


fgspec<-c("Pseudomys_novaehollandiae_ABTC08140", "Pseudomys_delicatulus_U1509", "Notomys_alexis_U1308", "Notomys_fuscus_M22830", "Notomys_mitchellii_M21518", "Zyzomys_pedunculatus_Z34925", "Bandicota_indica_ABTC119185", "Nesokia_indica_ABTC117074", "Hyomys_goliath_ABTC42697", "Pseudomys_shortridgei_Z25113", "Paruromys_dominator_JAE4870", "Eropeplus_canus_NMVZ21733")


newPerms = getPermsBinary(100, fgspec, NA, myTrees$masterTree$tip.label[1], myRER, myTrees, myTrees$masterTree, permmode = "cc")
legacyPerms = getPermsBinary(100, fgspec, NA, myTrees$masterTree$tip.label[1], myRER, myTrees, myTrees$masterTree, permmode = "ccLegacy")

numPerms = 100
trees = myTrees
permTrees<-list()
permFG_list<-list()
inFG<-c()
Sys.time()
pdf("Results/Emily/emilyNewPermulationPlot.pdf", onefile=TRUE, height=8.5, width=11)
for(i in 1:numPerms){
   
  phenotypeVector = rep(0, length(trees$masterTree$tip.label))
  names(phenotypeVector) = trees$masterTree$tip.label
  phenotypeVector[names(phenotypeVector) %in% fgspec] = 1
  
  permulationData = categoricalPermulations(trees, phenotypeVector, rm = "ER", rp = "auto", ntrees = 1)
    
    
    permulatedTrees = lapply(permulationData$trees, function(x) {
      tr = trees$masterTree
      tr$edge.length = c(x$tips, x$nodes)[tr$edge[,2]]
      tr$edge.length = tr$edge.length-1
      names(tr$edge.length) = NULL
      #tree2Paths(tr, treesObj, categorical = TRUE, useSpecies = names(phenvals))
      tr
    })
  permTree = permulatedTrees[[1]]
  #permTree<-getPermsBinary(1, fgspec, NA, myTrees$masterTree$tip.label[1], myRER, myTrees, myTrees$masterTree, permmode = "cc")
  permTrees[[i]]<-permTree
  fgEdges<-permTree$edge[which(permTree$edge.length==1),2]
  permFgs<-permTree$tip.label[fgEdges]
  permFG_list[[i]]<-permFgs
  inFG<-cbind(inFG, unlist(sapply(myTrees$masterTree$tip.label, function(x) if(x %in% permFgs){1} else{0})))
}
dev.off()
Sys.time()
inFG_sums<-rowSums(inFG)
inFG_props<-inFG_sums/numPerms
shortnames<-unlist(sapply(names(inFG_sums), function(x) paste(strsplit(x, "_")[[1]][1:2], collapse="_")))
names(inFG_sums)<-shortnames
names(inFG_props)<-shortnames
pdf("Results/Emily/newCCBarplot.pdf", onefile=TRUE, height=8.5, width=11)
barplot(height=inFG_sums, main=paste("Number of times each species appeared in the foreground out of", numPerms, "permulations\nccNew Permulations"), las=2, cex.names=0.5)
barplot(height=inFG_props, main=paste("Proportion of times each species appeared in the foreground out of", numPerms, "permulations\nccNew Permulations"), las=2, cex.names=0.5, ylim=c(0,0.8))
dev.off()

# ----
permTrees<-list()
permFG_list<-list()
inFG<-c()
sisters_list<-list("clade1"=c("Pseudomys_novaehollandiae_ABTC08140","Pseudomys_delicatulus_U1509"), "clade2"=c("Notomys_alexis_U1308","Notomys_fuscus_M22830"), "clade3"=c("clade2","Notomys_mitchellii_M21518"), "clade4"=c("Bandicota_indica_ABTC119185","Nesokia_indica_ABTC117074"), "clade5"=c("Paruromys_dominator_JAE4870", "Eropeplus_canus_NMVZ21733"))
fg_vec = fgspec
root_sp = trees$masterTree$tip.label[[1]]
Sys.time()
pdf("Results/Emily/emilyLegacyPermulationPlot.pdf", onefile=TRUE, height=8.5, width=11)
for(i in 1:numPerms){
  
  permulated.binphens = generatePermulatedBinPhen(trees$masterTree, 1, trees, root_sp, fg_vec, sisters_list, pathvec, permmode="cc")
  permTree = permulated.binphens[[1]][[1]]
  #permTree<-getPermsBinary(1, fgspec, NA, myTrees$masterTree$tip.label[1], myRER, myTrees, myTrees$masterTree, permmode = "cc")
  permTrees[[i]]<-permTree
  fgEdges<-permTree$edge[which(permTree$edge.length==1),2]
  permFgs<-permTree$tip.label[fgEdges]
  permFG_list[[i]]<-permFgs
  inFG<-cbind(inFG, unlist(sapply(myTrees$masterTree$tip.label, function(x) if(x %in% permFgs){1} else{0})))
}
dev.off()
Sys.time()
#inFGNew <- apply(inFG, 2, as.numeric)
#rownames(inFGNew) = rownames(inFG)
inFG_sums<-rowSums(inFG)
inFG_props<-inFG_sums/numPerms
shortnames<-unlist(sapply(names(inFG_sums), function(x) paste(strsplit(x, "_")[[1]][1:2], collapse="_")))
names(inFG_sums)<-shortnames
names(inFG_props)<-shortnames
pdf("Results/Emily/legacyCCBarPlot.pdf", onefile=TRUE, height=8.5, width=11)
barplot(height=inFG_sums, main=paste("Number of times each species appeared in the foreground out of", numPerms, "permulations\nccLegacy Permulations"), las=2, cex.names=0.5)
barplot(height=inFG_props, main=paste("Proportion of times each species appeared in the foreground out of", numPerms, "permulations\nccLegacy Permulations"), las=2, cex.names=0.5, ylim=c(0,0.8))
dev.off()



# ------------------------------------------------------------------
# --- Looking into differecen between liam and non-liam results  ----- 
# ------------------------------------------------------------------





nonLiamResults = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsvertivoreTreecombinedGeneResults.rds")
liamResults = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsvertivoreTreeLiamInferencecombinedGeneResults.rds")

which(!rownames(liamResults) %in% rownames(nonLiamResults))

liamHI = rownames(liamResults)[which(liamResults$`HI-p.adj` <0.05)]
nonliamHI = rownames(nonLiamResults)[which(nonLiamResults$`HI-p.adj` <0.05)]

length(which(liamHI %in% nonliamHI))


nonLiamResultsGO = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsvertivoreTreecombinedGOResults-KeggReactome.rds")
liamResultsGO = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsvertivoreTreeLiamInferencecombinedGOResults-KeggReactome.rds")



liamHIGO = rownames(liamResultsGO)[which(liamResultsGO$`HI-p.adj` <0.1)]
nonliamHIGO = rownames(nonLiamResultsGO)[which(nonLiamResultsGO$`HI-p.adj` <0.1)]


length(which(liamHIGO %in% nonliamHIGO))


names(liamResults) = paste0("liam-", names(liamResults))
names(liamResultsGO) = paste0("liam-", names(liamResultsGO))



which(liamResults$`HI-significant`)


# Get significance values for genes
significanceColumns = names(liamResults)[grep("significant", names(liamResults))]
geneSignificanceResults = liamResults[, names(liamResults) %in% significanceColumns]
colSums(geneSignificanceResults, na.rm = T)

significanceColumns = names(nonLiamResults)[grep("significant", names(nonLiamResults))]
geneSignificanceResultsNL = nonLiamResults[, names(nonLiamResults) %in% significanceColumns]
colSums(geneSignificanceResultsNL, na.rm = T)
#Remove the ch column from the non-liam results
geneSignificanceResultsNL = geneSignificanceResultsNL[,-1]


matchGeneSignificant = geneSignificanceResults == geneSignificanceResultsNL
matchGeneSignificant[!geneSignificanceResults & !geneSignificanceResultsNL] = NA
totalSignificant = colSums(!is.na(matchGeneSignificant))
sharedSignificant = colSums(matchGeneSignificant, na.rm = T)

sharedSignificant/totalSignificant 



significanceColumnsGO = names(liamResultsGO)[grep("significant", names(liamResultsGO))]
geneSignificanceResultsGO = liamResultsGO[, names(liamResultsGO) %in% significanceColumnsGO]
colSums(geneSignificanceResultsGO, na.rm = T)

significanceColumnsGO = names(nonLiamResultsGO)[grep("significant", names(nonLiamResultsGO))]
geneSignificanceResultsGONL = nonLiamResultsGO[, names(nonLiamResultsGO) %in% significanceColumnsGO]
colSums(geneSignificanceResultsGONL, na.rm = T)

geneSignificanceResultsGONL = geneSignificanceResultsGONL[,c(1,4,2,3,5,6)]

matchGeneSignificantGO = geneSignificanceResultsGO == geneSignificanceResultsGONL
matchGeneSignificantGO[!geneSignificanceResultsGO & !geneSignificanceResultsGONL] = NA
totalSignificantGO = colSums(!is.na(matchGeneSignificantGO))
sharedSignificantGO = colSums(matchGeneSignificantGO, na.rm = T)

sharedSignificantGO/totalSignificantGO 



# ------------------------------------------------------------------
# --- Making liam vs non-liam tree plot  ----- 
# ------------------------------------------------------------------

liamTree = readRDS("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalInsVertivoreTreeLiamInferenceCategoricalTree.rds")
nonliamTree = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalTree.rds")

all.equal(liamTree$edge, nonliamTree$edge)
which(! liamTree$edge == nonliamTree$edge)



which(liamTree$edge.length != nonliamTree$edge.length)
length(which(liamTree$edge.length != nonliamTree$edge.length))

mainTrees = readRDS("Data/zoonomiaAllMammalsTrees.rds")
commonMainTrees = mainTrees
source("Src/Reu/ZoonomTreeNameToCommon.R")
nameColumn = "ZoonomiaTip"
spreadSheetLocation = "Data/MergedData.csv" 
commonMainTrees$masterTree = ZoonomTreeNameToCommon(commonMainTrees$masterTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)
commonLiamTree = ZoonomTreeNameToCommon(liamTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)
commonNonliamTree = ZoonomTreeNameToCommon(nonliamTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)

equalBranchMaster = commonMainTrees$masterTree
equalBranchMaster$edge.length = rep(1, length(equalBranchMaster$edge.length ))

edgelabelcolor = rep("black", length(liamTree$edge.length))
edgelabelcolor[which(liamTree$edge.length != nonliamTree$edge.length)] = "purple"

palette(c( "darkgreen", "darkblue","black", "red"))

pdf("results/CompareLiamTrees.pdf", height = length(liamTree$tip.label)/10, width = 20)     
par(mfrow = c(1,2))
plotTreeCategorical(commonNonliamTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = equalBranchMaster)
edgelabels(bg=NULL, frame = "none", col = edgelabelcolor, cex = 0.3)

plotTreeCategorical(commonLiamTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = equalBranchMaster)
edgelabels(bg=NULL, frame = "none", col = edgelabelcolor, cex = 0.3)
 


plotTreeCategorical(commonNonliamTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = commonMainTrees$masterTree)
edgelabels(bg=NULL, frame = "none", col = edgelabelcolor, cex = 0.3)

plotTreeCategorical(commonLiamTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = commonMainTrees$masterTree)
edgelabels(bg=NULL, frame = "none", col = edgelabelcolor, cex = 0.3)


dev.off()  



# ------------------------------------------------------------------
# --- Looking into tree sizes  ----- 
# ------------------------------------------------------------------

mainTrees = readRDS("data/zoonomiaAllMammalsTrees.rds")
mainTrees$masterTree

unprunedTree = readRDS("Output/Categorical4CategoryUnprunedTree/Categorical4CategoryUnprunedTreeCategoricalTree.rds")
prunedTree = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalTree.rds")

prunedTree

?fastwilcoxGMTall 

function (vals, annotList, alternative = "two.sided", ...) 
{
  reslist = list()
  for (n in names(annotList)) {
    reslist[[n]] = fastwilcoxGMT(vals, annotList[[n]], alternative = alternative, 
                                 ...)
    message(paste0(nrow(reslist[[n]]), " results for annotation set ", 
                   n))
  }
  reslist
}

# ------------------------------------------------------------------
# --- using liam's zero-length-added-tip ancestral infrence method  ----- 
# ------------------------------------------------------------------

# --------------
# -- section which is actually run during phenotype generation --- 

masterTree = mainTrees$masterTree

nodesToAdd = c(455, 457, 471, 650, 492)
names(nodesToAdd) = c("Mammalia", "Marsupalia", "Placentalia", "Chiroptera", "Primates")

masterTreeAdded = masterTree
for(i in 1:length(nodesToAdd)){
  M<-matchNodes(masterTree,masterTreeAdded,method="distances")
  masterTreeAdded<-bind.tip(masterTreeAdded,names(nodesToAdd)[i],edge.length=0,
                              where=M[which(M[,1]==as.numeric(nodesToAdd[i])),2])
}
mainTrees$masterTree = masterTreeAdded

phenToAdd = c("Insectivore", "Insectivore", "Insectivore", "Insectivore", "Omnivore")
names(phenToAdd) = names(nodesToAdd)

phenotypeVector = append(phenotypeVector, phenToAdd)
speciesFilter = append(speciesFilter, names(phenToAdd))


# - Make common name versions of objects (used in visualization) - 
commonMainTrees = mainTrees
commonMainTrees$masterTree = ZoonomTreeNameToCommon(commonMainTrees$masterTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)
commonPhenotypeVector = phenotypeVector
names(commonPhenotypeVector) = ZonomNameConvertVectorCommon(names(commonPhenotypeVector), annotationLocation = spreadSheetLocation, tipColumn = nameColumn)
commonSpeciesFilter = ZonomNameConvertVectorCommon(speciesFilter, annotationLocation = spreadSheetLocation, tipColumn = nameColumn)

# - Categorical Tree - 
treeImageFilename = paste(outputFolderName, filePrefix, "CategoricalTree.pdf", sep="") #make a filename based on the prefix
palette(c( "darkgreen", "darkblue","black", "red"))

pdf(treeImageFilename, height = length(phenotypeVector)/18, width = 10)                     #make a pdf to store the plot, sized based on tree size
commonCategoricalTree = char2TreeCategorical(commonPhenotypeVector, commonMainTrees, commonSpeciesFilter, model = modelType, anctrait = ancestralTrait, plot = F)
categoricalTree = char2TreeCategorical(phenotypeVector, mainTrees, speciesFilter, model = modelType, anctrait = ancestralTrait, plot = F) #use the phenotype vector to make a tree

commonCategoricalTreeExtraTip = commonCategoricalTree
categoricalTreeExtraTip = categoricalTree

commonCategoricalTree = drop.tip(commonCategoricalTree, names(nodesToAdd))
categoricalTree = drop.tip(categoricalTree, names(nodesToAdd))
mainTrees$masterTree = drop.tip(mainTrees$masterTree, names(nodesToAdd))
commonMasterAdded = commonMainTrees$masterTree
commonMainTrees$masterTree = drop.tip(commonMainTrees$masterTree, names(nodesToAdd))

plotTreeCategorical(commonCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = commonMainTrees$masterTree)
plotTreeCategorical(categoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = mainTrees$masterTree)

plotTreeCategorical(commonCategoricalTreeExtraTip, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = commonMasterAdded)
plotTreeCategorical(categoricalTreeExtraTip, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = masterTreeAdded)
dev.off()  
# ----------------------------
# comaprision of this with the primary analysis

primaryCategoricalTree = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalTree.rds")


length(primaryCategoricalTree$edge)
length(categoricalTree$edge)
all.equal(primaryCategoricalTree$edge,categoricalTree$edge)

length(which(!primaryCategoricalTree$edge == categoricalTree$edge))

all.equal(primaryCategoricalTree$Nnode,categoricalTree$Nnode)
all.equal(primaryCategoricalTree$node.label,categoricalTree$node.label)

all.equal(primaryCategoricalTree$tip.label,categoricalTree$tip.label)

primaryCategoricalTree$tip.label %in% categoricalTree$tip.label
# What this is saying is that the tips are the same, but are stored in a different order. 
# Suggested that save-loading it will reload ordering potentially casued by tip changes) 
write.tree(categoricalTree, "Output/CategoricalInsVertivoreTreeLiamInference/CategoricalTreeRaw.tree")
categoricalTreeNew = read.tree("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalTreeRaw.tree")
all.equal(categoricalTreeNew$tip.label,categoricalTree$tip.label)
#no, that didn't do it. 

categoricalTreeNew = reorder.phylo(categoricalTree)
primaryTreeNew = reorder.phylo(primaryCategoricalTree)
#also no effect

match(categoricalTree$tip.label, primaryCategoricalTree$tip.label)

all.equal(primaryCategoricalTree$edge.length,categoricalTree$edge.length)

library(gridExtra)

liamTree = plotTreeCategorical(categoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = mainTrees$masterTree)
primaryTree = plotTreeCategorical(primaryCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = mainTrees$masterTree)

sideBySide = grid.arrange(liamTree, primaryTree)

pdf(height = length(phenotypeVector)/18, width = 20) 
par(mfrow=c(1,2))
liamTree = plotTreeCategorical(categoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = mainTrees$masterTree)
primaryTree = plotTreeCategorical(primaryCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = mainTrees$masterTree)

dev.off()
# ----------
#Code writteen to learn how to do the above 

#Make a plot of the master tree with node labels to find fing nodes
pdf("Output/CategoricalInsVertivoreTree/MasterTreeWithNodeLables.pdf", width = 12, height = 24,)
plotTree(commonMaster,offset=1,direction="rightwards",
         lwd=1)

nodelabels(commonMaster$node.label, cex=0.6, col= "red", frame="none")
nodelabels(cex=0.6, col= "green", frame="none", adj = c(0, 2))

dev.off()

#Nodes to set 
nodesToAdd = c(455, 457, 471, 650, 492)
names(nodesToAdd) = c("Mammalia", "Marsupalia", "Placentalia", "Chiroptera", "Primates")

commonMasterAdded = commonMaster
for(i in 1:length(nodesToAdd)){
  M<-matchNodes(commonMaster,commonMasterAdded,method="distances")
  commonMasterAdded<-bind.tip(commonMasterAdded,names(nodesToAdd)[i],edge.length=0,
                   where=M[which(M[,1]==as.numeric(nodesToAdd[i])),2])
}

pdf("Output/CategoricalInsVertivoreTreeLiamInference/MasterTreeWithAddedNodes.pdf", width = 24, height = 24,)
par(mfrow=c(1,2))
plotTree(commonMaster,fsize=0.8,lwd=1)
nodelabels(cex=0.6)
plotTree(commonMasterAdded,fsize=0.8,lwd=1)
dev.off()

# Liam's code

## load phytools
library(phytools)
## simulate a small tree
tree<-pbtree(n=26,tip.label=LETTERS,scale=1)
## set a value of Q
Q<-matrix(c(-1,1,0,1,-2,1,0,1,-1),3,3,
          dimnames=list(letters[1:3],letters[1:3]))
## simulate data
xx<-sim.Mk(tree,Q,internal=TRUE)
xx

plotTree(tree,offset=1,direction="upwards",
         lwd=1)
pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
cols<-setNames(palette.colors(3,"Polychrome 36"),
               letters[1:3])
points(pp$xx,pp$yy,pch=16,col=cols[xx],cex=1.5)

nn<-xx[sort(sample(1:tree$Nnode+Ntip(tree),10))]
nn

x<-xx[tree$tip.label]
x

nntree<-tree
for(i in 1:length(nn)){
  M<-matchNodes(tree,nntree,method="distances")
  nntree<-bind.tip(nntree,names(nn)[i],edge.length=0,
                   where=M[which(M[,1]==as.numeric(names(nn)[i])),2])
}

par(mfrow=c(1,2))
plotTree(tree,fsize=0.8,lwd=1)
nodelabels(cex=0.6)
plotTree(nntree,fsize=0.8,lwd=1)

# ------------------------------------------------------------------
# --- looking into data for methods section  ----- 
# ------------------------------------------------------------------
maintrees = readRDS("Data/zoonomiaAllMammalsTrees.rds")
masterTree = maintrees$masterTree

otherCategoicalTree = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalTree.rds")

categoricalTree$tip.label[which(!categoricalTree$tip.label %in% otherCategoicalTree$tip.label)]

all.equal(categoricalTree, otherCategoicalTree)

length(droppedTips)

# ------------------------------------------------------------------
# --- Getting named paths for emily ----- 
# ------------------------------------------------------------------

paths = readRDS(pathsFilename)
paths = paths-1
all.equal(binaryPaths, paths)

paths != binaryPaths

# ------------------------------------------------------------------
# --- Looking into I-V signifiacnt results ----- 
# ------------------------------------------------------------------

geneSet = "KeggReactome"

combinedDataFilename = paste0(outputFolderName, filePrefix, "combinedGeneResults.rds")
combinedResults = readRDS(combinedDataFilename)

combinedGODataFilename = paste0(outputFolderName, filePrefix, "combinedGOResults-", geneSet, ".rds")
GoCombinedResults = readRDS(combinedGODataFilename)

length(which(combinedResults$`IV-significant`))

IVSigGenes = rownames(combinedResults)[which(combinedResults$`IV-significant`)]

length(which(GoCombinedResults$`IV-significant`))

GoCombinedResults[which(GoCombinedResults$`IV-significant`),c(29,30)]


inverseDNARepairGenes = which(GoCombinedResults$`HI-HV-Overlap` & !GoCombinedResults$`HI-HV-CH-Overlap`)
inverseOlfactionGene = which(rownames(GoCombinedResults) == "REACTOME_OLFACTORY_SIGNALING_PATHWAY")
inverseGenes = append(inverseDNARepairGenes, inverseOlfactionGene)
rownames(GoCombinedResults)[inverseGenes]

GoCombinedResults[inverseGenes,c(29,30)]

length(which(GoCombinedResults$`OV-significant`))
GoCombinedResults[which(GoCombinedResults$`OV-significant`),c(35,36)]

vertCandidates = which(rownames(GoCombinedResults) %in% c("REACTOME_GLUCOCORTICOID_BIOSYNTHESIS", "REACTOME_METABOLISM_OF_STEROID_HORMONES"))

sigcols = grep("significant", colnames(GoCombinedResults))

GoCombinedResults[vertCandidates,sigcols]

# ------------------------------------------------------------------
# --- Making RERPlots fro the oxidation genes using the new maintrees ----- 
# ------------------------------------------------------------------
library(data.table)
source("Src/Reu/ZoonomTreeNameToCommon.R")
mainTrees = readRDS("data/categoricalInsVertivoreMaintrees.rds")
mainTrees$masterTree

filePrefix = "CategoricalSlimMainInsVertivoreTree"
outputFolderName = "Output/CategoricalSlimMainInsVertivoreTree/"
cat4phenotypeSet = c("Herbivore", "Insectivore",  "Omnivore", "Vertivore")
cat4colorset = c( "darkgreen", "darkblue","black", "red")

RERFileName = paste(outputFolderName, filePrefix, "RERFile.rds", sep= "")       #Set a filename for the RERs based on the prefix
cat4RERObject = readRDS(RERFileName)                                              #Use the existing ones
PathsFilename = paste(outputFolderName, filePrefix, "CategoricalPathsFile.rds", sep= "")  
pathObject = readRDS(PathsFilename)

commonRERs = cat4RERObject
colnames(commonRERs) = ZonomNameConvertVectorCommon(colnames(commonRERs), tipColumn = "ZoonomiaTip")


palette(c( "darkgreen", "darkblue","black", "red"))


genesOfNote = c("EHHADH", "ECI2", "ACADM", "ACOT12", "ACAD11", "ACOT13", "ACAA2")
k=1
plotRers(commonRERs, genesOfNote[k], pathObject)

pdf("results/temp.pdf", 20, 20)
treePlotRers(mainTrees, cat4RERObject, genesOfNote[k], phenv = pathObject)
dev.off()


mainTrees$trees$EHHADH$tip.label



categoricalTree = readRDS(paste(outputFolderName, filePrefix, "CategoricalTree.rds", sep= ""))

trimmedMainTrees = mainTrees
trimmedMasterTree = drop.tip(mainTrees$masterTree, which(!mainTrees$masterTree$tip.label %in% categoricalTree$tip.label))

trimmedMainTrees$masterTree = trimmedMasterTree
i = 1
for(i in 1:length(trimmedMainTrees$trees)){
  currentTree = trimmedMainTrees$trees[[i]]
  trimmedCurrentTree = drop.tip(currentTree, which(!currentTree$tip.label %in% categoricalTree$tip.label))
  trimmedMainTrees$trees[[i]] = trimmedCurrentTree
}


RERTree = returnRersAsTree(trimmedMainTrees, cat4RERObject, genesOfNote[k])


# -- Switching to making a new maintrees object fro returnRERs that's pruned to the right size. 

write.tree(trimmedMasterTree, "Results/InsVertMasterTree.tree")


# ------------------------------------------------------------------
# ---  -Making scatterplots of the DNA repair genes---- 
# ------------------------------------------------------------------
#Run alongside the MakeOVerlapFigure main script 
combinedDataFilename = paste0(outputFolderName, filePrefix, "combinedGeneResults.rds")
combinedResults = readRDS(combinedDataFilename)

combinedGODataFilename = paste0(outputFolderName, filePrefix, "combinedGOResults-", geneSet, ".rds")
GoCombinedResults = readRDS(combinedGODataFilename)

DNARepairGenesFilename = "Results/dnaRepairGenes.rds"
genesets = readRDS(DNARepairGenesFilename)


rhoValuesSpecific = rhoValues
specificGenes = which(rownames(rhoValues) %in% genesets$combined)

rhoValuesSpecific = rhoValues[specificGenes, ]

library(ggrepel)

for(i in 1:length(rhoValuesSpecific)){
  xName = names(rhoValuesSpecific)[i]
  if(i+1 <= length(rhoValuesSpecific)){
    for(j in (i+1):length(rhoValuesSpecific)){
      yName = names(rhoValuesSpecific)[j]
      
      rhoCorrellPlot = ggplot(rhoValuesSpecific, aes(x = .data[[xName]], y = .data[[yName]])) + 
        geom_point() + geom_pointdensity() + scale_color_viridis()
      
      
      denstiyScaleValue = ggplot_build(rhoCorrellPlot)$plot$scales$scales[[1]]$get_limits()[2]
      densityScaleSet = append(densityScaleSet, denstiyScaleValue)
      rm(rhoCorrellPlot)
    }
  }
}
densityScale = c(1, max(densityScaleSet))


rhoPlotSet = list()
netIndex= 0
for(i in 1:length(rhoValuesSpecific)){
  xName = names(rhoValuesSpecific)[i]
  if(i <= length(rhoValuesSpecific)){
    if(bothAxis){jStart = 1}else{jStart = i+1}
    for(j in (jStart):length(rhoValuesSpecific)){
      yName = names(rhoValuesSpecific)[j]
      yLabel =  paste0(replacePrefixWithName(addDashes(gsub(corrleationColumnType, "", yName))), " Dunn Z Statistic")
      xLabel =  paste0(replacePrefixWithName(addDashes(gsub(corrleationColumnType, "", xName))), " Dunn Z Statistic")
      
      rhoCorrellPlot = ggplot(rhoValuesSpecific, aes(x = .data[[xName]], y = .data[[yName]])) + 
        geom_point() + geom_pointdensity() + scale_color_viridis(name = "Density of genes", limits = densityScale) + 
        stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),formula = y ~ x,parse = TRUE, size = 6) +
        theme_classic()+
        geom_text_repel(aes(label = rownames(rhoValuesSpecific)), size = 3)+
        geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # vertical line at x=0
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # horizontal line at y=0
        xlab(xLabel) + ylab(yLabel)+
        theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))
      
      netIndex = netIndex +1
      rhoPlotSet[[netIndex]] = rhoCorrellPlot
      names(rhoPlotSet)[netIndex] = paste(xName, yName, sep="-")
      rm(rhoCorrellPlot)
    }
  }
}

rhoPlotSet$`HI-Rho-HV-Rho`

png("Results/DNARepairGenes.png", 1200,1200)
rhoPlotSet$`HI-Rho-HV-Rho`
dev.off()
# ------------------------------------------------------------------
# --- Making RERPlots fro the oxidation genes ----- 
# ------------------------------------------------------------------
library(data.table)
source("Src/Reu/ZoonomTreeNameToCommon.R")
mainTrees = readRDS("data/zoonomiaAllMammalsTrees.rds")


filePrefix = "CategoricalInsVertivoreTree"
outputFolderName = "Output/CategoricalInsVertivoreTree/"
cat4phenotypeSet = c("Herbivore", "Insectivore",  "Omnivore", "Vertivore")
cat4colorset = c( "darkgreen", "darkblue","black", "red")

RERFileName = paste(outputFolderName, filePrefix, "RERFile.rds", sep= "")       #Set a filename for the RERs based on the prefix
cat4RERObject = readRDS(RERFileName)                                              #Use the existing ones
PathsFilename = paste(outputFolderName, filePrefix, "CategoricalPathsFile.rds", sep= "")  
pathObject = readRDS(PathsFilename)

commonRERs = cat4RERObject
colnames(commonRERs) = ZonomNameConvertVectorCommon(colnames(commonRERs), tipColumn = "ZoonomiaTip")


palette(c( "darkgreen", "darkblue","black", "red"))
plotRers(commonRERs, "ACAA2", pathObject)

genesOfNote = c("EHHADH", "ECI2", "ACADM", "ACOT12", "ACAD11", "ACOT13", "ACAA2")

k=1

pdf("results/temp.pdf", 20, 20)
treePlotRers(mainTrees, cat4RERObject, genesOfNote[k], phenv = pathObject)
dev.off()


categoricalTree = readRDS(paste(outputFolderName, filePrefix, "CategoricalTree.rds", sep= ""))

trimmedMainTrees = mainTrees
trimmedMasterTree = drop.tip(mainTrees$masterTree, which(!mainTrees$masterTree$tip.label %in% categoricalTree$tip.label))

trimmedMainTrees$masterTree = trimmedMasterTree
i = 1
for(i in 1:length(trimmedMainTrees$trees)){
  currentTree = trimmedMainTrees$trees[[i]]
  trimmedCurrentTree = drop.tip(currentTree, which(!currentTree$tip.label %in% categoricalTree$tip.label))
  trimmedMainTrees$trees[[i]] = trimmedCurrentTree
}


RERTree = returnRersAsTree(trimmedMainTrees, cat4RERObject, genesOfNote[k])


# -- Switching to making a new maintrees object fro returnRERs that's pruned to the right size. 

write.tree(trimmedMasterTree, "Results/InsVertMasterTree.tree")

# ------------------------------------------------------------------
# ---  Getting gene list of the pathways in opposite directions  ----- 
# ------------------------------------------------------------------
library(RERconverge)
library("tools")
source("Src/Reu/cmdArgImport.R")
library(xlsx)

pathwayNames = c("REACTOME_HDR_THROUGH_SINGLE_STRAND_ANNEALING_SSA", "REACTOME_HDR_THROUGH_HOMOLOGOUS_RECOMBINATION_HRR", "REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR", "REACTOME_HOMOLOGY_DIRECTED_REPAIR", "REACTOME_DISEASES_OF_DNA_REPAIR")

gmtFileName = "Data/KeggReactome.gmt"
gmtFile = read.gmt(gmtFileName)

which(gmtFile$geneset.names %in% pathwayNames)

genesets = gmtFile$genesets[which(gmtFile$geneset.names %in% pathwayNames)]
names(genesets) = gmtFile$geneset.names[which(gmtFile$geneset.names %in% pathwayNames)]

unlist(genesets[1:5])

genesets$combined = list(unlist(genesets[1:5]))
genesets$combined = unlist(genesets[1:5])
genesets$combined = unique(genesets$combined)

DNARepairGenesFilename = "Results/dnaRepairGenes.rds"
saveRDS(genesets, DNARepairGenesFilename)

lines <- sapply(genesets, function(x) paste(x, collapse = "\t"))

writeLines(lines, "Results/dnaRepairGenes.txt")

names(genesets)

which(is.numeric(names(droppedTips)))
which(!is.na(as.numeric(names(droppedTips))))
length(which(!is.na(as.numeric(names(droppedTips)))))
# ------------------------------------------------------------------
# ---  Make plots using the branchlength removed mastertrees  ----- 
# ------------------------------------------------------------------
stableMaintrees = mainTrees 
mainTrees = stableMaintrees
stableMaintrees = readRDS(mainTreesLocation)

stableCommonMainTrees = stableMaintrees
stableCommonMainTrees$masterTree = ZoonomTreeNameToCommon(stableCommonMainTrees$masterTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)


mainTrees$masterTree$edge.length[1:length(mainTrees$masterTree$edge.length)] = 1

pdf(treeImageFilename, height = length(phenotypeVector)/14, width = 10)     
plotTreeCategorical(commonCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = stableCommonMainTrees$masterTree)

plotTreeCategorical(categoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = stableMaintrees$masterTree)
dev.off()  

categoricalCommonTreeFilename = paste(outputFolderName, filePrefix, "CategoricalCommonTree.rds", sep="") #make a filename based on the prefix
saveRDS(commonCategoricalTree, categoricalCommonTreeFilename)

plotTreeCategorical(categoricalTree, c("Carnivore", "Herbivore", "Omnivore"), master = stableMaintrees$masterTree)

plotTreeCategorical(commonCategoricalTree, c("Carnivore", "Herbivore", "Omnivore"), master = stableCommonMainTrees$masterTree)

pdf(treeImageFilename, height = length(phenotypeVector)/14, width = 10)     
plotTreeCategorical(commonCategoricalTree, c("Background", "Carnivore"), master = stableCommonMainTrees$masterTree)

plotTreeCategorical(categoricalTree, c("Background", "Carnivore"), master = stableMaintrees$masterTree)
dev.off()  



# ------------------------------------------------------------------
# ---  Make subset tree for tyler ----- 
# ------------------------------------------------------------------

tylerSpecies = read.csv("Results/toMichael.csv")


tipsToKeep = tylerSpecies$fa


commonCategoricalTree

tipsToKeep = ZonomNameConvertVectorCommon(tipsToKeep, tipColumn = "ZoonomiaTip")
tipsToKeep = commonCategoricalTree$tip.label[which(1:64 %% 2 ==0)]
tipsToDrop = commonCategoricalTree$tip.label[!commonCategoricalTree$tip.label %in% tipsToKeep]
commonCategoricalTreePruned = drop.tip(commonCategoricalTree, tipsToDrop)

commonCategoricalTree = commonCategoricalTreePruned

ggTreeOut = ggtree(commonCategoricalTree) +scale_color_manual(values=palette()) 
ggTreeOut = ggTreeOut %<+% edge + aes(color=CategorylengthChar)
ggTreeOut = ggTreeOut %<+% tip_data 
ggTreeOut$data$label = paste(ggTreeOut$data$label, "-", ggTreeOut$data$node, sep="")
ggTreeOut = ggTreeOut + geom_tiplab()
#ggTreeOut = ggTreeOut + geom_tiplab(geom = "phylopic", aes(image = uuid))
#ggTreeOut + geom_phylopic(aes(uuid = uuid), color = "black", alpha = 1, size = 0.08)
ggTreeOut




# ------------------------------------------------------------------
# ---  Confirm Rho meaning for Categorical Diet ----- 
# ------------------------------------------------------------------


results = readRDS("OUtput/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCombinedCategoricalCorrelationFile.rds")
RERObject = readRDS("Output/CategoricalBinaryCarnivoreTree/CategoricalBinaryCarnivoreTreeRERFile.rds")
phenotypeVector = readRDS("Output/CategoricalBinaryCarnivoreTree/CategoricalBinaryCarnivoreTreeCategoricalPhenotypeVector.rds")
pathsObject = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalPathsFile.rds")
mainTrees = readRDS("Data/zoonomiaAllMammalsTrees.rds")
phenotypeSet = c("Herbivore", "Insectivore", "Omnivore", "Vertivore")
colorset = c( "darkgreen", "darkblue", "black", "red")

plotRers(RERObject, "CPB1" , pathsObject,)
source("Src/Reu/rerViolinPlot.R")
rerViolinPlot(mainTrees, RERObject, pathsObject, phenotypeSet , geneOfInterest = "CPB1", colorScale = colorset)

pairwiseResults = results[[2]]

genelist = c("CPB1")
rowindex = which(rownames(pairwiseResults$`1 - 2`)%in% genelist)

pairwiseResults$`1 - 2`[rowindex,]

# ------------------------------------------------------------------
# ---  Making Demos of the plots or other main-repo addable functions ----- 
# ------------------------------------------------------------------

#MakeMasterAndGeneTreePlot 

source("Src/Reu/makeMasterAndGeneTreePlots.R")

mainTrees = readRDS("Data/zoonomiaAllMammalsTrees.rds")
RERObject = readRDS("Output/CategoricalBinaryCarnivoreTree/CategoricalBinaryCarnivoreTreeRERFile.rds")
phenotypeVector = readRDS("Output/CategoricalBinaryCarnivoreTree/CategoricalBinaryCarnivoreTreeCategoricalPhenotypeVector.rds")
geneInQuestion = "CFTR"
dropTips = T
convertNames = T
tipColumn = "ZoonomiaTip"

png("Results/MasterAndGeneExmaple.png", width = 1200, height = 1200)
makeMasterAndGeneTreePlots(mainTrees, "M6PR", RERObject, tipColumn = "ZoonomiaTip", phenotypeVector = phenotypeVector, fgcols = "orange", bgcolor = "darkgreen")
dev.off()


#RER violin plot 
source("Src/Reu/rerViolinPlot.R")

RerData = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeRERFile.rds")
PathsData = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalPathsFile.rds")
PhenotypeSet = c("Herbivore", "Insectivore",  "Omnivore", "Vertivore")
Colorset = c( "darkgreen", "darkblue","black", "red")
geneOfInterest = "CYP1A1"

png("Results/RERViolinPlotExample.png", 600, 600)
rerViolinPlot(mainTrees, RerData, pathsObject = PathsData, phenotypeSet = PhenotypeSet, colorScale = Colorset, geneOfInterest = geneOfInterest)
dev.off()

# ------------------------------------------------------------------
# --- Checking if mergedata has all hiller speices   ----- 
# ------------------------------------------------------------------
mergeData = read.csv("Data/mergedData.csv")
mainTrees = readRDS('data/zoonomiaAllMammalsTrees.rds')

mainTrees$masterTree$tip.label[which(!mainTrees$masterTree$tip.label %in% mergeData$ZoonomiaTip)]

which(mergeData)


# ------------------------------------------------------------------
# --- Making a column in mergeData that matchs the InsVertivore classification   ----- 
# ------------------------------------------------------------------

substitutions = list(
  c("C-Invertebrate-eater", "Insectivore"), c("C-InsVertivore-Insectivore", "Insectivore"),
  c("C-Herpetivore", "Vertivore"),
  c("C-Piscivore", "Vertivore"), c("C-InsVertivore-Piscivore", "Vertivore"),
  c("C-Endotherm-Carnivore", "Vertivore"), c("C-Scavenger", "Vertivore"), c("C-Nonspecific-Vertebrate-eater", "Vertivore"),
  c("C-Terrestrial-vertebrates-eater", "Vertivore"), c("C-All-vertebrate-eater", "Vertivore"), c("C-InsVertivore-Carnivore", "Vertivore"),
  c("C-InsVertivore-Mixed", "Omnivore"), 
  c("O-For Examination", "Omnivore"), c("O-Scavenger", "Omnivore"),
  c("H-Frugivore", "Herbivore"), 
  c("H-Nectarivore", "Herbivore"), 
  c("H-High-sugar-plants-Eater", "Herbivore"),
  c("H-Granivore", "Herbivore"), c("H-Nonspecific-Herbivore", "Herbivore"), 
  c("H-Low-sugar-plants-Eater", "Herbivore"), c("H-All-plants-Eater", "Herbivore"),
  c("O-Generalist", "Omnivore")
)
mergeData = read.csv("Data/mergedData.csv")

Dietvalues = mergeData$DerekDietClassification90InsVertivoreSorting

  for( i in 1:length(substitutions)){
    substitutePhenotypes = substitutions[[i]]
    message(paste("replacing", substitutePhenotypes[1], "with", substitutePhenotypes[2]))
    Dietvalues = gsub(substitutePhenotypes[1], substitutePhenotypes[2], Dietvalues)
  }

mergeData$insVertivoreDiet = Dietvalues

write.csv(mergeData, "Data/mergedData.csv", row.names = F)
# ------------------------------------------------------------------
# --- Making Updated RERConverge Explanation slide  ----- 
# ------------------------------------------------------------------
library(RERconverge)
source("Src/Reu/ZonomNameConvertMatrixCommon.R")
mainTrees = readRDS("Data/RemadeTreesAllZoonomiaSpecies.rds")
mainTrees = readRDS("data/zoonomiaAllMammalsTrees.rds")

CVHRERs = readRDS("Output/Old/CVHRemake/CVHRemakeRERFile.rds")
foregroundSpecies = readRDS("Output/Old/CVHRemake/CVHRemakeBinaryTreeForegroundSpecies.rds")
CVHPaths = readRDS("Output/Old/CVHRemake/CVHRemakePathsFile.rds")
commonRERs = ZonomNameConvertMatrixCommon(CVHRERs)


source("Src/Reu/makeMasterAndGeneTreePlots.R")

# this one is correctly sized for the most part, but the tip labels are too small (especially given the lighter orange). I don't know how to fix that -- changing the obvious values in the plotting funciton had no effect. 
png("Results/tempMasterandGenePlot.png", 575, 575)
makeMasterAndGeneTreePlots(mainTrees,"IQANK1", CVHRERs,  foregroundSpecies, correlationPlot = F, tipColumn = "manualAnnotations_FaName", fgcols = "orange", bgcolor = "darkgreen")
dev.off()

png("Results/tempCorrelationPlot.png", 420, 420)
makeMasterAndGeneTreePlots(mainTrees,"IQANK1", CVHRERs,  foregroundSpecies, correlationPlot = T, tipColumn = "manualAnnotations_FaName", fgcols = "orange", bgcolor = "darkgreen")
dev.off()


plotRers(commonRERs, "IQANK1", CVHPaths, sort = F)
plotRersNew = function (rermat = NULL, index = NULL, phenv = NULL, rers = NULL, method = "k", xlims = NULL, plot = 1, xextend = 0.2, sortrers = F) {
  {
    if (!is.null(phenv) && length(unique(phenv[!is.na(phenv)])) > 
        2) {
      categorical = TRUE
      if (method != "aov") {
        method = "kw"
      }
    }
    else {
      categorical = FALSE
    }
    if (is.null(rers)) {
      e1 = rermat[index, ][!is.na(rermat[index, ])]
      colids = !is.na(rermat[index, ])
      e1plot <- e1
      if (exists("speciesNames")) {
        names(e1plot) <- speciesNames[names(e1), ]
      }
      if (is.numeric(index)) {
        gen = rownames(rermat)[index]
      }
      else {
        gen = index
      }
    }
    else {
      e1plot = rers
      gen = "rates"
    }
    names(e1plot)[is.na(names(e1plot))] = ""
    if (!is.null(phenv)) {
      phenvid = phenv[colids]
      if (categorical) {
        fgdcor = getAllCor(rermat[index, , drop = F], phenv, 
                           method = method)[[1]]
      }
      else {
        fgdcor = getAllCor(rermat[index, , drop = F], phenv, 
                           method = method)
      }
      plottitle = paste0(gen, ": rho = ", round(fgdcor$Rho, 
                                                4), ", p = ", round(fgdcor$P, 4))
      if (categorical) {
        n = length(unique(phenvid))
        if (n > length(palette())) {
          pal = colorRampPalette(palette())(n)
        }
        else {
          pal = palette()[1:n]
        }
      }
      if (categorical) {
        df <- data.frame(species = names(e1plot), rer = e1plot, 
                         stringsAsFactors = FALSE) %>% mutate(mole = as.factor(phenvid))
      }
      else {
        df <- data.frame(species = names(e1plot), rer = e1plot, 
                         stringsAsFactors = FALSE) %>% mutate(mole = as.factor(ifelse(phenvid > 
                                                                                        0, 2, 1)))
      }
    }
    else {
      plottitle = gen
      df <- data.frame(species = names(e1plot), rer = e1plot, 
                       stringsAsFactors = FALSE) %>% mutate(mole = as.factor(ifelse(0, 
                                                                                    2, 1)))
    }
    if (sortrers) {
      df = filter(df, species != "") %>% arrange(desc(rer))
    }
    if (is.null(xlims)) {
      ll = c(min(df$rer) * 1.1, max(df$rer) + xextend)
    }
    else {
      ll = xlims
    }
  }
  if (categorical) {
    g <- ggplot(df, aes(x = rer, y = factor(species, levels = unique(ifelse(rep(sortrers, 
                                                                                nrow(df)), species[order(rer)], sort(unique(species))))), 
                        col = mole, label = species)) + scale_size_manual(values = c(1, 
                                                                                     1, 1, 1)) + geom_point(aes(size = mole)) + scale_color_manual(values = pal) + 
      scale_x_continuous(limits = ll) + geom_text(hjust = 1, 
                                                  size = 2) + ylab("Branches") + xlab("relative rate") + 
      ggtitle(plottitle) + geom_vline(xintercept = 0, linetype = "dotted") + 
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
            legend.position = "none", panel.background = element_blank(), 
            axis.text = element_text(size = 18, face = "bold", 
                                     colour = "black"), axis.title = element_text(size = 24, 
                                                                                  face = "bold"), plot.title = element_text(size = 24, 
                                                                                                                            face = "bold")) + theme(axis.line = element_line(colour = "black", 
                                                                                                                                                                             size = 1)) + theme(axis.line.y = element_blank())
  }
  else {
    g <- ggplot(df, 
                aes(x = rer, 
                    y = factor(species, levels = unique(ifelse(rep(sortrers, nrow(df)), species[order(rer)], sort(unique(species))))), 
                    col = mole, 
                    label = species
                )
    ) + 
      scale_size_manual(values = c(1, 1, 1, 1)) + 
      geom_point(aes(size = mole)) + 
      scale_color_manual(values = c("black", "blue")) + 
      scale_x_continuous(limits = ll) + 
      geom_text(hjust = "center", size = 3) + 
      ylab("Branches") + 
      xlab("relative rate") + 
      ggtitle(plottitle) + 
      geom_vline(xintercept = 0, linetype = "dotted") + 
      theme(
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
        legend.position = "none", 
        panel.background = element_blank(), 
        axis.text = element_text(size = 18, face = "bold", colour = "black"), 
        axis.title = element_text(size = 24, face = "bold"), 
        plot.title = element_text(size = 24, face = "bold")) + 
      theme(axis.line = element_line(colour = "black", size = 1)) + 
      theme(axis.line.y = element_blank())
  }
  if (plot) {
    print(g)
  }
  else {
    g
  }
}



mainTrees



# ------------------------------------------------------------------
# --- Getting trees for jack  ----- 
# ------------------------------------------------------------------
categoricalTreeFilename

categoricalTestTree = readRDS(categoricalTreeFilename)

masterTreeWithBranchLengths = stableMaintrees$masterTree
treeInMasterWithoutPhenotype = masterTreeWithBranchLengths$tip.label[!masterTreeWithBranchLengths$tip.label %in% categoricalTestTree$tip.label]
masterTreeWithBranchLengthsPruned = drop.tip(masterTreeWithBranchLengths, treeInMasterWithoutPhenotype)


togaTree = read.newick("Data/TogaTree.nwk")
treeInHillerWithoutPhenotype = togaTree$tip.label[!togaTree$tip.label %in% categoricalTestTree$tip.label]
hillerTreePruned = drop.tip(togaTree, treeInHillerWithoutPhenotype)

hillerTreePruned$edge.length
masterTreeWithBranchLengthsPruned$edge.length


write.tree(categoricalTestTree, paste0(outputFolderName, "ZoonomiaMaximalTreeCategoryBranchLengths.nwk"))
write.tree(masterTreeWithBranchLengthsPruned, paste0(outputFolderName, "ZoonomiaMaximalTreeAlignmentBranchLengths.nwk"))
write.tree(hillerTreePruned, paste0(outputFolderName, "ZoonomiaMaximalTreeHillerBranchLengths.nwk"))
saveRDS(categoricalTestTree, paste0(outputFolderName, "ZoonomiaMaximalTreeCategoryBranchLengths.rds"))
saveRDS(masterTreeWithBranchLengthsPruned, paste0(outputFolderName, "ZoonomiaMaximalTreeAlignmentBranchLengths.rds"))
saveRDS(hillerTreePruned, paste0(outputFolderName, "ZoonomiaMaximalTreeHillerBranchLengths.rds"))


# ------------------------------------------------------------------
# --- Making violin plots from the categorical data for presentaiton  ----- 
# ------------------------------------------------------------------

library(tools)
library(RERconverge)
filePrefix = "CategoricalInsVertivoreTree"
outputFolderName = "Output/CategoricalInsVertivoreTree/"
phenotypeStyle = "Categorical"
mainTreesLocation = "data/zoonomiaAllMammalsTrees.rds"
cat4phenotypeSet = c("Herbivore", "Insectivore",  "Omnivore", "Vertivore")
cat4colorset = c( "darkgreen", "darkblue","black", "red")

RERFileName = paste(outputFolderName, filePrefix, "RERFile.rds", sep= "")       #Set a filename for the RERs based on the prefix
cat4RERObject = readRDS(RERFileName)                                              #Use the existing ones

pathsFileName = paste(outputFolderName, filePrefix, phenotypeStyle, "PathsFile.rds", sep= "") #Set a filename for the pathss based on the prefix and style
cat4pathsObject = readRDS(pathsFileName)                                          #If the file already exists, use the existing one.

combinedDataFilename = paste0(outputFolderName, filePrefix, "combinedGeneResults.rds")
combinedResults = readRDS(combinedDataFilename)

#This one is for a gene that's significnat in all three predatory diets, it was used fof the volution talk 
#overlap = combinedResults[combinedResults$`CH-HI-HV-Overlap` & !is.na(combinedResults$`CH-HI-HV-Overlap`) & !combinedResults$`HO-significant`,]
#overlapGenes = rownames(overlap[order(overlap$`HI-p.adj`),])


#This one is for a gene that's significant in the two specific predatory diets but not combined carnivore. It was used in the fellowship proposal 
overlap = combinedResults[!combinedResults$`CH-HI-Overlap` & 
                            !combinedResults$`HI-HV-Overlap` & 
                            !is.na(combinedResults$`HI-HV-Overlap`) & 
                            !is.na(combinedResults$`CH-HI-HV-Overlap`) & 
                            !combinedResults$`HO-significant` &
                            combinedResults$`HI-significant`,]
overlapGenes = rownames(overlap[order(overlap$`HI-p.adj`),])


overlap2 = combinedResults[!combinedResults$`CH-HV-Overlap` & 
                            !combinedResults$`HI-HV-Overlap` & 
                            !is.na(combinedResults$`HI-HV-Overlap`) & 
                            !is.na(combinedResults$`CH-HI-HV-Overlap`) & 
                            !combinedResults$`HO-significant` &
                            combinedResults$`HV-significant`,]
overlapGenes2 = rownames(overlap2[order(overlap$`HV-p.adj`),])



filePrefix = "CategoricalPrunedCarnivoreTree"
outputFolderName = "Output/CategoricalPrunedCarnivoreTree/"
cat3phenotypeSet = c("Carnivore", "Herbivore", "Omnivore")
cat3colorset = c( "orange", "darkgreen", "black")

RERFileName = paste(outputFolderName, filePrefix, "RERFile.rds", sep= "")       #Set a filename for the RERs based on the prefix
cat3RERObject = readRDS(RERFileName)                                              #Use the existing ones

pathsFileName = paste(outputFolderName, filePrefix, phenotypeStyle, "PathsFile.rds", sep= "") #Set a filename for the pathss based on the prefix and style
cat3pathsObject = readRDS(pathsFileName)                                          #If the file already exists, use the existing one.



if(file_ext(mainTreesLocation) == "rds"){
  if(!exists("mainTrees")){mainTrees = readRDS(mainTreesLocation)}
}else{
  if(!exists("mainTrees")){mainTrees = readTrees(mainTreesLocation)} 
}





source("Src/Reu/rerViolinPlot.R")
library(gridExtra)

i = 1
{
currentGene = overlapGenes[i]
cat4Plot = rerViolinPlot(mainTrees, cat4RERObject, cat4pathsObject, cat4phenotypeSet , geneOfInterest = currentGene, colorScale = cat4colorset)
cat4Plot = cat4Plot + xlab(c("Herbivore", "Invertivore", "Omnivore", "Vertivore"))
cat3Plot = rerViolinPlot(mainTrees, cat3RERObject, cat3pathsObject, cat3phenotypeSet , geneOfInterest = currentGene, colorScale = cat3colorset)
comboPlot = grid.arrange(cat4Plot, cat3Plot, ncol = 2)
comboPlot
i = i+1
}
goodGenes = c(1, 8, 10)

png("Output/CategoricalInsVertivoreTree/Visualizations/InvertivoryUniqueGene.png", height = 382, width = 742)
print(comboPlot)
dev.off()
# ------------------------------------------------------------------
# --- Work on making a tree figure for lalitha  ----- 
# ------------------------------------------------------------------
lalithaData = read.csv("C:/Users/mit221/Downloads/SpeciesNamesAndPhenosForDE72.csv")


laltihaSpecies = lalithaData[,1]

#laltihaSpecies = unlist(laltihaSpecies)
#laltihaSpecies = laltihaSpecies[-1]
#laltihaSpecies = gsub("[0-9]+$", "", laltihaSpecies)
#laltihaSpecies = unique(laltihaSpecies)

mergeData = read.csv("Data/mergedData.csv")

for(i in 1:length(laltihaSpecies)){
laltihaSpecies[i] = gsub(" ", "_", laltihaSpecies[i]) #replace spaces with underscores 
laltihaSpecies[i] = sub('^([^_]+_[^_]+).*', '\\1', laltihaSpecies[i]) #remove anything  after a second underscore
laltihaSpecies[i] = tolower(laltihaSpecies[i])
}

lalithaData[,1] = laltihaSpecies
write.csv(lalithaData, "C:/Users/mit221/Downloads/SpeciesNamesAndPhenosForDE72.csv")

sum(!laltihaSpecies %in% mergeData$Scientific_Binomial)


ggTreeOut = ggtree(commonCategoricalTree) +scale_color_manual(values=palette()) 




# ------------------------------------------------------------------
# --- Assess driving diet of unqieu carnivory results  ----- 
# ------------------------------------------------------------------

combinedData = readRDS(paste0(combinedDataFilename, ".rds"))

carnivoryUNique = combinedData[which(combinedData$`CH-significant` & !combinedData$`HI-significant` & !combinedData$`HV-significant`),]


table(carnivoryUNique$`CH-Driver`)
table(sign(carnivoryUNique$`CH-Rho`))
# Right. I don't current have driver information because I haven't run the driver analysis for the 
# I need to see about making that. 

# quick fix for the driver analysis to retarget because this is not technically a part of the main analysis 

# ------------------------------------------------------------------
# --- Get set of genes most different between HI and HV  ----- 
# ------------------------------------------------------------------

combinedResults # from MakeOverlapFigure
combinedResultsModified = combinedResults


getComparisionDifference = function(dataframe, colOne, colTwo){
  colOneIndex = names(dataframe)[which(names(dataframe) == colOne)]
  colTwoIndex = names(dataframe)[which(names(dataframe) == colTwo)]
  
  distanceFromEqual = abs(dataframe[colOneIndex] - dataframe[colTwoIndex]) / sqrt(2)
  distanceFromEqual
}



getComparisionDifference(combinedResultsModified, "HI-Rho", "HV-Rho")
combinedResults$`HI-HV-Delta` = getComparisionDifference(combinedResultsModified, "HI-Rho", "HV-Rho")


test = combinedResultsModified[order(combinedResultsModified$`HI-HV-Delta`, decreasing = T),]


df$distance_from_y_eq_x <- abs(df$x - df$y) / sqrt(2)




# ------------------------------------------------------------------
# --- examining the GO sets in specific overlap sections ----- 
# ------------------------------------------------------------------
GoSignificanceResults # from MakeOverlapFigure

View(GoSignificanceResults)

# get pathways that show up in compionents but not carnivory 
IVnCpathways = rownames(GoSignificanceResults[which(GoSignificanceResults$`HI-significant` & GoSignificanceResults$`HV-significant` & !GoSignificanceResults$`CH-significant`),])

which(rownames(GoCombinedResults) %in% IVnCpathways)
IVnCResults = GoCombinedResults[which(rownames(GoCombinedResults) %in% IVnCpathways), ]
View(IVnCResults)
cat(IVnCpathways)

# get the pathways which are unique to vertivory 
VonlyPathways = rownames(GoSignificanceResults[which(!GoSignificanceResults$`HI-significant` & GoSignificanceResults$`HV-significant` & !GoSignificanceResults$`CH-significant`),])
IonlyPathways = rownames(GoSignificanceResults[which(GoSignificanceResults$`HI-significant` & !GoSignificanceResults$`HV-significant` & !GoSignificanceResults$`CH-significant`),])


# ------------------------------------------------------------------
# --- getting permualtion RERResult data fro rho plot creation ----- 
# ------------------------------------------------------------------

permulationIntermediate = readRDS(paste0(outputFolderName, "CategoricalInsVertivoreTreeCategoricalPermulationsIntermediates101.rds"))


permulationHIStatsValues = permulationIntermediate$Peffsize$`1 - 3`
permulationHVStatsValues = permulationIntermediate$Peffsize$`1 - 4`

permulationHIStatsCol = permulationHIStatsValues[,1]
permulationHVStatsCol = permulationHVStatsValues[,1]


rhoValuesPerm = data.frame(permulationHIStatsCol, permulationHVStatsCol)

for(i in 1:length(rhoValuesPerm)){
  xName = names(rhoValuesPerm)[i]
  if(i+1 <= length(rhoValuesPerm)){
    for(j in (i+1):length(rhoValuesPerm)){
      yName = names(rhoValuesPerm)[j]
      yLabel =  paste0(replacePrefixWithName(addDashes(gsub(corrleationColumnType, "", yName))), " Stat")
      xLabel =  paste0(replacePrefixWithName(addDashes(gsub(corrleationColumnType, "", xName))), " Stat")
      
      rhoCorrellPlot = ggplot(rhoValuesPerm, aes(x = .data[[xName]], y = .data[[yName]])) + 
        geom_point() + geom_pointdensity() + scale_color_viridis(name = "Number of nearby genes", limits = densityScale) + 
        stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),formula = y ~ x,parse = TRUE) +
        theme_classic()+
        xlab(xLabel) + ylab(yLabel)
    
      
      netIndex = netIndex +1
      rhoPlotSet[[netIndex]] = rhoCorrellPlot
      names(rhoPlotSet)[netIndex] = paste(xName, yName, sep="-")
    }
  }
}

rhoDfList = list()
permulationRhoDataframes = for(i in 1:2500){
  rhoDf = data.frame(permulationHIStatsValues[,i], permulationHVStatsValues[,i])
  rhoDfList[[i]] = rhoDf
}

rSquaredSet = NULL
for(k in 1:length(rhoDfList)){
  rhoDf = rhoDfList[[k]]
  names(rhoDf) = c("HI", "HV")
  model <- lm(HV ~ HI, data = rhoDf)
  rSquared = summary(model)$r.squared
  rSquaredSet = append(rSquaredSet, rSquared)
  #cat(rSquared, "\n")
}
pdf()
hist(rSquaredSet)
dev.off()

mean(rSquaredSet)
