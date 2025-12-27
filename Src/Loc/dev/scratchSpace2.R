a = b #This is to prevent accidental ful runs 


library(RERconverge)
?readTrees


mainTrees = readRDS("Data/zoonomiaAllMammalsTrees.rds")

togaTree = read.tree("Data/togaTree.nwk")
plotTree(togaTree)

allZoonomiaAlignments = read.table("Data/zoonomiaAllMammalsTrees.tsv")

togaTree$tip.label

mainTrees = readTrees(mainTreesLocation, masterTree = togaTree)
#saveRDS(mainTrees, "Data/zoonomiaAllMammalsTrees.rds")


mainTreesTrees = mainTrees[[1]]

alltips = character()

for(i in 1:length(mainTreesTrees)){
  currentTree = mainTreesTrees[[i]]
  currentTips = currentTree$tip.label
  
  if(i%%100 == 0){
    message(i)
  }
  
  tipsToAdd = currentTips[!currentTips %in% alltips]
  
  alltips = append(alltips, tipsToAdd)
}

alltipsUnique = unique(alltips)

allTips = character()
for(i in 1:length(mainTreesTrees)){
  currentTips = mainTreesTrees[[i]]$tip.label
  
  tipsToAdd = currentTips[!currentTips %in% allTips]
  
  allTips = append(allTips, tipsToAdd)
}
length(allTips)

mainTreesTrees[[1]]

treesTips = 


newickSet = allZoonomiaAlignments[,2]
newickSet[1]
read.newick(text = newickSet[1])$tip.label


allTipsNewick = character()
for(i in 1:length(newickSet)){
  currentTips = read.newick(text = newickSet[i])$tip.label
  
  tipsToAdd = currentTips[!currentTips %in% allTipsNewick]
  
  allTipsNewick = append(allTipsNewick, tipsToAdd)
  
  if(i%%100 == 0){
    message(i)
  }
}
saveRDS(allTipsNewick, "Results/allZoonomiaAlignmentTips.rds")

length(allTipsNewick)

tipsOnlyInMasterTree = togaTree$tip.label[!togaTree$tip.label %in% allTipsNewick]

togaTreeTrimmed = drop.tip(togaTree, tipsOnlyInMasterTree)

plot.phylo(togaTreeTrimmed)


?readTrees









demoTree = read.tree("data/nickDemoTree.nwk")



mainTrees = list()

mainTrees$masterTree = demoTree

saveRDS(mainTrees, "Data/batDemoMaintrees.rds")


ageTree = readRDS("Data/Output/Old")

names(enrichment1[1])
enrichment1[1]

library(xlsx)
write.xlsx(enrichment1[1], file="Output/MaturityLIfespanPercent/GOResults.xlsx", sheetName=names(enrichment1[1]), row.names=FALSE)
write.xlsx(enrichment2[1], file="Output/MaturityLIfespanPercent/GOResults.xlsx", sheetName=names(enrichment2[1]), append = T, row.names=FALSE)
write.xlsx(enrichment3[1], file="Output/MaturityLIfespanPercent/GOResults.xlsx", sheetName=names(enrichment3[1]), append = T, row.names=FALSE)
write.xlsx(enrichment4[1], file="Output/MaturityLIfespanPercent/GOResults.xlsx", sheetName=names(enrichment4[1]), append = T, row.names=FALSE)
write.xlsx(enrichment5[1], file="Output/MaturityLIfespanPercent/GOResults.xlsx", sheetName=names(enrichment5[1]), append = T, row.names=FALSE)



?readTrees

<<<<<<< HEAD
=======
togaTree = read.tree("Data/togaTree.nwk")
plot.phylo(togaTree)

mainTrees = readTrees(mainTreesLocation, masterTree = togaTree)

>>>>>>> 708e40e56dcfa8d258cfdd49b8ee5d6ac6f3167b


which(!relevantSpeciesNames %in% mainTrees$masterTree$tip.label)

f2tInputList

foreground2Tree()


is.binary(mainTrees$masterTree)


mainTrees$masterTree$tip.label %in% "vs_HLloxAfr4"

char2TreeCategorical

tipvals = commonPhenotypeVector; treesObj = commonMainTrees; useSpecies = commonSpeciesFilter; model = modelType; anctrait = ancestralTrait; plot = T
tipvals = phenotypeVector; treesObj = mainTrees; useSpecies = speciesFilter; model = modelType; anctrait = ancestralTrait; plot = T

scientificStates = intlabels$mapped_states
commonSates = intlabels$mapped_states


droppedtips = c("tupChi1","HLfukDam1","HLmarMar1","HLdipOrd2","HLbalMys1","mm10","HLenhLut1","HLdesRot1","echTel2","dasNov3")

droppedtips[7]

plotTree(mainTrees$masterTree)

ZonomNameConvertVectorCommon(droppedtips, tipColumn = "HillerName")
ZonomNameConvertVectorCommon
scientificTips = tipvals
commontips = tipvals
function (tipvals, treesObj, useSpecies = NULL, model = "ER", 
          root_prior = "auto", plot = FALSE, anctrait = NULL) 
{
  mastertree = treesObj$masterTree
  if (!is.null(useSpecies)) {
    sp.miss = setdiff(mastertree$tip.label, useSpecies)
    if (length(sp.miss) > 0) {
      message(paste0("Species from master tree not present in useSpecies: ", 
                     paste(sp.miss, collapse = ",")))
    }
    useSpecies = intersect(mastertree$tip.label, useSpecies)
    mastertree = pruneTree(mastertree, useSpecies)
    mastertree = unroot(mastertree)
  }
  else {
    mastertree = pruneTree(mastertree, intersect(mastertree$tip.label, 
                                                 names(tipvals)))
    mastertree = unroot(mastertree)
  }
  if (is.null(anctrait)) {
    tipvals <- tipvals[mastertree$tip.label]
    intlabels <- map_to_state_space(tipvals)
    print("The integer labels corresponding to each category are:")
    print(intlabels$name2index)
    ancliks = getAncLiks(mastertree, intlabels$mapped_states, 
                         rate_model = model, root_prior = root_prior)
    states = rep(0, nrow(ancliks))
    for (i in 1:length(states)) {
      states[i] = which.max(ancliks[i, ])
    }
    states = c(intlabels$mapped_states, states)
    tree = mastertree
    tree$edge.length = states[tree$edge[, 2]]
    if (length(unique(tipvals)) == 2) {
      if (sum(!unique(tipvals) %in% c(TRUE, FALSE)) > 0) {
        message("Returning categorical tree for binary phenotype because phenotype values are not TRUE/FALSE")
      }
      else {
        tree$edge.length = ifelse(tree$edge.length == 
                                    2, 1, 0)
        print("There are only 2 categories: returning a binary phenotype tree.")
        if (plot) {
          plotTree(tree)
        }
        return(tree)
      }
    }
    if (plot) {
      plotTreeCategorical(tree, category_names = intlabels$state_names, 
                          master = mastertree, node_states = states)
    }
    return(tree)
  }
  else {
    if (length(unique(tipvals)) <= 2) {
      fgspecs <- names(tipvals)[tipvals != anctrait]
      res <- foreground2Tree(fgspecs, treesObj, plotTree = plot, 
                             clade = "terminal", useSpecies = useSpecies)
      print("There are only 2 categories: returning a binary phenotype tree.")
      if (plot) {
        plotTree(res)
      }
      return(res)
    }
    else {
      tipvals <- tipvals[mastertree$tip.label]
      intlabels <- map_to_state_space(tipvals)
      j <- which(intlabels$state_names == anctrait)
      if (length(j) < 1) {
        warning("The ancestral trait provided must match one of the traits in the phenotype vector.")
      }
      res = mastertree
      res$edge.length <- rep(j, length(res$edge.length))
      traits <- intlabels$state_names
      for (trait in traits) {
        if (trait == anctrait) {
          next
        }
        i <- which(intlabels$state_names == trait)
        res$edge.length[nameEdges(res) %in% names(tipvals)[tipvals == 
                                                             trait]] = i
      }
      names(res$edge.length) = nameEdges(res)
      if (plot) {
        states = res$edge.length[order(res$edge[, 2])]
        states = c(j, states)
        plotTreeCategorical(res, category_names = traits, 
                            master = treesObj$masterTree, node_states = states)
      }
      print("Category names are mapped to integers as follows:")
      print(intlabels$name2index)
      return(res)
    }
  }
}








if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("qvalue")

RERObject[which(!is.na(RERObject))]

sapply(RERObject, is.numeric)

str(RERObject)

is.numeric(RERObject)

testMaster = readNexus("Data/output.nex")

testTree = read.tree("Data/nickDemoTree.nwk")
plotTree(testTree)

masterTree = test$masterTree

testMaster

masterTree$tip.label

df <- data.frame(
  A = c(1, 2, NA, 4),
  B = c(NA, 5, 6, NA),
  C = c(7, 8, 9, 10)
)

nonNA <- rowSums(!is.na(RERObject))

nonNA = nonNA[order(nonNA, decreasing = T)]

which(rownames(RERObject) == "UTP20")

which(!is.na(RERObject[1984,]))

trimRERs = RERObject[,which(!is.na(RERObject[1984,]))]

View(trimRERs)


nonNaVals = colSums(!is.na(trimRERs))

hist(nonNaVals)

manyNaSpecies = nonNaVals[order(scale(nonNaVals))]
manyNaSpecies = manyNaSpecies[-which(is.na(names(manyNaSpecies)))]
manyNaSpeciesCommon = manyNaSpecies

View(permulationData)

names(manyNaSpeciesCommon) = ZonomNameConvertVectorCommon(names(manyNaSpecies), tipColumn = "Zoonomia")

categoricalPermulations

getPermsBinaryFudged

foreground2Tree


function (foreground, treesObj, plotTree = T, clade = c("ancestral", 
                                                        "terminal", "all"), weighted = F, transition = "unidirectional", 
          useSpecies = NULL) 
{
  clade <- match.arg(clade)
  res = treesObj$masterTree
  if (!is.null(useSpecies)) {
    sp.miss = setdiff(res$tip.label, useSpecies)
    if (length(sp.miss) > 0) {
      message(paste0("Species from master tree not present in useSpecies: ", 
                     paste(sp.miss, collapse = ",")))
    }
    useSpecies = intersect(useSpecies, res$tip.label)
    res = pruneTree(res, useSpecies)
  }
  else {
    useSpecies = res$tip.label
  }
  foreground = intersect(foreground, useSpecies)
  res$edge.length <- rep(0, length(res$edge.length))
  if (clade == "terminal") {
    res$edge.length[nameEdges(res) %in% foreground] = 1
    names(res$edge.length) = nameEdges(res)
  }
  else if (clade == "ancestral") {
    weighted = F
    if (transition == "bidirectional") {
      res <- inferBidirectionalForegroundClades(res, foreground, 
                                                ancestralOnly = T)
    }
    else {
      res <- inferUnidirectionalForegroundClades(res, foreground, 
                                                 ancestralOnly = T)
    }
  }
  else {
    if (transition == "bidirectional") {
      res <- inferBidirectionalForegroundClades(res, foreground, 
                                                ancestralOnly = F)
    }
    else {
      res <- inferUnidirectionalForegroundClades(res, foreground, 
                                                 ancestralOnly = F)
    }
  }
  if (weighted) {
    if (clade == "all") {
      tobeweighted <- rep(TRUE, length(res$edge.length))
      tobeweighted[res$edge.length == 0] <- FALSE
      while (sum(tobeweighted) > 0) {
        edgetodo <- which(tobeweighted == T)[1]
        clade.down.edges = getAllCladeEdges(res, edgetodo)
        if (length(clade.down.edges) > 1) {
          clade.edges = c(clade.down.edges, edgetodo)
          clade.edges.toweight <- clade.edges[res$edge.length[clade.edges] == 
                                                1]
          res$edge.length[clade.edges.toweight] <- 1/(length(clade.edges.toweight))
          tobeweighted[clade.edges] <- FALSE
        }
        else {
          tobeweighted[clade.down.edges] <- FALSE
        }
      }
    }
    else if (clade == "terminal") {
      tobeweightededgeterminalnode <- unique(res$edge[(res$edge[, 
                                                                2] %in% c(1:length(res$tip.label))), 1])
      tobeweighted <- setdiff(match(tobeweightededgeterminalnode, 
                                    res$edge[, 2]), NA)
      for (edgetodo in tobeweighted) {
        clade.down.edges = getAllCladeEdges(res, edgetodo)
        if (all(res$edge.length[clade.down.edges] == 
                1)) {
          res$edge.length[clade.down.edges] <- 0.5
        }
      }
    }
  }
  if (plotTree) {
    res2 = res
    mm = min(res2$edge.length[res2$edge.length > 0])
    res2$edge.length[res2$edge.length == 0] = max(0.02, mm/20)
    plot(res2, main = paste0("Clade: ", clade, "\nTransition: ", 
                             transition, "\nWeighted: ", weighted), cex = 0.5)
    if (weighted) {
      labs <- round(res$edge.length, 3)
      labs[labs == 0] <- NA
      edgelabels(labs, col = "black", bg = "transparent", 
                 adj = c(0.5, -0.5), cex = 0.4, frame = "n")
    }
  }
  res
}


?plotTreeHighlightBranches

hist(nonNA)

arf5 = RERObject[1,]
length(which(!is.na(arf5)))

tipCol = "Zoonomia"

length(which(is.na(RERObject))) / length(RERObject) 

workingTree

masterTree = workingTree

View(mainTrees$report)

colSums(mainTrees$report)

hist(colSums(mainTrees$report))

zScores = scale(colSums(mainTrees$report))[,1]


which(zScores < -3)

mean(scale(colSums(mainTrees$report)))

which(colSums(mainTrees$report) < 12000)

ZonomNameConvertVectorCommon(names(which(colSums(mainTrees$report) < 12000)), tipColumn = )

dev.off()
dev.new()
mergedData = read.csv("Data/mergedData.csv")

which(!is.na(mergedData$HillerName) & !is.na(mergedData$Zoonomia))
mergedData$HillerZoonomiaOverlap = 0
mergedData$HillerZoonomiaOverlap[which(!is.na(mergedData$HillerName) & !is.na(mergedData$Zoonomia))] = 1

mergedData$HillerZoonomiaOverlapLogical = F
mergedData$HillerZoonomiaOverlapLogical[which(!is.na(mergedData$HillerName) & !is.na(mergedData$Zoonomia))] = T


#write.csv(mergedData, "Data/mergedData.csv", row.names = F)
#

binaryTree = readTest 
originalMasterTree = mainTrees$masterTree
masterTree = ZoonomTreeNameToCommon(masterTree, tipCol = "Zoonomia")

namesToKeep = binaryTree$tip.label
prunedMaster = drop.tip(masterTree, which(!masterTree$tip.label %in% namesToKeep))

masterTree= prunedMaster
masterTree = ZoonomTreeNameToCommon(masterTree, tipCol = "Zoonomia")





hist(prunedMaster$edge.length, breaks=100)

?hist

# -- Autopruner --
returnEdgeTable = T
procedurePlot = T
tipsToKeep = "Domestic goat"

overlapData = mergedData[mergedData$HillerZoonomiaOverlap ==1, ]

overlapTips = overlapData$Zoonomia














maturityCorrelations = readRDS("Output/MaturityLifespanPercent/MaturityLifespanPercentCorrelationFile.rds")

orderedMaturityCors = maturityCorrelations[order(maturityCorrelations$P),]

plotTree(binaryForegroundTreeOutput)

binaryForegroundTreeOutput


manualAnnots = read.csv("Data/mergedData.csv")

which(manualAnnots$Meyer.Lab.Classification)

manualAnnots$HerbBinary = rep(0)
manualAnnots$HerbBinary[grep("Herb", manualAnnots$Meyer.Lab.Classification)] = 1

manualAnnots$CarnBinary = rep(0)
manualAnnots$CarnBinary[grep("Carn", manualAnnots$Meyer.Lab.Classification)] = 1
manualAnnots$CarnBinary[grep("Pisc", manualAnnots$Meyer.Lab.Classification)] = 1


manualAnnots$CarnInsectBinary = rep(0)
manualAnnots$CarnInsectBinary[grep("Carn", manualAnnots$Meyer.Lab.Classification)] = 1
manualAnnots$CarnInsectBinary[grep("Pisc", manualAnnots$Meyer.Lab.Classification)] = 1
manualAnnots$CarnInsectBinary[grep("Insect", manualAnnots$Meyer.Lab.Classification)] = 1

#write.csv(manualAnnots, "Data/mergedData.csv", row.names = F)


testTree = mainTrees$masterTree
testTree = binaryForegroundTreeOutput


tipToDrop = testTree$tip.label[which(!testTree$tip.label %in% relevantSpeciesNames)]

trimTree = drop.tip(testTree, tipToDrop)
trimTree = drop.tip(trimTree, "mm10")

plotTree(trimTree)

binaryForegroundTreeOutput = trimTree

unique(orderedMaturityCors$p.adj)

source("Src/Reu/RERConvergeFunctions.r")

ZoonomTreeNameToCommon(binaryForegroundTreeOutput, tipCol = "Zoonomia")

tree = binaryForegroundTreeOutput
namesVector = tipNames

permdata = readRDS("Output/CategoricalDiet3Phen/CategoricalDiet3PhenPermulationsData1.rds")

permdtaaTrim = permdata

permdtaaTrim[[2]] = permdtaaTrim[[2]][1:3]
permdtaaTrim[[3]] = permdtaaTrim[[3]][1:3]

permDataTrim2 = permdata
permDataTrim2[[2]] = permDataTrim2[[2]][4:6]
permDataTrim2[[3]] = permDataTrim2[[3]][4:6]

permDataTrim3 = permdata
permDataTrim3[[2]] = permDataTrim3[[2]][7:9]
permDataTrim3[[3]] = permDataTrim3[[3]][7:9]

permDataTrim4 = permdata
permDataTrim4[[2]] = permDataTrim4[[2]][10:12]
permDataTrim4[[3]] = permDataTrim4[[3]][10:12]

saveRDS(permdtaaTrim, "Output/CategoricalDiet3Phen/CategoricalDiet3PhenPermulationsDataDev.rds")
rm(permdata)
rm(permdtaaTrim)

library(data.table)
source("Src/Reu/CategoricalPermulationsParallelFunctions.R")

permCorrelations = getPermPvalsCategorical3(correlationsObject, permulationsData$trees, phenotypeVector, mainTrees, RERObject)
permCorrelations = getPermPvalsCategorical2(correlationsObject, permulationsData$trees, phenotypeVector, mainTrees, RERObject)

permulationsDataSaving = permulationsData

permulationsData = permulationsDataSaving
permCorrelationsNew = CategoricalPermulationGetCor(correlationsObject, permulationsData$trees, phenotypeVector, mainTrees, RERObject, report=T)


permCorrelationsNew2 = CategoricalPermulationGetCor(correlationsObject, permDataTrim2$trees, phenotypeVector, mainTrees, RERObject, report=T)
permCorrelationsCombine = permCorrelationsNew


intermediate1 = permCorrelationsNew
intermediate2 = permCorrelationsNew2
i=1

combineCategoricalPermulationIntermediates = function(intermediate1, intermediate2){
  intermediatesCombined = intermediate1
  intermediatesCombined[[1]] = cbind(intermediatesCombined[[1]], intermediate2[[1]])
  for(i in 1:length(intermediatesCombined[[2]])){
    intermediatesCombined[[2]][[i]] = cbind(intermediatesCombined[[2]][[i]], intermediate2[[2]][[i]])
  }
  intermediatesCombined[[3]] = cbind(intermediatesCombined[[3]], intermediate2[[3]])
  for(i in 1:length(intermediatesCombined[[4]])){
    intermediatesCombined[[4]][[i]] = cbind(intermediatesCombined[[4]][[i]], intermediate2[[4]][[i]])
  }
  return(intermediatesCombined)
}

paste(basePermulationsFilename, 1, ".rds", sep= "")

permCorrelationsNew3 = CategoricalPermulationGetCor(correlationsObject, permDataTrim3$trees, phenotypeVector, mainTrees, RERObject, report=T)
permCorrelationsNew4 = CategoricalPermulationGetCor(correlationsObject, permDataTrim4$trees, phenotypeVector, mainTrees, RERObject, report=T)


saveRDS(permCorrelationsNew, paste(basePermulationsFilename, 1, ".rds", sep= ""))
saveRDS(permCorrelationsNew2, paste(basePermulationsFilename, 2, ".rds", sep= ""))
saveRDS(permCorrelationsNew3, paste(basePermulationsFilename, 3, ".rds", sep= ""))
saveRDS(permCorrelationsNew4, paste(basePermulationsFilename, 4, ".rds", sep= ""))


combinedPerms = combineCategoricalPermulationIntermediates(intermediate1, intermediate2)



ist1 = permCorrelationsNew
list2 = permCorrelationsNew2
listIndex = 1
matrix1 = list1[[listIndex]]
matrix2 = list2[[listIndex]]

matrixCombine = cbind(matrix1, matrix2)
matrix1Length = ncol(matrix1)
matrix2Length = ncol(matrix2)

matrixCombine[,4] = matrix2[,1]


permCorrelationsCombine[[1]][,matrix1Length+1:matrix1Length+matrix2Length] = permCorrelationsNew2[[1]][,1:matrix2Length]

permCorrelationsCombine[[1]][,matrix1Length+1:matrix1Length+matrix2Length]


time1 = Sys.time()
time2 = Sys.time()
timeDif = time2 - time1
str(timeDif)
attr(timeDif, "units")


permPval = CategoricalCalculatePermulationPValues(correlationsObject, permCorrelationsNew, report=F)
permPval = CategoricalCalculatePermulationPValues(correlationsObject, combinedPerms, report=F)



?correlateWithCategoricalPhenotype

all.equal(paths, pathsObject)
modelType ="ER"

#categoricalCorrelationOld = categoricalCorrelation

all.equal(categoricalCorrelation[[2]][1], categoricalCorrelationOld[[2]][1])


pairwiseTableNames = names(correlationsObject[[2]])                               #Prepare to repalce the number-number titles with phenotype-phenotype titles
for(i in 1:length(categoryNames)){                                            #for each phenotype
  pairwiseTableNames= gsub(i, names(categoryNames)[i], pairwiseTableNames)                        #replace the number with the phenotype name  
}
names(correlationsObject[[2]]) = pairwiseTableNames                               #update the dataframe titles

?install.packages()


source("Src/Reu/ZoonomTreeNameToCommon.R")
source("Src/Reu/ZonomNameConvertVector.R")
masterTree = mainTrees$masterTree

pdf(height = 25)
ZoonomTreeNameToCommon(readTest, isForegroundTree = T)


dev.off()

mainTips = masterTree$tip.label
commonTips = ZonomNameConvertVectorCommon(mainTips)
scietificTips = ZonomNameConvertVectorCommon(mainTips, common = F)
bothTips = append(scietificTips, commonTips)
write.csv(bothTips, file="Results/SpeciesNames.csv")


CVHRERs = readRDS("Output/CVHRemake/CVHRemakeRERFile.rds")
CVHCorrelation = readRDS("Output/CVHRemake/CVHRemakeCorrelationFile.rds")
CVHCorrelation[order(abs(CVHCorrelation$Rho), decreasing = F),]
source("Src/Reu/RERViolinPlot.R")
mainTrees = readRDS("Data/RemadeTreesAllZoonomiaSpecies.rds")
phenotypeTree = readRDS("Output/CVHRemake/CVHRemakeBinaryForegroundTRee.rds")
fgSpecies = readRDS("Output/CVHRemake/CVHRemakeBinaryTreeForegroundSpecies.rds")
rerViolinPlot()
quickViolin = function(gene){
  rerViolinPlot(mainTrees = mainTrees, CVHRERs, phenotypeTree = phenotypeTree, foregroundSpecies = foregroundSpecies, geneOfInterest = gene, foregroundName = "Carnivore", backgroundName = "Herbivore", backgroundColor = "darkgreen")
}
quickViolin("SLC47A1")
quickViolin("SPEGNB")
quickViolin("SPECC1L")
quickViolin("CRB1")
quickViolin("SDS")

quickViolin("PROM2")


source("Src/Reu/ZoonomTreeNameToCommon.R")
pdf(height = 25)
ZoonomTreeNameToCommon(binaryForegroundTreeOutput)
dev.off()

source("Src/Reu/RERConvergeFunctions.R")

realCors = correlationsObject
intermediateList = permCorrelations
gene = 1
testOut = list(res = realCors, pvals = list(corsMatPvals, Ppvals), effsize = list(corsMatEffSize, Peffsize))
testOutTrim = testOut$res


categoricalPerms = readRDS("Output/CategoricalDiet3Phen/CategoricalDiet3PhenPermulationsPValueCorrelations.rds")
categoricalNoPerms = readRDS("Output/CategoricalDiet3Phen/CategoricalDiet3PhenCombinedCategoricalCorrelationFile.rds")
categoricalPairwise = readRDS("Output/CategoricalDiet3Phen/CategoricalDiet3PhenPairwiseCorrelationFile.rds")


                            #select the group of pairwise comparisons

phenotypeVectorFilename = paste(outputFolderName, filePrefix, "CategoricalPhenotypeVector.rds",sep="") #select the phenotype vector based on prefix
phenotypeVector = readRDS(phenotypeVectorFilename)                            #load in the phenotype vector 
categories = map_to_state_space(phenotypeVector)                              #and use it to connect branch lengths to phenotype name
categoryNames = categories$name2index                                         #store the length-phenotype connection

pairwiseTableNames = names(categoricalPerms[[2]])                               #Prepare to replace the number-number titles with phenotype-phenotype titles
for(i in 1:length(categoryNames)){                                            #for each phenotype
  pairwiseTableNames= gsub(i, names(categoryNames)[i], pairwiseTableNames)    #replace the number with the phenotype name  
}
names(categoricalPerms[[2]]) = pairwiseTableNames                               #update the dataframe titles

permulationsPValuesFilename = paste(outputFolderName, filePrefix, "PermulationsPValueCorrelations.rds", sep= "")
saveRDS(categoricalPerms, permulationsPValuesFilename)

permulationsPValuesOutput = categoricalPerms

permulationsPValuesOverallFilename = paste(outputFolderName, filePrefix, "PermulationsOverallCorrelations.rds", sep= "")
saveRDS(permulationsPValuesOutput[[1]], permulationsPValuesOverallFilename)

for(i in 1:length(pairwiseTableNames)){
  pairwiseTableNames= gsub(" ", "", pairwiseTableNames)
  permulationsPValuesPairFilename = paste(outputFolderName, filePrefix, "PermulationsCorrelations", pairwiseTableNames[i],".rds", sep= "")
  saveRDS(permulationsPValuesOutput[[2]][i], permulationsPValuesPairFilename)
}

#---

outputSubdirectoryNoslash = paste(outputFolderName, "Overall", sep = "")
if(!dir.exists(outputSubdirectoryNoslash)){                       #create that directory if it does not exist
  dir.create(outputSubdirectoryNoslash)
}
outputSubdirectory = paste(outputSubdirectoryNoslash, "/", sep="")

permulationsPValuesOverallFilename = paste(outputSubdirectory, filePrefix, "PermulationsOverallCorrelations.rds", sep= "")
saveRDS(permulationsPValuesOutput[[1]], permulationsPValuesOverallFilename)

for(i in 1:length(pairwiseTableNames)){
  pairwiseTableNames= gsub(" ", "", pairwiseTableNames)
  
  outputSubdirectoryNoslash = paste(outputFolderName, pairwiseTableNames[i], sep = "")
  if(!dir.exists(outputSubdirectoryNoslash)){                       #create that directory if it does not exist
    dir.create(outputSubdirectoryNoslash)
  }
  outputSubdirectory = paste(outputSubdirectoryNoslash, "/", sep="")
  
  permulationsPValuesPairFilename = paste(outputSubdirectory, filePrefix, "PermulationsCorrelations", pairwiseTableNames[i],".rds", sep= "")
  saveRDS(permulationsPValuesOutput[[2]][i], permulationsPValuesPairFilename)
}


nonPermCorrelations = readRDS("Output/CategoricalDiet3Phen/CategoricalDiet3PhenCombinedCategoricalCorrelationFile.rds")


nonPermCorrelations[[1]]$Rho

permulationsPValuesOutput=categoricalPerms

i=1

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

test1 = permulationsPValuesOutput[[2]][[i]]


correlData$p.adj  

nonPermCorrelations[[2]]$`1 - 3`$p.adj
#saveRDS(mainTrees, "Data/FirstExpressionTrees.rds")
write.csv(mainTrees$masterTree$tip.label, "test.csv")

which(!is.na(RERObject[,1]))

mainTrees$masterTree$tip.label %in% speciesFilter


RERObjectNew = getAllResiduals(mainTrees, useSpecies = speciesFilter, plot = F, min.sp = 3, transform = "none", maxT = 5)

all(is.na(RERObjectNew))



?getAllResiduals()

function (treesObj, cutoff = NULL, transform = "sqrt", weighted = T, 
          useSpecies = NULL, min.sp = 10, scale = T, doOnly = NULL, 
          maxT = NULL, scaleForPproj = F, mean.trim = 0.05, plot = T) 
{
  if (is.null(cutoff)) {
    cutoff = quantile(treesObj$paths, 0.05, na.rm = T)
    message(paste("cutoff is set to", cutoff))
  }
  if (weighted) {
    weights = computeWeightsAllVar(treesObj$paths, transform = transform, 
                                   plot = plot)
    residfunc = fastLmResidMatWeighted
  }
  else {
    residfunc = fastLmResidMat
  }
  if (is.null(useSpecies)) {
    useSpecies = treesObj$masterTree$tip.label
  }
  if (is.null(maxT)) {
    maxT = treesObj$numTrees
  }
  if (transform != "none") {
    transform = match.arg(transform, c("sqrt", "log"))
    transform = get(transform)
  }
  else {
    transform = NULL
  }
  cm = intersect(treesObj$masterTree$tip.label, useSpecies)
  sp.miss = setdiff(treesObj$masterTree$tip.label, useSpecies)
  if (length(sp.miss) > 0) {
    message(paste0("Species from master tree not present in useSpecies: ", 
                   paste(sp.miss, collapse = ",")))
  }
  rr = matrix(nrow = nrow(treesObj$paths), ncol = ncol(treesObj$paths))
  maxn = rowSums(treesObj$report[, cm])
  if (is.null(doOnly)) {
    doOnly = 1
  }
  else {
    maxT = 1
  }
  skipped = double(nrow(rr))
  skipped[] = 0
  for (i in doOnly:(doOnly + maxT - 1)) {
    if (sum(!is.na(rr[i, ])) == 0 && !skipped[i] == 1) {
      tree1 = treesObj$trees[[i]]
      both = intersect(tree1$tip.label, cm)
      if (length(both) < min.sp) {
        next
      }
      tree1 = unroot(pruneTree(tree1, both))
      allreport = treesObj$report[, both]
      ss = rowSums(allreport)
      iiboth = which(ss == length(both))
      if (length(iiboth) < 2) {
        message(paste("Skipping i =", i, "(no other genes with same species set)"))
        next
      }
      nb = length(both)
      ai = which(maxn[iiboth] == nb)
      message(paste("i=", i))
      if (T) {
        ee = edgeIndexRelativeMaster(tree1, treesObj$masterTree)
        ii = treesObj$matIndex[ee[, c(2, 1)]]
        allbranch = treesObj$paths[iiboth, ii]
        if (is.null(dim(allbranch))) {
          message(paste("Issue with gettiing paths for genes with same species as tree", 
                        i))
          return(list(iiboth = iiboth, ii = ii))
        }
        if (weighted) {
          allbranchw = weights[iiboth, ii]
        }
        if (scaleForPproj) {
          nv = apply(scaleMatMean(allbranch), 2, mean, 
                     na.rm = T, trim = mean.trim)
        }
        else {
          nv = apply(allbranch, 2, mean, na.rm = T, trim = mean.trim)
        }
        iibad = which(allbranch < cutoff)
        if (!is.null(transform)) {
          nv = transform(nv)
          allbranch = transform(allbranch)
        }
        allbranch[iibad] = NA
        if (!scale) {
          if (!weighted) {
            proj = residfunc(allbranch[ai, , drop = F], 
                             model.matrix(~1 + nv))
          }
          else {
            proj = residfunc(allbranch[ai, , drop = F], 
                             model.matrix(~1 + nv), allbranchw[ai, , 
                                                               drop = F])
          }
        }
        else {
          if (!weighted) {
            proj = residfunc(allbranch[, , drop = F], 
                             model.matrix(~1 + nv))
          }
          else {
            proj = residfunc(allbranch[, , drop = F], 
                             model.matrix(~1 + nv), allbranchw)
          }
          proj = scale(proj, center = F)[ai, , drop = F]
        }
        rr[iiboth[ai], ii] = proj
      }
    }
  }
  message("Naming rows and columns of RER matrix")
  rownames(rr) = names(treesObj$trees)
  colnames(rr) = namePathsWSpecies(treesObj$masterTree)
  rr
}

correlData[order(correlData$p.adj),]

upham = ReadTntTree("Data/PrunedUphamTree65LiverTPMSpecies.tre")
??Tnt

library(TreeTools)

upham = read.tree("Data/PrunedUphamTree65LiverTPMSpecies.tre")

masterTree = mainTrees$masterTree
plotTree(upham)
plotTree(masterTree)

unrootupham = UnrootTree(upham)

#saveRDS(mainTrees, "Data/NoSignFirstExpressionTrees.rds")

mainTrees$masterTree = unrootupham
#saveRDS(mainTrees, "Data/NoSignFirstExpressionTreesNewMaster.rds")

oldRERs = readRDS("Output/CVHRemake/CVHRemakeRERFile.rds")


cat4Crrels = readRDS("Output/CategoricalDiet4Phen/CategoricalDiet4PhenPairwiseCorrelationFile.rds")
pairwiseTableNames = names(cat4Crrels)

cat5Crrels = readRDS("Output/CategoricalDiet5Phen/CategoricalDiet5PhenPairwiseCorrelationFile.rds")
pairwiseTableNames = names(cat5Crrels)

cat4Crrels[[1]]

for(i in 1:length(pairwiseTableNames)){
  pairwiseTableNames= gsub(" ", "", pairwiseTableNames)
  
  outputSubdirectoryNoslash = paste(outputFolderName, pairwiseTableNames[i], sep = "")
  if(!dir.exists(outputSubdirectoryNoslash)){                       #create that directory if it does not exist
    dir.create(outputSubdirectoryNoslash)
  }
  outputSubdirectory = paste(outputSubdirectoryNoslash, "/", sep="")
  
  permulationsPValuesPairFilename = paste(outputSubdirectory, filePrefix, pairwiseTableNames[i], "CorrelationFile",".rds", sep= "")
  saveRDS(cat5Crrels[[i]], permulationsPValuesPairFilename)
}
j=1

correlation[order(correlation$p.adj),]

saveRDS(mainTrees, "Data/NoSignExpressionTreesRound3.rds")

signmulitmat = readRDS("Data/SignMultiplierExpressionTreesRound3.rds")

unsignedRER = RERObject
signedRER = unsignedRER * signmulitmat

RERObject = signedRER

all.equal(signedRER, unsignedRER)
all.equal(signedRER, RERObject)


# ----------# 


foreground = permulatedForeground; treesObj = mainTrees; plotTree=F; clade="all"; transition="bidirectional"; useSpecies=speciesFilter; weighted = F; 

source("Src/Reu/RERCOnvergeFunctions.R")
foreground2TreeDebug = function (foreground, treesObj, plotTree = T, clade = c("ancestral", "terminal", "all"), weighted = F, transition = "unidirectional", useSpecies = NULL) 
{
  res = treesObj$masterTree
  if (!is.null(useSpecies)) {
    sp.miss = setdiff(res$tip.label, useSpecies)
    if (length(sp.miss) > 0) {
      message(paste0("Species from master tree not present in useSpecies: ", 
                     paste(sp.miss, collapse = ",")))
    }
    useSpecies = intersect(useSpecies, res$tip.label)
    res = pruneTree(res, useSpecies)
  }
  foreground = intersect(foreground, useSpecies)
  res$edge.length <- rep(0, length(res$edge.length))
  {
    if (transition == "bidirectional") {
      res <- inferBidirectionalForegroundCladesDebug(res, foreground, 
                                                ancestralOnly = F)
    }
  }
  res
}


treeinput = res; foreground = foreground; ancestralOnly = F; i=1

inferBidirectionalForegroundCladesDebug = function(treeinput, foreground = NULL, ancestralOnly = F){
  tree <- treeinput
  tip.vals=rep(0, length(tree$tip.label))
  names(tip.vals)=tree$tip.label
  tip.vals[foreground]=1
  tmp=cbind(as.character(tip.vals))
  rownames(tmp)=names(tip.vals)
  tip.vals=tmp
  #Add option to function for "type" within ancestral.pars
  ancres=ancestral.pars(tree, df<-as.phyDat(tip.vals, type="USER", levels=unique(as.character(tip.vals))),type="ACCTRAN" )
  ancres=unlist(lapply(ancres, function(x){x[2]}))
  internalVals=ancres
  #evals=matrix(nrow=nrow(treesObj$masterTree$edge), ncol=2)
  evals=matrix(nrow=nrow(tree$edge), ncol=2)
  eres=ancres
  #evals[,1]=eres[treesObj$masterTree$edge[,1]]
  evals[,1]=eres[tree$edge[,1]]
  #evals[,2]=eres[treesObj$masterTree$edge[,2]]
  evals[,2]=eres[tree$edge[,2]]
  tree$edge.length=evals[,2]-evals[,1]
  #res$edge.length[res$edge.length<1]=0
  if(!ancestralOnly){
    edgeIndex=which(tree$edge.length>0)
    edgeIndexNeg=which(tree$edge.length<0)
    edgeIndexAll = c(edgeIndex,edgeIndexNeg)
    edgeDirection = c(rep(1, length(edgeIndex)),rep(-1, length(edgeIndexNeg)))
    edgedf = data.frame(edgeIndexAll,edgeDirection)
    edgedf = edgedf[order(edgedf$edgeIndexAll),]
    clade.edges=NA
    clade.lengths=NA
    cladedf = data.frame(clade.edges,clade.lengths)
    for(i in 1:nrow(edgedf)) { #Does this go from ancestral to terminal?
      #save the clade until the edges no longer overlap
      clade.edges=getAllCladeEdgesDebug(tree, edgedf$edgeIndexAll[i])
      clade.edges=unique(c(edgedf$edgeIndexAll[i], clade.edges))
      if (any(clade.edges %in% cladedf$clade.edges)==F) {
        tree$edge.length[cladedf$clade.edges[which(cladedf$clade.lengths==1)]]=1
        if (edgedf$edgeDirection[i] == 1) {
          clade.lengths = c(rep(1,length(clade.edges)))
        } else {
          clade.lengths = c(rep(0,length(clade.edges)))
        }
        cladedf = data.frame(clade.edges,clade.lengths)
      } else {
        #update df lengths
        if (edgedf$edgeDirection[i] == 1) {
          cladedf$clade.lengths[which(cladedf$clade.edges %in% clade.edges)] = 1
        } else {
          cladedf$clade.lengths[which(cladedf$clade.edges %in% clade.edges)] = 0
        }
      }
    }
    #update edge lengths from the final clade
    tree$edge.length[cladedf$clade.edges[which(cladedf$clade.lengths==1)]]=1
    tree$edge.length[tree$edge.length<0]=0
  }
  tree$edge.length[tree$edge.length<0]=0
  tree
}

tree= tree; AncEdge = edgedf$edgeIndexAll[i]

getAllCladeEdgesDebug = function(tree, AncEdge){
  node=tree$edge[AncEdge,2]
  #get descendants
  iid=getDescendantsDebug(tree, node)
  #find their edges
  iim=match(iid, tree$edge[,2])
  iim
}

tree= tree; node = node; curr = NULL

getDescendantsDebug = function (tree, node, curr = NULL) 
{
  
  daughters <- tree$edge[which(tree$edge[, 1] == node), 2]
  curr <- daughters
  if(is.na(length(curr)) || is.null(length(curr)) || is.na(Ntip(tree)) || is.null(Ntip(tree)) || is.na(node) ||is.null(node)){
    print(node)
    print(curr)
    print(tree)
  }
  if (length(curr) == 0 && node <= Ntip(tree)) 
    curr <- node
  w <- which(daughters > Ntip(tree))
  if (length(w) > 0) 
    for (i in 1:length(w)) curr <- getDescendants(tree, daughters[w[i]], 
                                                  curr)
  return(curr)
}


readRDS("Output/LiverExpression3/LiverExpression3CombinedPrunedFastAppendedPermulationsPValue.rds")

permVal = readRDS("Output/LiverExpression3/LiverExpression3CombinedPrunedFastAppendedPermulationsPValue.rds")
permValNew = readRDS("Output/LiverExpression3/LiverExpression3CombinedPrunedFastAppendedPermulationsPValue.rds")

which(is.null(names(permVal)))
tail(names(permVal))
cleanedPermVal = a
  
which(!names(permVal) %in% unique(names(permVal)))

rownames(correlData)
which(!names(permVal) %in% rownames(correlData))

length(rownames(correlData))
length(names(permVal))

which(!rownames(correlData) %in% names(permVal))

cleanedPermVals = unique(permVal)

length(unique(names(permVal)))

which(duplicated(names(permVal)))
length(which(duplicated(names(permVal))))
4250-3401

permSub1 = permVal[1:850]
permSub2 = permVal[3401:4250]

names(permSub2)

all(names(permSub1) %in% names(permSub2))

all.equal(permSub1, permSub2)

permValClean = permVal[-c(3401:4250)]

all.equal(rownames(correlData), names(permValClean))
match(rownames(correlData), names(permValClean), nomatch = T)

names = data.frame(rownames(correlData)[c(1:17000)], names(permValNew))

unMatched = subset(names, names[[1]] != names[[2]])

all.equal(rownames(correlData), names(permValNew))

tail(rownames(unMatched))


pValset1= readRDS("Output/LiverExpression3/LiverExpression3CombinedPrunedFast1-3400PermulationsPValue.rds")
pValset2= readRDS("Output/LiverExpression3/LiverExpression3CombinedPrunedFast3401-6800PermulationsPValue.rds")
pValset3= readRDS("Output/LiverExpression3/LiverExpression3CombinedPrunedFast6801-10200PermulationsPValue.rds")
pValset4= readRDS("Output/LiverExpression3/LiverExpression3CombinedPrunedFast10201-13600PermulationsPValue.rds")
pValset5= readRDS("Output/LiverExpression3/LiverExpression3CombinedPrunedFast13601-17000PermulationsPValue.rds")


pvalCombined = append(pValset1, pValset2)
pvalCombined = append(pvalCombined, pValset3)
pvalCombined = append(pvalCombined, pValset4)
pvalCombined = append(pvalCombined, pValset5)

names = data.frame(rownames(correlData)[c(1:17000)], names(pvalCombined))
all.equal(rownames(correlData), names(pvalCombined))
unMatched = subset(names, names[[1]] != names[[2]])

saveRDS(pvalCombined, file = "Output/LiverExpression3/LiverExpression3CombinedPrunedFastAllPermulationsPValue.rds")
correlData = correlData[1:17000,]

#correlDataBackup = correlData
correlData = correlDataBackup

correlData = correlData[-which(correlData$permPValue ==0),]
correlDataPositive = correlDataPositive[which(correlDataPositive$permPValue==0)]


#saveRDS(correlData, "Output/LiverExpression3/LiverExpression3PermulatedCorrelations.rds")
#write.csv(correlData, "Output/LiverExpression3/LiverExpression3PermulatedCorrelations.CSV")

nameConvert = read.csv("Data/mart_export_convert_ensid_genename.csv")

oldNames = rownames(correlData)
convertTable = data.frame(oldNames)
convertTable$stoneIndex = match(oldNames, nameConvert$Gene.stable.ID)

convertTable$newName = NA
for(i in 1:length(convertTable$newName)){
  convertTable$newName[i] = nameConvert$Gene.name[convertTable$stoneIndex[i]]
}

correlDataFixedNames = correlData



length(which(duplicated(nameConvert$Gene.name)))

convertTable$newName[(which(duplicated(convertTable$newName)))]
convertTable$newNameFixed = convertTable$newName
convertTable$newNameFixed[(which(duplicated(convertTable$newNameFixed)))]

i=41
problemPositions = (which(duplicated(convertTable$newNameFixed)))
replacements = c("NOGENENAME1", "HLA-DQB1-DUPLICATE1", "HLA-DRB1-DUPLICATE2", "MISSINGCONVERSION", "NOGENENAME14", "HLA-DRB1-DUP1", "HLA-DRB1-DUP2", "NOGENENAME2", "TAP2-DUP1", "HLA-C-DUP1", "DDR1-DUP1", "ATP6V1G2-DUP1", "TNXB-DUP1", "EHMT2-DUP1", "NOGENENAME3", "GART-DUP1", "NOGENENAME4", "NOGENENAME5", "NOGENENAME6", "NOGENENAME7", "NOGENENAME8", "NOGENENAME9", "UPK3B-DUP1", "NOGENENAME10", "KIR3DL1-DUP1", "UNC79-DUP1","GTF2H2-DUP1", "OPRL1-DUP1", "PRODH-DUP1", "NOGENENAME11", "CCL3L3-DUP1", "MUC20-DUP1", "MUC4-DUP1", "NOGENENAME12", "NOGENENAME13", "MATR3-DUP1", "NXN-DUP1", "CABIN1-DUP1", "SEPTIN9-DUP1", "CCL4L2-DUP1", "NOGENENAME15")
for(i in 1:41){
  convertTable$newNameFixed[problemPositions[i]] = replacements[i] 
}
convertTable$newNameFixed[which(is.na(convertTable$newNameFixed))] = "MISSINGCONVERSION3"
rownames(correlDataFixedNames) = convertTable$newNameFixed

#
convertTable$newNameFixed[which(duplicated(convertTable$newNameFixed))]
length(replacements)
length(problemPositions)

problemNames = convertTable$newNameFixed[(which(duplicated(convertTable$newNameFixed)))]
convertTable$newNameFixed[problemPositions]

which(is.na(convertTable$newNameFixed))
#

all.equal(correlDataFixedNames, correlData)

#rownames(correlDataFixedNames)[13976] = "NOGENENAME16"

#saveRDS(correlDataFixedNames, "Output/LiverExpression3/LiverExpression3CorrelationDataPermulatedNamesConverted.rds")
correlData = correlDataFixedNames
i=4
correlationData = correlDataFixedNames



geneSetNames = unlist(annotationsList$GO_Biological_Process_2023$genesets[c(1:5407)])

length(which(convertTable$newNameFixed %in% geneSetNames))
length(convertTable$newNameFixed)
rerStats[which(rerStats == -Inf)] = 0

which(is.na(rerStats))
which(names(rerStats) == "")
which(rownames(correlationData) == "")
which(is.null(names(correlationData)))


vals=rerStats
gmt = annotationsList[[1]]

fastwilcoxGMT=function(vals, gmt, simple=T, use.all=F, num.g=10,genes=NULL, outputGeneVals=T, order=F,
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
alt = alternative

simpleAUCgenesRanks=function(pos, neg, alt = "two.sided"){
  
  posn=length(pos)
  negn=length(neg)
  posn=as.numeric(posn)
  negn=as.numeric(negn)
  stat=sum(pos)-posn*(posn+1)/2
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

which(names(stat) =="")
which(rownames(res) =="")
which(rownames(correlationData) =="")
res=correlationData
function(res){
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

convertTable[13976,]
nameConvert[124243,]
correlData = correlationData
#write.csv(correlDataFixedNames, "Output/LiverExpression3/LiverExpression3CorrelationDataPermulatedNamesConverted.csv")





##  ------ NAME CONVERSION ----- ##

nameConvert = read.csv("Data/mart_export_convert_ensid_genename.csv")

oldNames = rownames(correlData)
convertTable = data.frame(oldNames)
convertTable$stoneIndex = match(oldNames, nameConvert$Gene.stable.ID)

convertTable$newName = NA
for(i in 1:length(convertTable$newName)){
  convertTable$newName[i] = nameConvert$Gene.name[convertTable$stoneIndex[i]]
}

correlDataFixedNames = correlData

convertTable$newNameFixed = convertTable$newName
convertTable$newNameFixed[(which(duplicated(convertTable$newNameFixed)))]

problemPositions = (which(duplicated(convertTable$newNameFixed)))
replacements = c("NOGENENAME1", "HLA-DQB1-DUPLICATE1", "HLA-DRB1-DUPLICATE2", "MISSINGCONVERSION", "NOGENENAME14", "HLA-DRB1-DUP1", "HLA-DRB1-DUP2", "NOGENENAME2", "TAP2-DUP1", "HLA-C-DUP1", "DDR1-DUP1", "ATP6V1G2-DUP1", "TNXB-DUP1", "EHMT2-DUP1", "NOGENENAME3", "GART-DUP1", "NOGENENAME4", "NOGENENAME5", "NOGENENAME6", "NOGENENAME7", "NOGENENAME8", "NOGENENAME9", "UPK3B-DUP1", "NOGENENAME10", "KIR3DL1-DUP1", "UNC79-DUP1","GTF2H2-DUP1", "OPRL1-DUP1", "PRODH-DUP1", "NOGENENAME11", "CCL3L3-DUP1", "MUC20-DUP1", "MUC4-DUP1", "NOGENENAME12", "NOGENENAME13", "MATR3-DUP1", "NXN-DUP1", "CABIN1-DUP1", "SEPTIN9-DUP1", "CCL4L2-DUP1", "NOGENENAME15")

for(i in 1:41){
  convertTable$newNameFixed[problemPositions[i]] = replacements[i] 
}
convertTable$newNameFixed[which(is.na(convertTable$newNameFixed))] = "MISSINGCONVERSION3"
rownames(correlDataFixedNames) = convertTable$newNameFixed




# -- RER NAME CONVERSIONS --
RERsToTranslate = readRDS("Output/LiverExpression3/LiverExpression3RERFile.rds")

oldNames = rownames(RERsToTranslate)
convertTable = data.frame(oldNames)
convertTable$stoneIndex = match(oldNames, nameConvert$Gene.stable.ID)

convertTable$newName = NA
for(i in 1:length(convertTable$newName)){
  convertTable$newName[i] = nameConvert$Gene.name[convertTable$stoneIndex[i]]
}

RERsFixedNames = RERsToTranslate

convertTable$newNameFixed = convertTable$newName
convertTable$newNameFixed[(which(duplicated(convertTable$newNameFixed)))]


rownames(RERsFixedNames) = convertTable$newNameFixed

#saveRDS(RERsFixedNames, "Output/LiverExpression3/LiverExpression3RERsFileNamesFixed.rds")

saveRDS(correlData, file = "Output/CVHRemake/CVHRemakeCorrelationsFilePermulated.rds")

hist(correlData$P, breaks = 40)
?hist()


# ---------

realCors = correlationsObject; nullPhens =  permulationsData$trees; phenvals = phenotypeVector; treesObj = mainTrees; RERmat =  RERObject; method = "kw"; 
min.sp = 10; min.pos = 2; winsorizeRER = NULL; winsorizetrait = NULL; 
weighted = F; extantOnly = FALSE; report=T
saveRDS(nullPaths, "Results/nullPathDev.rds")
saveRDS(nullPhens, "Results/nullPhens.rds")
i =6

function (realCors, nullPhens, phenvals, treesObj, RERmat, method = "kw", 
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
  
  for (i in 7:17) {
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
  }
  output = list(corsMatEffSize, Peffsize, corsMatPvals, Ppvals)
  names(output) = c("corsMatEffSize", "Peffsize", "corsMatPvals", "Ppvals")
  return(output)
}


corrWithBinTest = function (RERmat, charP, min.sp = 10, min.pos = 2, weighted = "auto") 
{
  if (weighted == "auto") {
    if (any(charP > 0 & charP < 1, na.rm = TRUE)) {
      message("Fractional values detected, will use weighted correlation mode")
      weighted = T
    }
    else {
      weighted = F
    }
  }
  getAllCorTest(RERmat, charP, min.sp, min.pos, method = "k", weighted = weighted)
}

CVHRer = readRDS("Output/CVHRemake/CVHRemakeRERFile.rds")
CVHCharP = readRDS("Output/CVHRemake/CVHRemakePathsFile.rds")
length(unique(CVHCharP))

corrWithBinTest(CVHRer, CVHCharP, min.sp =10)



getAllCorTest = function (RERmat, charP, method = "auto", min.sp = 10, min.pos = 2, 
          winsorizeRER = NULL, winsorizetrait = NULL, weighted = F) 
{
  RERna = (apply(is.na(RERmat), 2, all))
  iicharPna = which(is.na(charP))
  if (!all(RERna[iicharPna])) {
    warning("Species in phenotype vector are a subset of the those used for RER computation. For best results run getAllResiduals with the useSpecies")
  }
  if (method == "auto") {
    lu = length(unique(charP))
    if (lu == 2) {
      method = "k"
      message("Setting method to Kendall")
    }
    else if (lu <= 5) {
      method = "s"
      message("Setting method to Spearman")
    }
    else {
      method = "p"
      message("Setting method to Pearson")
      if (is.null(winsorizeRER)) {
        message("Setting winsorizeRER=3")
        winsorizeRER = 3
      }
      if (is.null(winsorizetrait)) {
        message("Setting winsorizetrait=3")
        winsorizetrait = 3
      }
    }
  }
  win = function(x, w) {
    xs = sort(x[!is.na(x)], decreasing = T)
    xmax = xs[w]
    xmin = xs[length(xs) - w + 1]
    x[x > xmax] = xmax
    x[x < xmin] = xmin
    x
  }
  corout = matrix(nrow = nrow(RERmat), ncol = 3)
  rownames(corout) = rownames(RERmat)
  if (method == "aov" || method == "kw") {
    lu = length(unique(charP[!is.na(charP)]))
    n = choose(lu, 2)
    tables = lapply(1:n, matrix, data = NA, nrow = nrow(RERmat), 
                    ncol = 2, dimnames = list(rownames(RERmat), c("Rho", 
                                                                  "P")))
    if (method == "aov") {
      names(tables) = rep(NA, n)
    }
    else {
      for (i in 2:lu) {
        for (j in 1:(i - 1)) {
          index <- (i - 1) * (i - 2)/2 + j
          names(tables)[index] <- paste0(j, " - ", i)
        }
      }
    }
  }
  colnames(corout) = c("Rho", "N", "P")
  for (i in 1:nrow(corout)) {
    if (((nb <- sum(ii <- (!is.na(charP) & !is.na(RERmat[i, 
    ])))) >= min.sp)) {
      if (method == "kw" || method == "aov") {
        counts = table(charP[ii])
        num_groups = length(counts)
        if (num_groups < 2 || min(counts) < min.pos) {
          next
        }
      }
      else if (method != "p" && sum(charP[ii] != 0) < min.pos) {
        next
      }
      if (!weighted) {
        x = RERmat[i, ]
        indstouse = which(!is.na(x) & !is.na(charP))
        if (!is.null(winsorizeRER)) {
          x = win(x[indstouse], winsorizeRER)
        }
        else {
          x = x[indstouse]
        }
        if (!is.null(winsorizetrait)) {
          y = win(charP[indstouse], winsorizetrait)
        }
        else {
          y = charP[indstouse]
        }
        if (method == "aov") {
          yfacts = as.factor(y)
          df = data.frame(x, yfacts)
          colnames(df) = c("RER", "category")
          ares = aov(RER ~ category, data = df)
          sumsq = summary(ares)[[1]][1, 2]
          sumsqres = summary(ares)[[1]][2, 2]
          effect_size = sumsq/(sumsq + sumsqres)
          ares_pval = summary(ares)[[1]][1, 5]
          corout[i, 1:3] = c(effect_size, nb, ares_pval)
          tukey = TukeyHSD(ares)
          groups = rownames(tukey[[1]])
          unnamedinds = which(is.na(names(tables)))
          if (length(unnamedinds > 0)) {
            newnamesinds = which(is.na(match(groups, 
                                             names(tables))))
            if (length(newnamesinds) > 0) {
              names(tables)[unnamedinds][1:length(newnamesinds)] = groups[newnamesinds]
            }
          }
          for (k in 1:length(groups)) {
            name = groups[k]
            tables[[name]][i, "Rho"] = tukey[[1]][name, 
                                                  1]
            tables[[name]][i, "P"] = tukey[[1]][name, 
                                                4]
          }
        }
        else if (method == "kw") {
          yfacts = factor(y)
          kres = kwdunn.test(x, yfacts, ncategories = lu)
          effect_size = kres$kw$H/(nb - 1)
          corout[i, 1:3] = c(effect_size, nb, kres$kw$p)
          for (k in 1:length(kres$dunn$Z)) {
            tables[[k]][i, "Rho"] = kres$dunn$Z[k]
            tables[[k]][i, "P"] = kres$dunn$P.adjust[k]
          }
        }
        else {
          cres = cor.test(x, y, method = method, exact = F)
          corout[i, 1:3] = c(cres$estimate, nb, cres$p.value)
        }
      }
      else {
        charPb = (charP[ii] > 0) + 1 - 1
        weights = charP[ii]
        weights[weights == 0] = 1
        cres = wtd.cor(RERmat[i, ii], charPb, weight = weights, 
                       mean1 = F)
        corout[i, 1:3] = c(cres[1], nb, cres[4])
      }
    }
    else {
    }
  }
  corout = as.data.frame(corout)
  corout$p.adj = p.adjust(corout$P, method = "BH")
  if (method == "aov" || method == "kw") {
    for (i in 1:length(tables)) {
      tables[[i]] = as.data.frame(tables[[i]])
      tables[[i]]$p.adj = p.adjust(tables[[i]]$P, method = "BH")
    }
    return(list(corout, tables))
  }
  else {
    corout
  }
}

CVHRemakeBinTree = readRDS("OUtput/CVHRemake/CVHRemakeBinaryTree.rds")
CVHUpdateBinTree = readRDS("Output/CVHUpdate/CVHUpdateBinaryTree.rds")
all.equal(CVHRemakeBinTree, CVHUpdateBinTree)

all.equal(RERObject, CVHRer)
all.equal(pathsObject, CVHCharP)




function (realCors, nullPhens, phenvals, treesObj, RERmat, method = "kw", 
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
  }
  output = list(corsMatEffSize, Peffsize, corsMatPvals, Ppvals)
  names(output) = c("corsMatEffSize", "Peffsize", "corsMatPvals", "Ppvals")
  return(output)
}


> CategoricalCalculatePermulationPValues
function(realCors, intermediateList, start=1, end=NULL, report=F){
  {totalStart = Sys.time()}
  corsMatEffSize = intermediateList[[1]]
  Peffsize = intermediateList[[2]]
  corsMatPvals = intermediateList[[3]]
  Ppvals = intermediateList[[4]]
  message("Obtaining permulations p-values")
  N = nrow(realCors[[1]]) #
  #if(start = 1){ #Only do this if start = 1, because otherwise it's already made and you'll overwrite the old script's results 
  realCors[[1]]$permP = rep(NA, N) #Make a column for permP values in all of the dataframes 
  for (j in 1:length(realCors[[2]])) {
    realCors[[2]][[j]]$permP = rep(NA, N) #Make a column for permP values in all of the dataframes 
  }
  #}
  
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
      p = (sum(abs(signedMatEffSizes) > abs(realCors[[1]]$Rho[gene]), na.rm = TRUE)+1)/(sum(!is.na(signedMatEffSizes))+1)
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
        p = sum(abs(signedPeffValues) > abs(realValue), na.rm = TRUE)/ (sum(!is.na(signedPeffValues))+1)
      }
      realCors[[2]][[j]]$permP[gene] = p
    }
    if(report){geneEnd = Sys.time(); geneDuration = geneEnd - geneStart;message(paste("Completed Gene", gene, "Duration", geneDuration, attr(geneDuration, "units")))}
  }
  message("Done")
  {totalEnd = Sys.time(); totalDuration = totalEnd - totalStart;message(paste("Completed p-Values; Duration", totalDuration, attr(totalDuration, "units")))}
  return(list(res = realCors, pvals = list(corsMatPvals, Ppvals), effsize = list(corsMatEffSize, Peffsize)))}


permPValue = readRDS(permulationFileLocation) 
permPValue[16210:16310]
which(is.na(permPValue))  
plot(c(1:17000), match(names(permPValue), rownames(correlData)))

perm1 = readRDS("Output/EcholocationUpdate2/EcholocationUpdate2CombinedPrunedFast1-3400PermulationsPValue.rds")
perm2 = readRDS("Output/EcholocationUpdate2/EcholocationUpdate2CombinedPrunedFast3401-6800PermulationsPValue.rds")
perm3 = readRDS("Output/EcholocationUpdate2/EcholocationUpdate2CombinedPrunedFast6801-10200PermulationsPValue.rds")
perm4 = readRDS("Output/EcholocationUpdate2/EcholocationUpdate2CombinedPrunedFast10201-13600PermulationsPValue.rds")
perm5 = readRDS("Output/EcholocationUpdate2/EcholocationUpdate2CombinedPrunedFast13601-17000PermulationsPValue.rds")
tail(perm5)
perm5Trimmed = perm5[1:2609]
tail(perm5Trimmed)

permsAppended = append(perm1, perm2)
permsAppended = append(permsAppended, perm3)
permsAppended = append(permsAppended, perm4)
permsAppended = append(permsAppended, perm5Trimmed)

saveRDS(permsAppended, "Output/EcholocationUpdate2/EcholocationUpdate2CombinedPrunedFastAllPermulationsPValue.rds")
write.csv(correlData, "Output/EcholocationUpdate2/EcholocationUpdate2Correlations.csv")
enriment = readRDS("Output/EcholocationUpdate2/EcholocationUpdate2Enrichment-EnrichmentHsSymbolsFile2.rds")
geneset = "tissue_specific"
write.csv(enrichmentResultSets[5], paste("Output/EcholocationUpdate2/EcholocationUpdate", geneset, ".csv"))


# -- manually appending CVHRemake perm vals
perm1 = readRDS("Output/CVHRemake/CVHRemakeCombinedPrunedFast1-3242PermulationsPValue.rds")
perm2 = readRDS("Output/CVHRemake/CVHRemakeCombinedPrunedFast3243-6484PermulationsPValue.rds")
perm3 = readRDS("Output/CVHRemake/CVHRemakeCombinedPrunedFast6485-9726PermulationsPValue.rds")
perm4 = readRDS("Output/CVHRemake/CVHRemakeCombinedPrunedFast9727-12968PermulationsPValue.rds")
perm5 = readRDS("Output/CVHRemake/CVHRemakeCombinedPrunedFast12969-16210PermulationsPValue.rds")
autoAppended = readRDS("Output/CVHRemake/CVHRemakeCombinedPrunedFastAppendedPermulationsPValue.rds") 
tail(perm5)
perm5Trimmed = perm5[1:2609]
tail(perm5Trimmed)

permsAppended = append(perm1, perm2)
permsAppended = append(permsAppended, perm3)
permsAppended = append(permsAppended, perm4)
permsAppended = append(permsAppended, perm5)

all.equal(permsAppended, autoAppended)

tail(autoAppended)
tail(perm5)

saveRDS(permsAppended, "Output/EcholocationUpdate2/EcholocationUpdate2CombinedPrunedFastAllPermulationsPValue.rds")
autoAppended = readRDS("Output/CVHRemake/CVHRemakeCombinedPrunedFastAppendedPermulationsPValueOld.rds")
autoAppended = autoAppended[1:16209]



originalAppended = readRDS("Output/CVHRemake/CVHRemakeCombinedPrunedFastAllPermulationsPValue.rds") 
all.equal(names(originalAppended), names(appenedPermPValues))
all.equal(appenedPermPValues, autoAppended)


CVHCOrrelsPermAttached = readRDS("Output/CVHRemake/CVHRemakeCorrelationsFilePermulated.rds")
correlP = CVHCOrrelsPermAttached$permPValue  

which(is.na(correlP)) %in% which(is.na(autoAppended))
correlP[16209] = NA
autoAppended[16209] = NA
all.equal(correlP, autoAppended)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------



fgEdgeObjects = inputTree$edge[which(inputTree$edge.length>=1) ,]                                        #Make an object of the edges in the foreground. This is used as opposed to just referencing the tree directly to allow for "walking" in the final loop of the code
foregroundNodes = which(1:length(inputTree$tip.label) %in% as.vector(fgEdgeObjects))
foregroundSpecies = inputTree$tip.label[foregroundStartNodes]

#---- generate a phenotypeVector (named int of all tips with0/1 indicating foreground)-----
phenotypeVector = c(0,0);length(phenotypeVector) = length(inputTree$tip.label);phenotypeVector[] = 0 
names(phenotypeVector) = inputTree$tip.label
phenotypeVector[(names(phenotypeVector) %in% foregroundSpecies)] = 1
if(trimPhenotypeVector){
  phenotypeVector = phenotypeVector[names(phenotypeVector) %in% speciesFilter]
}
#Save the phenotypeVector
phenotypeVectorFilename = paste(outputFolderName, filePrefix, "phenotypeVector.rds", sep="")
saveRDS(phenotypeVector, file = phenotypeVectorFilename)


masterTree = mainTrees$masterTree

write_tree(masterTree, "Data/MasterTree.tree")


treePlotRers(mainTrees, RERObject, index = "CMBL", phenv = pathsObject)


source("Src/Reu/ZonomNameConvertMatrixCommon.R")
source("Src/Reu/ZonomNameConvertVector.R")
source("Src/Reu/ZoonomTreeNameToCommon.R")

commonRER = ZonomNameConvertMatrixCommon(RERObject)
commonPath  = ZonomNameConvertVectorCommon(pathsObject)

plotRers(commonRER, "DECR2", commonPath)
treePlotRers(mainTrees, RERObject, index = "GADD45GIP1", phenv = pathsObject)

MGICVHBin = enrichment5$MGI_Mammalian_Phenotype_Level_4
rownames(MGICVHBin)

rownames(MGICVHBin)[grep("intestinal", rownames(MGICVHBin))]

DisGenCVH = enrichment1$DisGeNET
rownames(DisGenCVH)
rownames(DisGenCVH)[grep("intestinal", rownames(DisGenCVH))]

GOProcessesCVH = enrichment3$GO_Biological_Process_2023
rownames(GOProcessesCVH)
rownames(GOProcessesCVH)[grep("digestive", rownames(GOProcessesCVH))]

rownames(GOProcessesCVH)[grep("digestive", rownames(GOProcessesCVH))]

CVHEnrichmentHs = enrichment2$EnrichmentHsSymbolsFile2
rownames(CVHEnrichmentHs)
CVHEnrichmentHs[grep("STARCH", rownames(CVHEnrichmentHs)),]


correlationData
correlationData[grep("ITBG4", rownames(correlationData))]
rownames(correlationData)
correlData[grep("ITGA6", rownames(correlData)),]


foregroundNames = readRDS(foregroundFilename)
makeMasterAndGeneTreePlots(mainTrees, "DNAH1", foregroundVector = foregroundNames)


library(readr)
scores = read_tsv("../../ZoonomiaSpecies_CarnivoryScores.tsv")



scoreValues = scores$`100`
scoreValues = append(100, scoreValues)
scoreNames = scores$Chrysochloris_asiatica
scoreNames = append("Chrysochloris_asiatica", scoreNames)
names(scoreValues) = scoreNames
saveRDS(scoreValues, "../../ZoonomiaSpecies_CanivoryScores.rds")
write.csv(scoreValues, "../../ZoonomiaSpecies_CanivoryScores.csv")

readRDS()

myVectorNmaemIchael = c(1,2,3,4,5,5)

vector = c(1,2,3,4)

vectorTheSecond = c(1,4,3,2)

vectorTheThird = c(1,2,3,4)

plot(vector, vectorTheSecond)


myVal = 5

myVal = c(1,5)


set.seed(25852)                             # Create example data
x <- rnorm(30)
y1 <- x + rnorm(30)
y2 <- x + rnorm(30, 5)

par(mar = c(5, 4, 4, 4) + 0.3)              # Additional space for second y-axis
plot(x, y1, pch = 16, col = 2)              # Create first plot
par(new = TRUE)                             # Add new plot
plot(x, y2, pch = 17, col = 3,              # Create second plot without axes
     axes = FALSE, xlab = "", ylab = "")
axis(side = 4, at = pretty(range(y2)))      # Add second axis
mtext("y2", side= 4, line = 3)             # Add second axis label


tipvals = phenotypeVector; treesObj= mainTrees; useSpecies = speciesFilter; model = modelType; anctrait = ancestralTrait; plot = T; root_prior = "auto"
source("Src/Reu/RERConvergeFunctions.R")
function (tipvals, treesObj, useSpecies = NULL, model = "ER", 
          root_prior = "auto", plot = FALSE, anctrait = NULL) 
{
  mastertree = treesObj$masterTree
  if (!is.null(useSpecies)) {
    sp.miss = setdiff(mastertree$tip.label, useSpecies)
    if (length(sp.miss) > 0) {
      message(paste0("Species from master tree not present in useSpecies: ", 
                     paste(sp.miss, collapse = ",")))
    }
    useSpecies = intersect(mastertree$tip.label, useSpecies)
    mastertree = pruneTree(mastertree, useSpecies)
    mastertree = unroot(mastertree)
  }
  else {
    mastertree = pruneTree(mastertree, intersect(mastertree$tip.label, 
                                                 names(tipvals)))
    mastertree = unroot(mastertree)
  }
  if (is.null(anctrait)) {
    tipvals <- tipvals[mastertree$tip.label]
    intlabels <- map_to_state_space(tipvals)
    print("The integer labels corresponding to each category are:")
    print(intlabels$name2index)
    ancliks = getAncLiks(mastertree, intlabels$mapped_states, 
                         rate_model = model, root_prior = root_prior)
    states = rep(0, nrow(ancliks))
    for (i in 1:length(states)) {
      states[i] = which.max(ancliks[i, ])
    }
    states = c(intlabels$mapped_states, states)
    tree = mastertree
    tree$edge.length = states[tree$edge[, 2]]
    if (length(unique(tipvals)) == 2) {
      if (sum(!unique(tipvals) %in% c(TRUE, FALSE)) > 0) {
        message("Returning categorical tree for binary phenotype because phenotype values are not TRUE/FALSE")
      }
      else {
        tree$edge.length = ifelse(tree$edge.length == 
                                    2, 1, 0)
        print("There are only 2 categories: returning a binary phenotype tree.")
        if (plot) {
          plotTree(tree)
        }
        return(tree)
      }
    }
    if (plot) {
      plotTreeCategorical(tree, category_names = intlabels$state_names, 
                          master = mastertree, node_states = states)
    }
    return(tree)
  }
  else {
    if (length(unique(tipvals)) <= 2) {
      fgspecs <- names(tipvals)[tipvals != anctrait]
      res <- foreground2Tree(fgspecs, treesObj, plotTree = plot, 
                             clade = "terminal", useSpecies = useSpecies)
      print("There are only 2 categories: returning a binary phenotype tree.")
      if (plot) {
        plotTree(res)
      }
      return(res)
    }
    else {
      tipvals <- tipvals[mastertree$tip.label]
      intlabels <- map_to_state_space(tipvals)
      j <- which(intlabels$state_names == anctrait)
      if (length(j) < 1) {
        warning("The ancestral trait provided must match one of the traits in the phenotype vector.")
      }
      res = mastertree
      res$edge.length <- rep(j, length(res$edge.length))
      traits <- intlabels$state_names
      for (trait in traits) {
        if (trait == anctrait) {
          next
        }
        i <- which(intlabels$state_names == trait)
        res$edge.length[nameEdges(res) %in% names(tipvals)[tipvals == 
                                                             trait]] = i
      }
      names(res$edge.length) = nameEdges(res)
      if (plot) {
        states = res$edge.length[order(res$edge[, 2])]
        states = c(j, states)
        plotTreeCategorical(res, category_names = traits, 
                            master = treesObj$masterTree, node_states = states)
      }
      print("Category names are mapped to integers as follows:")
      print(intlabels$name2index)
      return(res)
    }
  }
}

ERohenv = commonPhenotypeVector
erMainTrees = 
char2TreeCategorical(commonPhenotypeVector, commonMainTrees, commonSpeciesFilter, model = "ARD", anctrait = ancestralTrait, plot = T)


rm = matrix(c(1,2,3,2,4,5,6,7,8),3)
categoricalPath = char2PathsCategorical(phenotypeVector, mainTrees, speciesFilter, model = rm, anctrait = ancestralTrait, plot = T) #use the phenotype vector to make a tree

codeTable = readRDS("Data/zoonomiaToHillerCodesTable.rds")
manualAnnots = read.csv("Data/manualAnnotationsSheet.csv") 
zoonomNames = names(phenotypeVector) 
workingNames = zoonomNames


hillerNames = codeTable$TipLabel
workingNames = hillerNames
remainingIndexes = numeric()
phenotype = vector()
commonName = vector()
for(i in 1:length(hillerNames)){
  
  zonName = codeTable$zoonomiaCode[which(codeTable$TipLabel %in% hillerNames[i])]
  if(length(zonName)!= 0 & !is.na(zonName)){
    workingNames[i] = zonName
    annotRowNumber = which(manualAnnots$FaName %in% zonName)
    annotRow = manualAnnots[annotRowNumber,]
    phenotype[i] = annotRow$Meyer.Lab.Classification
    commonName[i] = annotRow$Common.Name.or.Group
  }else{
    remainingIndexes = append(remainingIndexes, i)
    phenotype[i] = NA
  }
  
}
workingNames
remainingIndexes
codeTable$`species binomial`[remainingIndexes]
length(which(is.na(phenotype)))

remainingIndexes2 = numeric()
for(i in remainingIndexes){
  hilName = hillerNames[i]
  scientificName = codeTable$ScientificName[i]
  annotRowNumber = which(manualAnnots$Species.Name %in% scientificName) 
  if(!length(annotRowNumber) == 0){
    annotRow = manualAnnots[annotRowNumber,]
    zoName = annotRow$FaName
    phenotype[i] = annotRow$Meyer.Lab.Classification
    commonName[i] = annotRow$Common.Name.or.Group
    if(!is.na(zoName) & length(zoName)!= 0 & zoName != ""){
      workingNames[i] = zoName
    }else{
      workingNames[i] = "NoZoonomiaName"
    }
  }else{
    remainingIndexes2 = append(remainingIndexes2, i)
  }
}
remainingIndexes2
length(which(is.na(phenotype)))

remainingIndexes3 = numeric()
for(i in remainingIndexes2){
  hilName = hillerNames[i]
  scientificName = codeTable$`species binomial`[i]
  annotRowNumber = which(tolower(manualAnnots$Tip_Label..Red.is.included.in.CMU.enhancer.dataset..but.missing.alignment.) %in% tolower(scientificName)) 
  if(!length(annotRowNumber) == 0){
    annotRow = manualAnnots[annotRowNumber,]
    zoName = annotRow$FaName
    phenotype[i] = annotRow$Meyer.Lab.Classification
    commonName[i] = annotRow$Common.Name.or.Group
    if(!is.na(zoName) & length(zoName)!= 0){
      workingNames[i] = zoName
    }else{
      remainingIndexes3 = append(remainingIndexes3, i)
      workingNames[i] = "NoZoonomiaName"
    }
  }else{
    remainingIndexes3 = append(remainingIndexes3, i)
  }
}
remainingIndexes2
remainingIndexes3
length(which(is.na(phenotype)))

workingNames[which(workingNames == hillerNames)] = NA

zoonomiaNames = workingNames
conversionTable = data.frame(hillerNames, zoonomiaNames, commonName, codeTable$`species binomial`, phenotype)
names(conversionTable) = c("Hiller", "Zoonomia", "common", "scientific", "phenotype")

#manually add the six missing phenotypes 
conversionTable$phenotype[which(is.na(conversionTable$phenotype))] = c("Insectivore", "Herbivore", "Piscivore", "Planktivore", NA, NA, "Herbivore", "Herbivore")
all.equal(conversionTable$phenotype[55], conversionTable$phenotype[57])
conversionTable$phenotype = trimws(conversionTable$phenotype)

saveRDS(conversionTable, "Results/HillerZoonomPhenotypeTable.rds")


conversionTable[which(is.na(conversionTable$phenotype)),]
conversionTable2 = conversionTable

dietPhenVector = conversionTable$phenotype
names(dietPhenVector) = conversionTable$Hiller

substitutions = list(c("Generalist","Anthropivore"))
if(!is.null(substitutions)){
  for( i in 1:length(substitutions)){
    substitutePhenotypes = substitutions[[i]]
    message(paste("replacing", substitutePhenotypes[1], "with", substitutePhenotypes[2]))
    dietPhenVector = gsub(substitutePhenotypes[1], substitutePhenotypes[2], dietPhenVector)
  }
}

saveRDS(dietPhenVector, "Output/OnetwentyWay6Phen/OnetwentyWay6PhenAllPhenotypesVector.rds")


if(!is.null(substitutions)){                                                    #Consider species with multiple combined categories as the merged category
  for( i in 1:length(substitutions)){                                           #Eg if [X] is replaced with [Y], [X/Y] becomes [Y]
    substitutePhenotypes = substitutions[[i]]
    message(paste("Combining", substitutePhenotypes[1], "/", substitutePhenotypes[2]))
    entriesWithPhen1 = grep(substitutePhenotypes[1], dietPhenVector)
    entriesWithPhen2 = grep(substitutePhenotypes[2], dietPhenVector)
    combineEntries = which(entriesWithPhen1 %in% entriesWithPhen2)
    combineIndexes = entriesWithPhen1[combineEntries]
    dietPhenVector[combineIndexes] = substitutePhenotypes[2]
  }
}

if(!is.null(mergeOnlys)){                                                    #Consider species with multiple combined categories as the merged category
  for( i in 1:length(mergeOnlys)){                                           #Eg if [X] is replaced with [Y], [X/Y] becomes [Y]
    substitutePhenotypes = mergeOnlys[[i]]
    message(paste("Merging Hybrids of", substitutePhenotypes[1], "/", substitutePhenotypes[2], "to", substitutePhenotypes[2]))
    entriesWithPhen1 = grep(substitutePhenotypes[1],  dietPhenVector)
    entriesWithPhen2 = grep(substitutePhenotypes[2],  dietPhenVector)
    combineEntries = which(entriesWithPhen1 %in% entriesWithPhen2)
    combineIndexes = entriesWithPhen1[combineEntries]
    dietPhenVector[combineIndexes] = substitutePhenotypes[2]
  }
}

if(!is.null(substitutions)){
  for( i in 1:length(substitutions)){
    substitutePhenotypes = substitutions[[i]]
    message(paste("replacing", substitutePhenotypes[1], "with", substitutePhenotypes[2]))
    dietPhenVector = gsub(substitutePhenotypes[1], substitutePhenotypes[2], dietPhenVector)
  }
}
saveRDS(dietPhenVector, phenotypeVectorFilename)


names(dietPhenVector)[which(is.na(dietPhenVector))]

dietPhenVector[which(is.na(dietPhenVector))] = c("Insectivore", "Herbivore", "Piscivore", "Planktivore", NA, NA, "Herbivore", "Herbivore")

conversionTable$phenotype[which(is.na(conversionTable$phenotype))] = c("Insectivore", "Herbivore", "Piscivore", "Planktivore", NA, NA, "Herbivore", "Herbivore")
all.equal(conversionTable$phenotype[55], conversionTable$phenotype[57])
conversionTable$phenotype = trimws(conversionTable$phenotype)

codeTable$`species binomial`[which(workingNames == "NoZoonomiaName")]
codeTable$`species binomial`[remainingIndexes3]




codeTable[remainingIndexes3,]
grep("ornAna", codeTable$TipLabel)



codeTable = readRDS("Data/zoonomiaToHillerCodesTable.rds")
zoonomNames = names(phenotypeVector) 
workingNames = zoonomNames  
for(i in 1:length(zoonomNames)){
  
  hilName = codeTable$TipLabel[which(codeTable$zoonomiaCode %in% zoonomNames[i])]
  if(length(hilName)!= 0){
    workingNames[i] = hilName
  }
}
remainingIndexes = which(workingNames == zoonomNames)
remainingIndexes2 = remainingIndexes
for(i in remainingIndexes){
  annotRow = manualAnnots[manualAnnots$FaName %in% zoonomNames[i],]
  scientificName = annotRow$Tip_Label..Red.is.included.in.CMU.enhancer.dataset..but.missing.alignment.
  hilName = codeTable$TipLabel[which(tolower(codeTable$`species binomial`) %in% tolower(scientificName))]
  if(length(hilName)!= 0){
    workingNames[i] = hilName
    message(i)
    remainingIndexes2 = remainingIndexes2[-which(remainingIndexes2 == i)]
  }
}

workingNames = workingNames
names(workingNames) = zoonomNames

hillerNames = workingNames
hillerNames[which(workingNames == names(workingNames))] = NA
length(hillerNames)
saveRDS(hillerNames, "Results/hillerToZoonomiaDietConversion.rds")


names(phenotypeVector) = hillerNames
phenotypeVector2 = phenotypeVector[-which(is.na(names(phenotypeVector)))]
length(phenotypeVector2)
phenotypeVector = phenotypeVector2


remainingIndexes2 

zoonomNames[remainingIndexes]



remaining120s = codeTable$TipLabel[!codeTable$TipLabel %in% workingNames]
codeTable$`species binomial`[codeTable$TipLabel %in% remaining120s]
remaining120s = remaining120s[-1]
remaining120s = remaining120s[-2]
remaining120s = remaining120s[-3]
remainingIndexes = remainingIndexes[-c(1:3)]


grep("dip", zoonomNames)
zoonomNames

remaining120s[which(!remaining120s %in% manualAnnots$Genome)]
remainingScNames = codeTable$`species binomial`[codeTable$TipLabel %in% remaining120s[(!remaining120s %in% manualAnnots$Genome)]]
possibleZonoms = manualAnnots$FaName[which(tolower(manualAnnots$Tip_Label..Red.is.included.in.CMU.enhancer.dataset..but.missing.alignment.) %in% tolower(remainingScNames))]

zoonomNames[remainingIndexes][which(zoonomNames[remainingIndexes] %in% possibleZonoms)]

which(manualAnnots$FaName %in% zoonomNames[remainingIndexes][which(zoonomNames[remainingIndexes] %in% possibleZonoms)])

which(zoonomNames == "vs_HLproCap3")





hillerNames = readRDS("Data/hillerZoonomiaCommonTable.rds")
hillerTable = data.frame(hillerNames, names(hillerNames))
names(hillerTable) = c("Hiller", "Zoonomia")
hillerTable$common = ZonomNameConvertVectorCommon(hillerTable$Zoonomia)

saveRDS(hillerTable, "Results/hillerZoonomiaCommonTable.rds")

hillerTable = readRDS("Data/HillerZoonomPhenotypeTable.rds")
commonPhentypeVector = readRDS("Results/OnetwentyWay6PhenCategoricalPhenotypeVector.rds")
commonTips = hillerTable$common[match(masterTips, hillerTable$Hiller)]
commonPhenotypeVecNames =  hillerTable$common[match(names(phenotypeVector), hillerTable$Hiller)]

phenotypeVector = commonPhenotypeVector
commonMainTrees = mainTrees
commonMainTrees$masterTree$tip.label = commonTips
commonPhenotypeVector = phenotypeVector
names(commonPhenotypeVector) = commonPhenotypeVecNames


treeImageFilename = paste(outputFolderName, filePrefix, "CategoricalTree.pdf", sep="") #make a filename based on the prefix
pdf(treeImageFilename, height = length(phenotypeVector)/18)                     #make a pdf to store the plot, sized based on tree size
char2TreeCategorical(commonPhenotypeVector, commonMainTrees, model = "ER", anctrait = ancestralTrait, plot = T)

categoricalTree = char2TreeCategorical(phenotypeVector, mainTrees, model = "ER", anctrait = ancestralTrait, plot = T) #use the phenotype vector to make a tree
dev.off()                                                                       #save the plot to the pdf

categoricalTreeFilename = paste(outputFolderName, filePrefix, "CategoricalTree.rds", sep="") #make a filename based on the prefix
saveRDS(categoricalTree, categoricalTreeFilename)  


phenotypeVector = hillerTable$phenotype
names(phenotypeVector) = hillerTable$Hiller


phenotypeVector  = readRDS("Data/CategoricalPermulationsTimingHillerPhenotypes.rds")

phenotypeVector= phenotypeVector[-which(phenotypeVector == "Insectivore")]




speciesFilter = NULL


vec = readRDS("Output/ThreePhenLikeihoodTest/ThreePhenLikeihoodTestCategoricalPhenotypeVector.rds")

catCVHDisGen = readRDS("Output/CategoricalDiet3Phen/Carnivore-Herbivore/CategoricalDiet3PhenCarnivore-HerbivoreEnrichment-DisGeNET.rds")

catCVHDisGen = catCVHDisGen$DisGeNET

catCVHDisGen = catCVHDisGen[order(catCVHDisGen$stat),]


# Make custom gene sets 

HumanLocalAdaptionDietAll = c("LCT", "FADS", "AS3MT", "DI2", "SelS", "GPX1", "GPX3", "CELF1", "SPS2", "SEPSECS", "HFE", "TRIP4", "TRVP6", "SLC30A9", "SLC39A8", "IBD5", "SLC22A4", "SLC22A5", "CREBRF")
NAFLDGWAS = c("PNPLA3", "TM6SF2", "APOE", "GCKR", "TRIB1", "GPAM", "MTARC1", "MTTP", "TOR1B", "ADH1B", "FTO", "COBLL1", "INSR", "DRG2", "GID4", "PTPRD", "PNPLA2" )
expressionDirectionalSelection = c("HLA-DQB1",
                                   "HLA-DRB1",
                                   "FADS1",
                                   "POU5F1",
                                   "HLA-DRB5",
                                   "KAT8",
                                   "HLA-DQA2",
                                   "LILRB1",
                                   "KHK",
                                   "TRIM40",
                                   "DEF8",
                                   "ZBTB12",
                                   "ZNF646",
                                   "SCAPER",
                                   "HLA-DQA1",
                                   "COMMD6",
                                   "FEN1",
                                   "SBK1",
                                   "ACO2",
                                   "ZNF668",
                                   "LY6K",
                                   "FADS3",
                                   "GSDMD",
                                   "ERBB2",
                                   "SULT1A2",
                                   "C18orf8",
                                   "LRRC61",
                                   "CLDN23",
                                   "TLR10",
                                   "LGALS2",
                                   "MAPT",
                                   "UBE2U",
                                   "NSL1",
                                   "TLR6",
                                   "SPG7",
                                   "ITGAM",
                                   "AXIN1",
                                   "PSCA",
                                   "ZNF19",
                                   "PCDHA4",
                                   "DPCR1",
                                   "STX1B",
                                   "LILRA4",
                                   "HLA-C",
                                   "HSD17B8")
genesets = list(HumanLocalAdaptionDietAll, NAFLDGWAS, expressionDirectionalSelection)
CustomGeneSet = list(genesets)
annotation

cat(HumanLocalAdaptionDietAll, sep="\t")
cat(NAFLDGWAS, sep="\t")
cat(expressionDirectionalSelection, sep="\t")




enrichmentBackup1 = enrichmentResult[[1]]
rownames(enrichmentBackup1) = paste(rownames(enrichmentBackup1), "-OvH", sep = "")

enrichmentBackup2 = enrichmentResult[[1]]
rownames(enrichmentBackup2) = paste(rownames(enrichmentBackup2), "-CvH", sep = "")

enrichmentBackup3 = enrichmentResult[[1]]
rownames(enrichmentBackup3) = paste(rownames(enrichmentBackup3), "-OvC", sep = "")

enrichmentBackup4 = enrichmentResult[[1]]
rownames(enrichmentBackup4) = paste(rownames(enrichmentBackup4), "-Overall", sep = "")

enrichmentResulta = rbind(enrichmentBackup1, enrichmentBackup2, enrichmentBackup3, enrichmentBackup4)

enrichmentResult = list(enrichmentResulta)

"Output/LiverExpression3/LiverExpression3CorrelationDataPermulatedNamesConverted.rds"


tissueEnrich = readRDS("Output/CVHRemake/CVHRemakeEnrichment-tissue_specific.rds")
write.csv(tissueEnrich, "Output/CVHRemake/CVHRemakeEnrichment-tissue_specific.csv")


oldDisG = read.gmt("Data/DisGeNET.gmt")
newDisG = read.gmt("Data/DisGeNETTest.gmt")
newerDisG = read.gmt("Data/DisGegene_associations.gmt")

all.equal(oldDisG, newDisG)


saveRDS(mainTrees, "Data/UNICORNsDemo.rds")
mainTrees$masterTree$tip.label



toothData = read.csv("Data/ToothData.csv")

toothData$FaName = toothData$Taxon

toothData$FaName = sub(" ", "_", toothData$FaName)

write.csv(toothData, "Results/ToothData.csv")



