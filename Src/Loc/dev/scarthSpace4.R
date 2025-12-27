a = b #this is to prevent accidental full runs

palette(c("yellowgreen", "darkgray", "yellow", "darkgreen", "darkblue", "lightblue", "gold", "black", "pink", "red"))
palette(c("yellowgreen", "yellow", "darkgreen", "darkblue", "lightblue", "gold", "black", "pink", "red"))


palette(c("yellow", "darkgreen", "darkblue", "lightblue", "black", "pink", "red"))
palette(c( "darkgreen", "darkblue", "lightblue", "black", "red"))
palette(c( "darkgreen", "darkblue", "lightblue", "black", "pink", "red"))
palette(c( "darkgreen", "darkblue", "black", "red"))
palette(c(  "red", "darkgreen", "black"))


palette(c( "darkgreen", "black", "darkblue", "red"))
palette(c( "darkgreen", "darkblue", "black", "red"))
palette(c("black", "darkblue"))

palette(c( "darkgreen", "darkblue", "black", "red", "gray"))
palette(c( "darkgreen", "blue", "pink", "red"))
palette(c( "darkgreen", "blue", "purple", "red"))
palette(c( "gray", "darkblue", "darkgreen",  "black", "red"))

library(RERconverge)
source("Src/Reu/ZoonomTreeNameToCommon.R")
library(data.table)
# ---------------------------------------
length(phenotypeVector)


?correlateWithCategoricalPhenotype
# -----------------------------------------
stableMaintrees = mainTrees
mainTrees = stableMaintrees
stableMaintrees = readRDS(mainTreesLocation)

mainTrees$masterTree$edge.length[1:length(mainTrees$masterTree$edge.length)] = 1

stableCommonMainTrees = stableMaintrees
stableCommonMainTrees$masterTree = ZoonomTreeNameToCommon(stableCommonMainTrees$masterTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)

?plotTreeCategorical
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



# ----- Determine gene in both H-V and H-I 

RERResults = readRDS("Output/CategoricalInsvertivoreTree/CategoricalInsVertivoreTreePairwiseCorrelationFile.rds")

HICorrelations = RERResults$`Herbivore - Insectivore`
HVCorrelations = RERResults$`Herbivore - Vertivore`
HOCorrelations = RERResults$`Herbivore - Omnivore`


HISignificantGenes = rownames(HICorrelations)[which(HICorrelations$p.adj < 0.05)]
HVSignificantGenes = rownames(HVCorrelations)[which(HVCorrelations$p.adj < 0.05)]
HOSignificantGenes = rownames(HOCorrelations)[which(HOCorrelations$p.adj < 0.05)]

sharedGenes = HISignificantGenes[which(HISignificantGenes %in% HVSignificantGenes)]

triSharedGenes = HOSignificantGenes[which(HOSignificantGenes %in% sharedGenes)]


HIDrivingData = readRDS("Output/CategoricalInsvertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreDirectionalityTable.rds")
HIDrivingData[which(rownames(HIDrivingData) %in% sharedGenes),]


HVDrivingData = read.csv("Output/CategoricalInsvertivoreTree/Herbivore-Vertivore/CategoricalInsVertivoreTreeHerbivore-VertivoreDirectionalityTable.csv")
rownames(HVDrivingData) = HVDrivingData$X
HVDrivingData[which(rownames(HVDrivingData) %in% sharedGenes),]

CombinedDrivingData = HIDrivingData
colnames(CombinedDrivingData)[10] = "HIDirectionality"
colnames(CombinedDrivingData)[11] = "HIDirectionalityNumeric"
CombinedDrivingData = cbind(CombinedDrivingData, HVDrivingData[,c(11,12)])
colnames(CombinedDrivingData)[12] = "HVDirectionality"
colnames(CombinedDrivingData)[13] = "HVDirectionalityNumeric"

all.equal(rownames(HVDrivingData), rownames(HIDrivingData))


SharedCombinedDrivingData = CombinedDrivingData[which(rownames(CombinedDrivingData) %in% sharedGenes),]

nrow(SharedCombinedDrivingData[which(SharedCombinedDrivingData$HIDirectionality == "Herbivore" & SharedCombinedDrivingData$HVDirectionality =="Herbivore"),])
nrow(SharedCombinedDrivingData[which(SharedCombinedDrivingData$HIDirectionality == "Herbivore" | SharedCombinedDrivingData$HVDirectionality =="Herbivore"),])


HISharedDriving = HIDrivingData[which(rownames(HIDrivingData) %in% sharedGenes),]
HVSharedDriving = HVDrivingData[which(rownames(HVDrivingData) %in% sharedGenes),]

rownames(HISharedDriving[which(HISharedDriving$directionality == "Herbivore"),])
rownames(HVSharedDriving[which(HVSharedDriving$directionality == "Herbivore"),])



length(which(rownames(HISharedDriving[which(HISharedDriving$directionality == "Herbivore"),]) %in% rownames(HVSharedDriving[which(HVSharedDriving$directionality == "Herbivore"),])))


# ---- Exmaine SDS and SDSL RER plots

library(RERconverge)
RERobject = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeRERFile.rds")
commonRERS = RERobject
colnames(commonRERS) = ZonomNameConvertVectorCommon(colnames(commonRERS), tipColumn = "ZoonomiaTip" )

pathsObject = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalPathsFile.rds")
mainTrees = readRDS("Data/zoonomiaAllMammalsTrees.rds")
phenotypeSet = c("Herbivore", "Insectivore", "Omnivore")
phenotypeSet = c("Herbivore", "Insectivore", "Omnivore", "Vertivore")
palette(c( "darkgreen", "darkblue", "black", "red"))

?plotRers
plotRers(RERobject, "BPIFB1", pathsObject)

source("Src/Reu/rerViolinPlot.R")
rerViolinPlot(mainTrees, RERobject, pathsObject, phenotypeSet , geneOfInterest = "SDS", colorScale = colorset)


rerViolinPlot(mainTrees, RERobject, pathsObject, phenotypeSet , geneOfInterest = "STXBP1", colorScale = colorset)
rerViolinPlot(mainTrees, RERobject, pathsObject, phenotypeSet , geneOfInterest = "GABRB3", colorScale = colorset)
rerViolinPlot(mainTrees, RERobject, pathsObject, phenotypeSet , geneOfInterest = "ALDH4A1", colorScale = colorset)


colorset = palette(c( "darkgreen", "darkblue", "white", "red"))
plotRers(commonRERS, "SDS", pathsObject)
plotRers(commonRERS, "SDSL", pathsObject)

speciesFilter = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeSpeciesFilter.rds")

?returnRersAsTree
rerTree = returnRersAsTree(mainTrees, commonRERS, "SDS", pathsObject)
rerTree = returnRersAsTree(mainTrees, RERobject, "SDS", pathsObject)

rerTree2 = drop.tip(rerTree, rerTree$tip.label[which(!rerTree$tip.label %in% speciesFilter)])

plot.phylo(rerTree2)
rerTree2$edge.length

plot(phenMasterTree)

# --- SDS and SDSL tree analysis 

source("Src/Reu/paths2Tree.R")
colorset = palette(c( "darkgreen", "darkblue", "black", "red"))
pathTree = paths2Tree(mainTrees, pathsObject, "SDS")

which(pathTree$tip.label %in% speciesFilter)

pathTree2 = drop.tip(pathTree, pathTree$tip.label[which(!pathTree$tip.label %in% speciesFilter)])

plot.phylo(pathTree2)

trgene = pathTree2
phenv = pathsObject
{
  par(mar = c(1, 1, 1, 0))
  edgcols <- rep("black", nrow(trgene$edge))
  edgwds <- rep(1, nrow(trgene$edge))
  plot.phylo(trgene, font = 2, edge.color = edgcols, edge.width = edgwds, 
             cex = tip.cex)
  rerlab <- round(rertree, 3)
  rerlab[is.na(rerlab)] <- nalab
  edgelabels(rerlab, bg = NULL, adj = c(0.5, 0.9), col = edgcols, 
             frame = "none", cex = rer.cex, font = 2)
}

plot(mainTrees$trees$SDS)
SDSPruned = drop.tip(mainTrees$trees$SDS, mainTrees$trees$SDS$tip.label[which(!mainTrees$trees$SDS$tip.label %in% speciesFilter)])
SDSPrunedCommon = ZoonomTreeNameToCommon(SDSPruned, tipCol = "ZoonomiaTip")

plot(SDSPrunedCommon)
plot(mainTrees$trees$SDSL)

# --- SDS and SDSL tree analysis 
source("Src/Reu/treeColorPlots.R")
source("Src/Reu/makePhenMasterTree.R")
palette(c( "darkgreen", "darkblue", "black", "red"))
palette(c( "darkblue", "darkgreen", "black", "red"))

testPhenMaster = makePhenMasterTree("SDS", "CategoricalInsVertivoreTree", convertToCommon = F, tipCol = "ZoonomiaTip")

pdf("Results/SDSTree.pdf", height = 20, width = 10)
png("Results/SDSTree.png", height = 1800, width = 900)
treeColorByLabel(testPhenMaster)
dev.off()

palette(c( "darkgreen", "darkblue", "black", "red"))
plotRers(commonRERS, "SDS", pathsObject)
plotRers(commonRERS, "SDSL", pathsObject)

carnivoraTips = c("vs_HLcryFer2", "cs_HLhyaHya1", "vs_HLpanOnc1", "vs_HLpanLeo1", "vs_HLaclJub2", "vs_ursMar1", "vs_lepWed1", "vs_HLphoVit1", "vs_HLcalUrs1", "vs_HLzalCal1", "vs_HLtaxTax1", "vs_HLneoVis1", "vs_HLmusPut1", "vs_HLpteBra1", "vs_HLlonCan1", "vs_HLlutLut1", "vs_HLvulVul1", "vs_HLlycPlc2")
rodentiaTips = c("vs_ochPri3", "vs_HLlepTim1", "vs_HLmusAve1", "vs_HLgilGil1", "vs_HLsciCar1", "vs_HLmarFla1", "vs_HLcynGun1", "vs_HLaplRuf1", "vs_HLpedCap1", "vs_jacJac1", "vs_HLrhlPru1", "vs_HLneoLep1", "vs_mesAur1", "vs_HLondZib1", "vs_HLmicAgr2", "vs_HLarvNll1", "vs_HLpsaObe1", "vs_HLrhoOpl1", "vs_HLcasCan3", "vs_HLperLonPac1", "vs_HLcteGun1", "vs_HLmyoCoy1", "vs_HLhydHyd1", "vs_HLcavTsc1")


library(seqinr)
SDSFasta = read.fasta("Results/ENST00000257549.SDS.filt.fa")
SDSLFasta = read.fasta("Results/ENST00000403593.SDSL.filt.fa")
sub("\\t.*", "", names(SDSFasta))
SDSCanivoraFasta = SDSFasta[which(sub("\\t.*", "", names(SDSFasta)) %in% carnivoraTips)]
SDSRondentiaFasta = SDSFasta[which(sub("\\t.*", "", names(SDSFasta)) %in% rodentiaTips)]

SDSLCanivoraFasta = SDSLFasta[which(sub("\\t.*", "", names(SDSLFasta)) %in% carnivoraTips)]
SDSLRondentiaFasta = SDSLFasta[which(sub("\\t.*", "", names(SDSLFasta)) %in% rodentiaTips)]

write.fasta(SDSCanivoraFasta, sub("\\t.*", "", names(SDSCanivoraFasta)), "Results/SDSCanivoraFasta.fasta")
write.fasta(SDSRondentiaFasta, sub("\\t.*", "", names(SDSRondentiaFasta)), "Results/SDSRondentiaFasta.fasta")
write.fasta(SDSLCanivoraFasta, sub("\\t.*", "", names(SDSLCanivoraFasta)), "Results/SDSLCanivoraFasta.fasta")
write.fasta(SDSLRondentiaFasta, sub("\\t.*", "", names(SDSLRondentiaFasta)), "Results/SDSLRondentiaFasta.fasta")

mainTrees$trees$SDS

CarnivoraSDSTree = drop.tip(mainTrees$trees$SDS, mainTrees$trees$SDS$tip.label[which(!mainTrees$trees$SDS$tip.label %in% carnivoraTips)])
CarnivoraSDSLTree = drop.tip(mainTrees$trees$SDSL, mainTrees$trees$SDSL$tip.label[which(!mainTrees$trees$SDSL$tip.label %in% carnivoraTips)])

RodentiaTree = drop.tip(mainTrees$trees$SDS, mainTrees$trees$SDS$tip.label[which(!mainTrees$trees$SDS$tip.label %in% rodentiaTips)])


CarnivoraSDSTree = drop.tip(CarnivoraSDSTree, "vs_HLtaxTax1")

write.tree(CarnivoraSDSTree, "Results/CarnivoraSDSTree.txt")
write.tree(RodentiaTree)


#------ Run directionality assessment code --- 

directionalityTable = AssessRERDirection("CategoricalInsVertivoreTree", "Herbivore-Insectivore", "CategoricalBinaryHerbivoreTree", "Herbivore", "CategoricalBinaryInsectivoreTree", "Insectivore")

write.csv(directionalityTable, "Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreDirectionalityTable.csv")
saveRDS(directionalityTable, "Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreDirectionalityTable.rds")


GODirecitonality = AssessGoCategoryDirection("CategoricalInsVertivoreTree", "Herbivore-Insectivore", "KeggReactome", 0.1, 1, F)
View(GODirecitonality)
write.csv(GODirecitonality, "Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreGoDirectionalityTable.csv")


#--- split GO directionality into positive and negative

GODirecitonality = read.csv("Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreGoDirectionalityTable.csv")


GoDirectionalityPositive = GODirecitonality[which(GODirecitonality$stat > 0),]
GoDirectionalityNegative = GODirecitonality[which(GODirecitonality$stat < 0),]

GODirecitonalityPositiveColored = GODirecitonality
GODirecitonalityPositiveColored$p.adj[GODirecitonalityPositiveColored$stat < 0] = 0.11
GODirecitonalityPositiveColored$pval[GODirecitonalityPositiveColored$stat < 0] = 0.11
write.table(GODirecitonalityPositiveColored, "Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/cytoscape/CategoricalInsVertivoreTreeHerbivore-InsectivoreGoDirectionalityPositiveColoredTable.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)




write.csv(GoDirectionalityPositive, "Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreGoDirectionalityPositiveTable.csv", row.names = F)
write.csv(GoDirectionalityNegative, "Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreGoDirectionalityNegativeTable.csv", row.names = F)
write.csv(GODirecitonalityPositiveColored, "Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreGoDirectionalityPositiveColoredTable.csv", row.names = F)



# ---- Convert a directionality table ot cytoscape format 

GODirecitonality

GOCytoscape = GODirecitonality[,c(1,7,3,4,2,6)]
colnames(GOCytoscape) = c("GO.ID", "Description", "p.adj", "DriverPackagedAsQval", "Phenotype", "Gene.vals")
GOCytoscape$Phenotype = sign(GOCytoscape$Phenotype)
GOCytoscape$Phenotype[GOCytoscape$Phenotype == 1] = "+1"
write.table(GOCytoscape, "Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/cytoscape/CytoscapeInput.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

readLines("Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/cytoscape/KeggReactome.gmt")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GSEABase")
library(GSEABase)
gmtData = getGmt("Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/cytoscape/KeggReactome.gmt")

which(names(gmtData) %in% GoDirectionalityPositive$X)
gmtList = as.list(gmtData)

gmtData[which(names(gmtData) %in% GoDirectionalityPositive$X)]


gmtNames = names(gmtData)
gmtDriver = rep(NA, length(gmtNames))
gmtUpdate = data.frame(gmtNames, gmtDriver)
gmtDirections = GODirecitonality$Directionality[match(gmtNames, GODirecitonality$X)]
write.csv(gmtDirections, "Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/cytoscape/gmtDirectionColumn.csv")


gmtNames[1]
length(GODirecitonality$X)

?match



# --- run Driver for Herbivore-Vertivore -----
source("Src/Reu/AssessRERDirection.R")
source("Src/Reu/AssessGoCategoryDirection.R")
directionalityTable = AssessRERDirection("CategoricalInsVertivoreTree", "Herbivore-Vertivore", "CategoricalBinaryHerbivoreTree", "Herbivore", "CategoricalBinaryVertivoreTree", "Vertivore")

write.csv(directionalityTable, "Output/CategoricalInsVertivoreTree/Herbivore-Vertivore/CategoricalInsVertivoreTreeHerbivore-VertivoreDriverTable.csv")
saveRDS(directionalityTable, "Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-VertivoreDriverTable.rds")
GODirecitonality = AssessGoCategoryDirection("CategoricalInsVertivoreTree", "Herbivore-Vertivore", "KeggReactome", 0.1, 1, F)
write.csv(GODirecitonality, "Output/CategoricalInsVertivoreTree/Herbivore-Vertivore/CategoricalInsVertivoreTreeHerbivore-VertivoreGoDriverTable.csv")
saveRDS(GODirecitonality, "Output/CategoricalInsVertivoreTree/Herbivore-Vertivore/CategoricalInsVertivoreTreeHerbivore-VertivoreGoDriverTable.rds")


# --- run Driver for Herbivore-Insectivore -----
directionalityTable = AssessRERDirection("CategoricalInsVertivoreTree", "Herbivore-Insectivore", "CategoricalBinaryHerbivoreTree", "Herbivore", "CategoricalBinaryInsectivoreTree", "Insectivore")

write.csv(directionalityTable, "Output/CategoricalInsVertivoreTree/Herbivore-Vertivore/CategoricalInsVertivoreTreeHerbivore-VertivoreDriverTable.csv")
saveRDS(directionalityTable, "Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-VertivoreDriverTable.rds")
GODirecitonality = AssessGoCategoryDirection("CategoricalInsVertivoreTree", "Herbivore-Vertivore", "KeggReactome", 0.1, 1, F)
write.csv(GODirecitonality, "Output/CategoricalInsVertivoreTree/Herbivore-Vertivore/CategoricalInsVertivoreTreeHerbivore-VertivoreGoDriverTable.csv")
saveRDS(GODirecitonality, "Output/CategoricalInsVertivoreTree/Herbivore-Vertivore/CategoricalInsVertivoreTreeHerbivore-VertivoreGoDriverTable.rds")





# -----------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GOSemSim")
library(GOSemSim)





# ----- Hyphy and other analysis work --- 

#------- getting the Rho values of gene across mutliple analyses --- 
HerbInsCorrelations = readRDS("Output/CategoricalInsvertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreCorrelationFile.rds")
binaryInsCorrelations = readRDS("Output/CategoricalBinaryInsectivoreTree/Background-Insectivore/CategoricalBinaryInsectivoreTreeBackground-InsectivoreCorrelationFile.rds")
binaryHerbCorrelations = readRDS("Output/CategoricalBinaryHerbivoreTree/Background-Herbivore/CategoricalBinaryHerbivoreTreeBackground-HerbivoreCorrelationFile.rds")

colnames(HerbInsCorrelations) = paste0("HI_", colnames(HerbInsCorrelations))
colnames(binaryInsCorrelations) = paste0("BI_", colnames(binaryInsCorrelations))
colnames(binaryHerbCorrelations) = paste0("BH_", colnames(binaryHerbCorrelations))

combinedData = cbind(HerbInsCorrelations, binaryInsCorrelations, binaryHerbCorrelations)

rhoVals = combinedData[,c(1,4,7)]

geneOfInterest = "SDK1"

rhoVals[which(rownames(rhoVals) == geneOfInterest),]


# --- getting number of significant genes in each comparison 

HiGeneData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreCorrelationFile.rds")
HvGeneData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Vertivore/CategoricalInsVertivoreTreeHerbivore-VertivoreCorrelationFile.rds")
HoGeneData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Omnivore/CategoricalInsVertivoreTreeHerbivore-OmnivoreCorrelationFile.rds")
IoGeneData = readRDS("Output/CategoricalInsVertivoreTree/Insectivore-Omnivore/CategoricalInsVertivoreTreeInsectivore-OmnivoreCorrelationFile.rds")
IvGeneData = readRDS("Output/CategoricalInsVertivoreTree/Insectivore-Vertivore/CategoricalInsVertivoreTreeInsectivore-VertivoreCorrelationFile.rds")
OvGeneData = readRDS("Output/CategoricalInsVertivoreTree/Omnivore-Vertivore/CategoricalInsVertivoreTreeOmnivore-VertivoreCorrelationFile.rds")
DhiGeneData = readRDS("Output/CategoricalDownsampledInsvertTree/Herbivore-Insectivore/CategoricalDownsampledInsvertTreeHerbivore-InsectivoreCorrelationFile.rds")

arrangeData = function(data, prefix){
  data$index = 1:nrow(data)
  data= data[order(data$p.adj),]
  data$rank = 1:nrow(data)
  data= data[order(data$index),]
  data$index=NULL
  colnames(data) = paste0(prefix, colnames(data))
  return(data)
}

HiGeneData = arrangeData(HiGeneData, "Hi_")
HvGeneData = arrangeData(HvGeneData, "Hv_")
HoGeneData = arrangeData(HoGeneData, "Ho_")
IoGeneData = arrangeData(IoGeneData, "Io_")
IvGeneData = arrangeData(IvGeneData, "Iv_")
OvGeneData = arrangeData(OvGeneData, "Ov_")
DhiGeneData = arrangeData(DhiGeneData, "Dhi_")


combinedData = cbind(HiGeneData, HvGeneData, HoGeneData, IoGeneData, IvGeneData, OvGeneData, DhiGeneData)
combinedData$index = 1:nrow(combinedData)


combinedDataPval = combinedData[,c(2,6,10,14,18,22,26)]
combinedDataPadjVal = combinedData[,c(3,7,11,15,19,23,27)]



length(which(combinedDataPval$Dhi_p.adj < 0.05))

apply(combinedDataPval, MARGIN = 2, mean)

apply(combinedDataPval, 2, function(column) length(which(column < 0.02)))

colnames(combinedDataPval) = c("Herbivore-Insectivore", "Herbivore-Vertivore", "Herbivore-Omnivore", "Insectivore-Omnivore", "Insectivore-Vertivore", "Omnivore-Vertivore", "Downsampled Insectivore-Herbivore")

colnames(combinedDataPval) = c("H-I", "H-V", "H-O", "I-O", "I-V", "O-V", "Downsampled H-I")
combinedDataPval$`Downsampled H-I` = NULL

colnames(combinedDataPadjVal) = c("H-I", "H-V", "H-O", "I-O", "I-V", "O-V", "Downsampled H-I")
combinedDataPadjVal$`Downsampled H-I` = NULL
sigGenesData =data.frame(category = colnames(combinedDataPadjVal), value = apply(combinedDataPadjVal, 2, function(column) length(which(column < 0.02)))) 


library(ggplot2)
library(ggpattern)

ggplot(sigGenesData, aes(x = category, y = value, fill = category, pattern)) +
  geom_bar(stat = "identity", color = "black", show.legend = FALSE) +  # Base bar color
  geom_bar_pattern(
    stat = "identity",
    pattern = "stripe",  # Options: "stripe", "crosshatch", "dots", etc.
    pattern_density = 0.25,
    pattern_fill = c("darkblue", "red", "black","black", "red", "red"),  # Pattern color
    aes(pattern = Category),  # Apply pattern per category
    show.legend = FALSE
  ) +
  theme_minimal() +
  labs(title = "Significant Genes per Pairwise Analysis",
       x = "Category", y = "Number of significant genes") +
  scale_fill_manual(values = c("darkgreen", "darkgreen", "darkgreen", "darkblue", "darkblue", "black"))+
  scale_pattern_fill_manual(values = c("darkblue", "black", "red", "black", "red", "red"))




?barplot
geom_ba
CarnivoreGeneData = readRDS("Output/CategoricalPrunedCarnivoreTree/Carnivore-Herbivore/CategoricalPrunedCarnivoreTreeCarnivore-HerbivoreCorrelationFile.rds")


InsectivoreGeneData$index = 1:nrow(InsectivoreGeneData)
VertivoreGeneData$index = 1:nrow(VertivoreGeneData)
CarnivoreGeneData$index = 1:nrow(CarnivoreGeneData)

InsectivoreGeneData = InsectivoreGeneData[order(InsectivoreGeneData$p.adj),]
VertivoreGeneData = VertivoreGeneData[order(VertivoreGeneData$p.adj),]
CarnivoreGeneData = CarnivoreGeneData[order(CarnivoreGeneData$p.adj),]

InsectivoreGeneData$rank = 1:nrow(InsectivoreGeneData)
VertivoreGeneData$rank = 1:nrow(VertivoreGeneData)
CarnivoreGeneData$rank = 1:nrow(CarnivoreGeneData)

InsectivoreGeneData = InsectivoreGeneData[order(InsectivoreGeneData$index),]
VertivoreGeneData = VertivoreGeneData[order(VertivoreGeneData$index),]
CarnivoreGeneData = CarnivoreGeneData[order(CarnivoreGeneData$index),]

colnames(InsectivoreGeneData) = paste0("I_", colnames(InsectivoreGeneData))
colnames(VertivoreGeneData) = paste0("V_", colnames(VertivoreGeneData))
colnames(CarnivoreGeneData) = paste0("C_", colnames(CarnivoreGeneData))

combinedData = cbind(InsectivoreGeneData, VertivoreGeneData, CarnivoreGeneData)


# --------- Combine the Hyphy results ----

hyphyDir = "Output/CategoricalInsVertivoreTree/Hyphy"
csvFileList = list.files(path = hyphyDir, pattern = "\\.csv$")
phenotypeList = c("H", "I", "O", "V" )
mainRowList = NULL
baseOutput = read.csv(paste0(hyphyDir, "/",csvFileList[1]), row.names = 1)
baseOutput = baseOutput[,2,drop=F]

library(stringr)
for(i in csvFileList){
  currentFileName = i 
  geneName = unlist(strsplit(currentFileName, "-"))[4]
  foreGroundNum = sub(".*_(.*?)\\..*", "\\1", currentFileName)[1]
  
  print(currentFileName)
  
  inputFile = read.csv(paste0(hyphyDir, "/",currentFileName), row.names = 1)
  inputSelected = inputFile[,2,drop=F]
  colnames(inputSelected) = paste0(geneName, "_", phenotypeList[foreGroundNum], "_correctedPValue")
  
  baseOutput = merge(baseOutput, inputSelected, by = 0, all = TRUE)
}

write.csv(baseOutput, "Output/CategoricalInsvertivoreTree/Hyphy/CombineHyphy.csv")





# --- make plots to demonstrate binary trees --- 

mainTrees = readRDS("data/zoonomiaAllMammalsTrees.rds")
mainCategoricalTree = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalTree.rds")

plotTreeCategorical(mainCategoricalTree, c("Carnivore", "Herbivore", "Omnivore"), master = stableMaintrees$masterTree)
stableCommonMainTrees = stableMaintrees
stableCommonMainTrees$masterTree = ZoonomTreeNameToCommon(stableCommonMainTrees$masterTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)

commonCategoricalTree = ZoonomTreeNameToCommon(mainCategoricalTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)
phenotypeVector = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalPhenotypeVector.rds")
pdf("Results/BinaryDemoTrees.pdf", height = length(phenotypeVector)/18, width = 10)     
palette(c( "darkgreen", "gray", "gray", "gray"))
plotTreeCategorical(commonCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = stableCommonMainTrees$masterTree)
palette(c( "gray", "darkblue", "gray", "gray"))
plotTreeCategorical(commonCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = stableCommonMainTrees$masterTree)
palette(c( "gray", "gray", "black", "gray"))
plotTreeCategorical(commonCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = stableCommonMainTrees$masterTree)
palette(c( "gray", "gray", "gray", "red"))
plotTreeCategorical(commonCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = stableCommonMainTrees$masterTree)
dev.off()


# ------ plot RERs of CategoricalINsvertivore to check rho direction meaning 
library(RERconverge)
RERobject = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeRERFile.rds")
pathsObject = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalPathsFile.rds")
mainTrees = readRDS("Data/zoonomiaAllMammalsTrees.rds")
phenotypeSet = c("Herbivore", "Insectivore", "Omnivore", "Vertivore")
colorset = c( "darkgreen", "darkblue", "black", "red")

?plotRers
plotRers(RERobject, "BPIFB1", pathsObject)

source("Src/Reu/rerViolinPlot.R")
rerViolinPlot(mainTrees, RERobject, pathsObject, phenotypeSet , geneOfInterest = "BPIFB1", colorScale = colorset)
rerViolinPlot()

# -- extract genes from top results -- 


lines <- readLines("Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/TopResults.txt")

# Filter lines that start with two tab characters
filtered_lines <- grep("^\\t\\t", lines, value = TRUE)

# Print or save the filtered lines
print(filtered_lines)

genesToRun = unique(filtered_lines)
genesToRun = gsub("\t", "", genesToRun)
HIgenesToRun = genesToRun

writeLines(genesToRun, "Results/hyphyGenesToRun.txt")

parseTopResultFileToGenes = function(file){
  lines <- readLines(file)
  
  # Filter lines that start with two tab characters
  filtered_lines <- grep("^\\t\\t", lines, value = TRUE)
  
  # Print or save the filtered lines
  #print(filtered_lines)
  
  genesToRun = unique(filtered_lines)
  genesToRun = gsub("\t", "", genesToRun)
  genesToRun
}

HItopResults = "Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/TopResults.txt"
HVtopResults = "Output/CategoricalInsVertivoreTree/Herbivore-Vertivore/TopResults.txt"
IVtopResults = "Output/CategoricalInsVertivoreTree/Insectivore-Vertivore/TopResults.txt"

testOut = parseTopResultFileToGenes(HItopResults)

HVgenesToRun = parseTopResultFileToGenes(HVtopResults)
IVgenesToRun = parseTopResultFileToGenes(IVtopResults)

VgenesToRun = append(HVgenesToRun, IVgenesToRun)
VgenesToRun[duplicated(VgenesToRun)]
VgenesToRun = unique(VgenesToRun)

VOonlyGenes = VgenesToRun[VgenesToRun %in% HIgenesToRun]
writeLines(VOonlyGenes, "Results/hyphyGenesToRunVO.txt")

HVOonlyGenes = VgenesToRun[!VgenesToRun %in% HIgenesToRun]
writeLines(HVOonlyGenes, "Results/hyphyGenesToRunHVO.txt")
#-----------------------------------


# -- exmaining the maturity results ---
??rer
rerTree = returnRersAsTree(mainTrees, RERObject, index = "PTCD1", phenv = pathsObject)
treePlotRers(mainTrees, RERObject, index = "PTCD1", phenv = pathsObject, type = "color")

source("Src/Reu/treeColorPlots.R")

treeColorByLabel(phenMasterTree)
nodelabels(frame="none")


# --- making nexus trees of genes of interest ------

EHHADHTree = mainTrees$trees$EHHADH

write.nexus(EHHADHTree, "Results/EHHADHTree.nex")
write.tree(EHHADHTree, "Results/EHHADHTree.tree")
writeNexus(EHHADHTree, "Results/EHHADHTree.nex")
EHHADHTree$node.label = NULL

masterTree = mainTrees$masterTree
writeNexus(EHHADHTree, "Results/masterTree.nex")

?writeNexus


# ----- Trim fasta file to master tree ------

fasta = read_fasta("Results/ENST00000231887.EHHADH.filt.fa")

fastaTipHeaders = fasta$headers
fastaTipHeaders = sub("\\t.*", "", fastaTipHeaders)
fastaTipHeaders[fastaTipHeaders == "REFERENCE"] = "vs_hg38"

which(!fastaTipHeaders %in% masterTree$tip.label)

fasta$headers = fasta$headers[-which(!fastaTipHeaders %in% masterTree$tip.label)]
fasta$sequences = fasta$sequences[-which(!fastaTipHeaders %in% masterTree$tip.label)]

read.dna("Results/ENST00000231887.EHHADH.filt.fa")
# -------- Redo with new package ---- 


??fasta
library(seqinr)
fastaLocation = "Results/ENST00000231887.EHHADH.filt.fa"
mainTreesLocation = 'data/zoonomiaAllMammalsTrees.rds'

fasta = read.fasta(fastaLocation)

fastaTipHeaders = names(fasta)
fastaTipHeaders = sub("\\t.*", "", fastaTipHeaders)
fastaTipHeaders[fastaTipHeaders == "REFERENCE"] = "vs_hg38"

if(!exists("mainTrees")){mainTrees = readRDS(mainTreesLocation)}
masterTree = mainTrees$masterTree
noDataTips = masterTree$tip.label[!masterTree$tip.label %in% fastaTipHeaders]
masterTree = drop.tip(masterTree, noDataTips)
masterTree$node.label = NULL


fastaToDrop = which(!fastaTipHeaders %in% masterTree$tip.label)

fasta = fasta[-fastaToDrop]
fastaTipHeaders = fastaTipHeaders[-fastaToDrop]

names(fasta) = fastaTipHeaders

write.fasta(fasta, names = names(fasta), file.out = "Results/outFasta.fa")
writeNexus(masterTree, "Results/masterTreeOut.nex")
write.tree(masterTree, "Results/masterTreeOut.tree")


masterTree$tip.label[order(masterTree$tip.label)]
fastaTipHeaders[order(fastaTipHeaders)]

#
masterTree2 = mainTrees$masterTree
masterTree2$tip.label
masterTree2$node.label = NULL

write.tree(masterTree2, "Results/fullMasterTreeOut.tree")


phenotypeTree = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalTree.rds")
masterTree2$tip.label
phenotypeTree$tip.label

phenMasterTree = masterTree2
phenMasterTree = drop.tip(phenMasterTree, phenMasterTree$tip.label[!phenMasterTree$tip.label %in% phenotypeTree$tip.label])

# ---- Remake file creation code from start cleanly -------
library(seqinr)
fastaLocation = "Results/ENST00000231887.EHHADH.filt.fa"
mainTreesLocation = 'data/zoonomiaAllMammalsTrees.rds'
foregroundCategory = "1"

useManualTree = F
filePrefix = "CategoricalInsVertivoreTree"
index = "EHHADH"
phenotypeTreeLocation = "Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalTree.rds"




if(!exists("mainTrees")){mainTrees = readRDS(mainTreesLocation)}
masterTree = mainTrees$masterTree
masterTree$node.label = NULL

if(useManualTree){
  phenotypeTree = readRDS(phenotypeTreeLocation)
}else{
  source("Src/Reu/paths2Tree.R")
  outputFolderName = paste("Output/",filePrefix,"/", sep = "")
  pathsFilename = paste(outputFolderName, filePrefix, "CategoricalPathsFile.rds", sep= "") #make a filename based on the prefix
  
  pathsObject = readRDS(pathsFilename)
  pathsTree = paths2Tree(mainTrees, pathsObject, index)
  #paths tree actually currently being unused because of how the master tree phenotype matching works. Because it's relying on the trees beingthe same shape and therefore having matching node numbers, I can't use the paths -- or, at least, it's very messy to try, so I'm not.
  
  
  phenotypeTreeCategoricalLocation = paste(outputFolderName, filePrefix, "CategoricalTree.rds", sep="") #make a filename based on the prefix
  phenotypeTreeBinaryLocation = paste(outputFolderName, filePrefix, "BinaryTree.rds", sep="") #make a filename based on the prefix
  if(file.exists(phenotypeTreeCategoricalLocation)){
    fullPhenotypeTree = readRDS(phenotypeTreeCategoricalLocation)
  }else if(file.exists(phenotypeTreeBinaryLocation)){
    fullPhenotypeTree = readRDS(phenotypeTreeBinaryLocation)
  }else{
    stop("The prefix has neither a categorical or binary phenotype tree")
  }
    
}
phenotypeTree = fullPhenotypeTree


phenMasterTree = masterTree
phenMasterTree = drop.tip(phenMasterTree, phenMasterTree$tip.label[!phenMasterTree$tip.label %in% fullPhenotypeTree$tip.label])
phenMasterTree = drop.tip(phenMasterTree, phenMasterTree$tip.label[!phenMasterTree$tip.label %in% phenotypeTree$tip.label])

#add category as label to nodes
allLabels = rep("", (length(phenMasterTree$tip.label)+phenMasterTree$Nnode))
for(i in 1:length(allLabels)){
  message(i)
  parentEdge = which(phenMasterTree$edge[,2]==i)
  if(!length(parentEdge)==0){
    allLabels[i] = parentEdge
    allLabels[i] = phenotypeTree$edge.length[parentEdge]
    }
}

allLabels[which(allLabels == foregroundCategory)] = "Foreground"

tipLabels = allLabels[c(1:length(phenotypeTree$tip.label))]
originalTipValues = phenMasterTree$tip.label
phenMasterTree$tip.label = paste0(phenMasterTree$tip.label, "{", tipLabels, "}")
internalLabels = allLabels[-c(1:length(phenotypeTree$tip.label))]
phenMasterTree$node.label = paste0("{", internalLabels, "}")



# - Read Fasta file - 
fasta = read.fasta(fastaLocation)

fastaTipHeaders = names(fasta)
fastaTipHeaders = sub("\\t.*", "", fastaTipHeaders)
fastaTipHeaders[fastaTipHeaders == "REFERENCE"] = "vs_hg38"

# trim files to match eachother 
noDataTips = phenMasterTree$tip.label[!originalTipValues %in% fastaTipHeaders]
phenMasterTree = drop.tip(phenMasterTree, noDataTips)

fastaToDrop = which(!fastaTipHeaders %in% originalTipValues)
fasta = fasta[-fastaToDrop]
fastaTipHeaders = fastaTipHeaders[-fastaToDrop]
names(fasta) = fastaTipHeaders

#Write output file
write.fasta(fasta, names = names(fasta), file.out = "Results/outFasta.fa")
write.tree(phenMasterTree, "Results/masterTreeOut.tree")

fastaOut = readLines("Results/outFasta.fa")
treeOut = readLines("Results/masterTreeOut.tree")
combinedContent = c(fastaOut, treeOut)
writeLines(combinedContent, "Results/hyphyOutput.fna")

# ------------------------------------------
treeColorPlot(phenMasterTree)
treeColorByLabel(phenMasterTree)
prunedPaths = drop.tip(pathsTree, pathsTree$tip.label[!pathsTree$tip.label %in% substr(phenMasterTree$tip.label, 1, nchar(phenMasterTree$tip.label)-3)])

treeColorByLabel(prunedPaths)
treeColorPlot(prunedPaths)
plot.phylo(prunedPaths)
prunedPaths$edge.length

fasta = fasta[1:74]

#




#
?correlateWithContinuousPhenotype

pathsVector = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalPathsFile.rds")
??paths
tree = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalTree.rds")
treesObj = mainTrees
categorical = T
tree2Paths
useSpecies = NULL
source("Src/Reu/RERConvergeFunctions.R")



?RERconverge
?plotTreeCategorical()
?treePlotNew
?plotTreeHighlightBranches()


pathsVector
length(pathsVector)
rerObject = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeRERFile.rds")
ncol(rerObject)

rerObjectTrimmed = (rerObject[1:2,])

test = rbind(pathsVector, pathsVector)

returnRersAsTree(mainTrees, test, 2)
returnRersAsTree(mainTrees, pathsVector, 1)

for(i in 1:10){
 outname = paste0("path", i)
 outPath = paths2Tree(mainTrees, pathsVector, i)
 assign(outname, outPath)
}
all.equal(path1, path2)

# -------------------------------

categoricalTree = readRDS()

??paths
??tree

rawMaturityVector = readRDS("Output/MaturityLogRaw/MaturityLogRawContinuousPhenotypeVector.rds")
percentMaturityVector = readRDS("Output/MaturityLifespanPercent/MaturityLifespanPercentContinuousPhenotypeVector.rds")
rawPaths = readRDS("Output/MaturityLogRaw/MaturityLogRawContinuousPathsFile.rds")
percentPaths = readRDS("Output/MaturityLifespanPercent/MaturityLifespanPercentContinuousPathsFile.rds")

percentRERs = readRDS("Output/MaturityLifespanPercent/MaturityLifespanPercentRERFile.rds")
rawMatRers = readRDS("Output/MaturityLogRaw/MaturityLogRawRERFile.rds")
all.equal(percentRERs, rawMatRers)



commonRERs = RERObject
colnames(commonRERs) = ZonomNameConvertVectorCommon(colnames(commonRERs), annotationLocation = spreadSheetLocation, tipCol = nameColumn)

targetGene = "PTCD1"
targetGene = "ACOT8"

targetGene = "RPS8"
targetGene = "RPS3"
targetGene = "RPS17"
targetGene = "RPL27"

targetGene = "ZIC1"
targetGene = "GTSF1"
targetGene = "OR51E1"

which(colnames(commonRERs) == "Chinese River Dolphin")


testMat = commonRERs[1:2, 1:3]

testDropRER = commonRERs[,-which(colnames(commonRERs) %in% c("Chinese River Dolphin", "Star-nosed mole", "Beaver"))]
testDropPaths = pathsObject[-which(colnames(commonRERs) %in% c("Chinese River Dolphin", "Star-nosed mole","Beaver"))]

testDropRER = commonRERs[,-which(colnames(commonRERs) %in% c("Chinese River Dolphin", "Star-nosed mole"))]
testDropPaths = pathsObject[-which(colnames(commonRERs) %in% c("Chinese River Dolphin", "Star-nosed mole"))]


commonRERs[which(rownames(commonRERs) == "GTSF1"),][order(commonRERs[which(rownames(commonRERs) == "GTSF1"),])]
{
x=pathsObject
y=commonRERs[targetGene,]
names(y)==namePathsWSpecies(mainTrees$masterTree)

plot(x,y, cex.axis=1, cex.lab=1, cex.main=1, xlab="Maturity Percentage Change",
     ylab="Evolutionary Rate", main=paste("Gene",targetGene,"Pearson Correlation"),
     pch=19, cex=1)
text(x,y, labels=names(y), pos=4)
abline(lm(y~x), col='red',lwd=3)
} # plot with percentage labels

{
  x=pathsObject
  y=commonRERs[targetGene,]
  names(y)==namePathsWSpecies(mainTrees$masterTree)
  
  plot(x,y, cex.axis=1, cex.lab=1, cex.main=1, xlab="Maturity Raw Change",
       ylab="Evolutionary Rate", main=paste("Gene",targetGene,"Pearson Correlation"),
       pch=19, cex=1, ylim=c(-1, 3))
  text(x,y, labels=names(y), pos=4)
  abline(lm(y~x), col='red',lwd=3)
} #plot with Raw labels
m <- lm(y ~ x)


{
  x=testDropPaths
  y=testDropRER[targetGene,]
  names(y)==namePathsWSpecies(mainTrees$masterTree)
  
  plot(x,y, cex.axis=1, cex.lab=1, cex.main=1, xlab="Maturity Raw Change",
       ylab="Evolutionary Rate", main=paste("Gene",targetGene,"Pearson Correlation"),
       pch=19, cex=1, ylim=c(-1, 3))
  text(x,y, labels=names(y), pos=4)
  abline(lm(y~x), col='red',lwd=3)
} #plot with Raw labels
m <- lm(y ~ x)

eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                 list(a = format(unname(coef(m)[1]), digits = 2),
                      b = format(unname(coef(m)[2]), digits = 2),
                      r2 = format(summary(m)$r.squared, digits = 3)))
as.character(as.expression(eq));

{
x=rawPaths
y=percentPaths
names(y)==namePathsWSpecies(mainTrees$masterTree)

plot(x,y, cex.axis=1, cex.lab=1, cex.main=1, xlab="Lifespan Percent Change",
     ylab="Log Raw change", main=paste("Comparission of percentage and raw change"),
     pch=19, cex=1)
text(x,y, labels=names(y), pos=4)
abline(lm(y~x), col='red',lwd=3)
}
# ---- Compare the ranking of the downsampled and non-downsampled results ---- 

mainIhCorrels = readRDS("Output/CategoricalInsvertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreCorrelationFile.rds")
mainVhCorrels = readRDS("Output/CategoricalInsvertivoreTree/Herbivore-Vertivore/CategoricalInsVertivoreTreeHerbivore-VertivoreCorrelationFile.rds")

downIhCorrels = readRDS("Output/CategoricalDownsampledInsvertTree/Herbivore-Insectivore/CategoricalDownsampledInsVertTreeHerbivore-InsectivoreCorrelationFile.rds")
downVhCorrels = readRDS("Output/CategoricalDownsampledInsvertTree/Herbivore-Vertivore/CategoricalDownsampledInsVertTreeHerbivore-VertivoreCorrelationFile.rds")

halfIhCorrels = readRDS("Output/CategoricalNoMegabranchInsvertTree/Herbivore-Insectivore/CategoricalNoMegabranchInsvertTreeHerbivore-InsectivoreCorrelationFile.rds")
halfVhCorrels = readRDS("Output/CategoricalNoMegabranchInsvertTree/Herbivore-Vertivore/CategoricalNoMegabranchInsvertTreeHerbivore-VertivoreCorrelationFile.rds")


mainIhCorrels$order = rep(1:nrow(mainIhCorrels))
mainVhCorrels$order = rep(1:nrow(mainVhCorrels))
downIhCorrels$order = rep(1:nrow(downIhCorrels))
downVhCorrels$order = rep(1:nrow(downVhCorrels))
halfIhCorrels$order = rep(1:nrow(halfIhCorrels))
halfVhCorrels$order = rep(1:nrow(halfVhCorrels))


mainIhCorrels = mainIhCorrels[order(mainIhCorrels$p.adj),]
mainVhCorrels = mainVhCorrels[order(mainVhCorrels$p.adj),]
downIhCorrels = downIhCorrels[order(downIhCorrels$p.adj),]
downVhCorrels = downVhCorrels[order(downVhCorrels$p.adj),]
halfIhCorrels = halfIhCorrels[order(halfIhCorrels$p.adj),]
halfVhCorrels = halfVhCorrels[order(halfVhCorrels$p.adj),]

mainIhCorrels$rank = rep(1:nrow(mainIhCorrels))
mainVhCorrels$rank = rep(1:nrow(mainVhCorrels))
downIhCorrels$rank = rep(1:nrow(downIhCorrels))
downVhCorrels$rank = rep(1:nrow(downVhCorrels))
halfIhCorrels$rank = rep(1:nrow(halfIhCorrels))
halfVhCorrels$rank = rep(1:nrow(halfVhCorrels))

mainIhCorrels = mainIhCorrels[order(mainIhCorrels$order),]
colnames(mainIhCorrels) = paste0("mih", colnames(mainIhCorrels))
mainVhCorrels = mainVhCorrels[order(mainVhCorrels$order),]
colnames(mainVhCorrels) = paste0("mvh", colnames(mainVhCorrels))
downIhCorrels = downIhCorrels[order(downIhCorrels$order),]
colnames(downIhCorrels) = paste0("dih", colnames(downIhCorrels))
downVhCorrels = downVhCorrels[order(downVhCorrels$order),]
colnames(downVhCorrels) = paste0("dvh", colnames(downVhCorrels))
halfIhCorrels = halfIhCorrels[order(halfIhCorrels$order),]
colnames(halfIhCorrels) = paste0("hih", colnames(halfIhCorrels))
halfVhCorrels = halfVhCorrels[order(halfVhCorrels$order),]
colnames(halfVhCorrels) = paste0("hvh", colnames(halfVhCorrels))

combinedData = cbind(mainIhCorrels, mainVhCorrels, downIhCorrels, downVhCorrels, halfIhCorrels, halfVhCorrels)


plot(combinedData$mihrank, combinedData$dihrank)
plot(combinedData$mvhrank, combinedData$dvhrank)
plot(combinedData$hvhrank, combinedData$dvhrank)
plot(combinedData$mvhrank, combinedData$hvhrank)

plot(combinedData$mihp.adj, combinedData$dihp.adj)
plot(combinedData$mvhp.adj, combinedData$dvhp.adj)
plot(combinedData$hvhp.adj, combinedData$dvhp.adj)


# --- Check downsampling RERS are equal ---- 

mainRERs = readRDS("Output/CategoricalInsvertivoreTree/CategoricalInsVertivoreTreeRERFile.rds")
halfDownRERs = readRDS("Output/CategoricalNoMegabranchInsvertTree/CategoricalNoMegabranchInsvertTreeRERFile.rds")
downRERs = readRDS("Output/CategoricalDownsampledInsvertTree/CategoricalDownsampledInsvertTreeRERFile.rds")

all.equal(mainRERs, halfDownRERs)
all.equal(downRERs, halfDownRERs)

mainFilter = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeSpeciesFilterJanThirty.rds")
halfdownFilter = readRDS("Output/CategoricalNoMegabranchInsvertTree/CategoricalNoMegabranchInsvertTreeSpeciesFilter.rds")
downFilter = readRDS("Output/CategoricalDownsampledInsvertTree/CategoricalDownsampledInsvertTreeSpeciesFilter.rds")

mainFilterOld = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeSpeciesFilter.rds")

mainFilter
halfdownFilter
downFilter
all.equal(halfdownFilter, downFilter)

all.equal(mainFilter, downFilter)
all.equal(mainFilter, mainFilterOld)

# ---- Make phylogeneticaly matching downsampled trees ----
insVertTree = readRDS("Output/CategoricalInsvertivoreTree/CategoricalInsVertivoreTreeCategoricalTree.rds")
insVertPhenv = readRDS("Output/CategoricalInsvertivoreTree/CategoricalInsVertivoreTreeCategoricalPhenotypeVector.rds")
table(insVertTree$edge.length)
table(insVertPhenv)

pdf(height = 18, width = 10)  
plotTreeCategorical(commonCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = stableCommonMainTrees$masterTree)
edgelabels(col = "darkgreen", frame = "none")
dev.off()

noMegaBranchesTree = insVertTree
noMegaBranchesTree$edge.length[c(1,2,367)] = NA

categoricalTree$edge.length[c(1,2)] = 5
commonCategoricalTree$edge.length[c(1,2)] = 5

pathsFilename = paste(outputFolderName, filePrefix, "CategoricalPathsFile.rds", sep= "") #make a filename based on the prefix
paths = tree2Paths(categoricalTree, mainTrees, useSpecies = speciesFilter, categorical = T)
#char2PathsCategorical(phenotypeVector, mainTrees, speciesFilter, model = modelType, anctrait = ancestralTrait) #make a path based on the phenotype vector
saveRDS(paths, file = pathsFilename)

?tree2Paths
categoricalTree$edge.length[c(389,388, 370, 366, 365, 359, 351, 352, 350, 347, 339, 337, 110, 336, 335, 334)] = 5
categoricalTree$edge.length[c(1,2, 364, 360, 358, 342, 343, 344, 345, 346, 348, 329, 330, 331, 332, 333, 338, 145, 3, 147)] = 5



commonCategoricalTree$edge.length[c(389,388, 370, 366, 365, 359, 351, 352, 350, 347, 339, 337, 110, 336, 335, 334)] = 5
commonCategoricalTree$edge.length[c(1,2, 364, 360, 358, 342, 343, 344, 345, 346, 348, 329, 330, 331, 332, 333, 338, 145, 3, 147)] = 5



#-------------------

maturityRER = readRDS("Output/MaturityLifespanPercent/MaturityLifespanPercentRERFile.rds")
maturityPath = readRDS("Output/MaturityLifespanPercent/MaturityLifespanPercentContinuousPathsFile.rds")
maturityMaintrees= readRDS("data/newHillerMainTrees.rds")

commonMaturityRER = maturityRER
colnames(commonMaturityRER)
colnames(commonMaturityRER) = ZonomNameConvertVectorCommon(colnames(commonMaturityRER), tipColumn = "ZoonomiaName")

returnRersAsTree(maturityMaintrees, maturityRER, "NDRG4", maturityPath)
plotRers(commonMaturityRER, "NDRG4", phenv = maturityPath)



#-----------------------------------

InsectivoreGoData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreEnrichment-GO_Biological_Process_2023.rds")
write.csv(InsectivoreGoData, "Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreEnrichment-GO_Biological_Process_2023.csv")

?correlateWithBinaryPhenotype
install.packages()

#----------------------------------------

InsectivoreGoData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreEnrichment-GO_Biological_Process_2023.rds")
VertivoreGoData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Vertivore/CategoricalInsVertivoreTreeHerbivore-VertivoreEnrichment-GO_Biological_Process_2023.rds")
CarnivoreGoData = readRDS("Output/CategoricalPrunedCarnivoryTree/Carnivore-Herbivore/CategoricalPrunedCarnivoreTreeCarnivore-HerbivoreEnrichment-GO_Biological_Process_2023.rds")


#----------------------------------------
library(ggvenn)
InsectivoreGoData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreEnrichment-GO_Biological_Process_2023.rds")[[1]]
VertivoreGoData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Vertivore/CategoricalInsVertivoreTreeHerbivore-VertivoreEnrichment-GO_Biological_Process_2023.rds")[[1]]
CarnivoreGoData = readRDS("Output/CategoricalPrunedCarnivoreTree/Carnivore-Herbivore/CategoricalPrunedCarnivoreTreeCarnivore-HerbivoreEnrichment-GO_Biological_Process_2023.rds")[[1]]

InsectivoreGoData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreEnrichment-KeggReactome.rds")[[1]]
VertivoreGoData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Vertivore/CategoricalInsVertivoreTreeHerbivore-VertivoreEnrichment-KeggReactome.rds")[[1]]
CarnivoreGoData = readRDS("Output/CategoricalPrunedCarnivoreTree/Carnivore-Herbivore/CategoricalPrunedCarnivoreTreeCarnivore-HerbivoreEnrichment-KeggReactome.rds")[[1]]



InsectivoreGoData = InsectivoreGoData[order(InsectivoreGoData$p.adj),]
VertivoreGoData = VertivoreGoData[order(VertivoreGoData$p.adj),]
CarnivoreGoData = CarnivoreGoData[order(CarnivoreGoData$p.adj),]


signficiantInsectivore = InsectivoreGoData[which(InsectivoreGoData$p.adj <0.05),]
signficiantVertivore = VertivoreGoData[which(VertivoreGoData$p.adj <0.05),]
signficiantCarnivore = CarnivoreGoData[which(CarnivoreGoData$p.adj <0.05),]


valuedInsectivore = InsectivoreGoData[which(InsectivoreGoData$pval <1),]
valuedVertivore = VertivoreGoData[which(VertivoreGoData$pval <1),]
valuedCarnivore = CarnivoreGoData[which(CarnivoreGoData$pval <1),]

vennData = list(
  Insectivore = rownames(valuedInsectivore),
  Vertivore = rownames(valuedVertivore),
  Carnivore = rownames(valuedCarnivore)
)

signficianterInsectivore = InsectivoreGoData[which(InsectivoreGoData$p.adj <0.05),]
signficianterVertivore = VertivoreGoData[which(VertivoreGoData$p.adj <0.05),]
signficianterCarnivore = CarnivoreGoData[which(CarnivoreGoData$p.adj <0.05),]

topInsectivore = InsectivoreGoData[1:100,]
topVertivore = VertivoreGoData[1:100,]
topCarnivore = CarnivoreGoData[1:100,]

vennData = list(
  Insectivore = rownames(signficiantInsectivore),
  Vertivore = rownames(signficiantVertivore),
  Carnivore = rownames(signficiantCarnivore)
)
ggvenn(vennData, fill_color = c("blue", "red", "orange"))


# --------

InsectivoreGeneData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreCorrelationFile.rds")
VertivoreGeneData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Vertivore/CategoricalInsVertivoreTreeHerbivore-VertivoreCorrelationFile.rds")
CarnivoreGeneData = readRDS("Output/CategoricalPrunedCarnivoreTree/Carnivore-Herbivore/CategoricalPrunedCarnivoreTreeCarnivore-HerbivoreCorrelationFile.rds")

InsectivoreGeneData$index = 1:nrow(InsectivoreGeneData)
VertivoreGeneData$index = 1:nrow(VertivoreGeneData)
CarnivoreGeneData$index = 1:nrow(CarnivoreGeneData)

InsectivoreGeneData = InsectivoreGeneData[order(InsectivoreGeneData$p.adj),]
VertivoreGeneData = VertivoreGeneData[order(VertivoreGeneData$p.adj),]
CarnivoreGeneData = CarnivoreGeneData[order(CarnivoreGeneData$p.adj),]

InsectivoreGeneData$rank = 1:nrow(InsectivoreGeneData)
VertivoreGeneData$rank = 1:nrow(VertivoreGeneData)
CarnivoreGeneData$rank = 1:nrow(CarnivoreGeneData)

InsectivoreGeneData = InsectivoreGeneData[order(InsectivoreGeneData$index),]
VertivoreGeneData = VertivoreGeneData[order(VertivoreGeneData$index),]
CarnivoreGeneData = CarnivoreGeneData[order(CarnivoreGeneData$index),]

colnames(InsectivoreGeneData) = paste0("I_", colnames(InsectivoreGeneData))
colnames(VertivoreGeneData) = paste0("V_", colnames(VertivoreGeneData))
colnames(CarnivoreGeneData) = paste0("C_", colnames(CarnivoreGeneData))

combinedData = cbind(InsectivoreGeneData, VertivoreGeneData, CarnivoreGeneData)
combinedData$V_index = NULL
combinedData$C_index = NULL

combinedData = combinedData[order(combinedData$I_p.adj),]


equationLinePlot = function(data, xIn, yIn){
 linearModel = lm(yIn ~ xIn, data = data) 
  
 equation = paste0("y=", round(coef(linearModel)[2], 2), "*x", round(coef(linearModel)[1], 2))
 rSquared = paste("R² =", round(summary(linearModel)$r.squared, 2))
 
  ggplot(data, aes(x = xIn, y = yIn)) + 
    geom_point()+
    geom_smooth(method = "lm")  
}

equationLinePlot(combinedData, "I_Rho", "C_Rho")

linearModel = lm(I_Rho ~ C_Rho, data = combinedData) 

library(gridExtra)



linearModel = lm(-C_Rho ~ I_Rho, data = combinedData) 
equation = paste0("y=", round(coef(linearModel)[2], 2), "*x", round(coef(linearModel)[1], 2))
rSquared = paste("R² =", round(summary(linearModel)$r.squared, 2))
plot1 = ggplot(combinedData, aes(x = I_Rho, y = -C_Rho)) + 
  geom_point()+
  geom_smooth(method = "lm")+
  annotate("text", x = 3, y = 9, label = paste(equation, rSquared, sep = "\n"), color = "blue", size = 10)+
  theme_minimal()

linearModel = lm(-C_Rho ~ V_Rho, data = combinedData) 
equation = paste0("y=", round(coef(linearModel)[2], 2), "*x", round(coef(linearModel)[1], 2))
rSquared = paste("R² =", round(summary(linearModel)$r.squared, 2))
plot2 = ggplot(combinedData, aes(x = V_Rho, y = -C_Rho)) + 
  geom_point()+
  geom_smooth(method = "lm")+
  annotate("text", x = 3, y = 9, label = paste(equation, rSquared, sep = "\n"), color = "blue", size = 10)+
  theme_minimal()


grid.arrange(plot1, plot2, ncol =2)



?lm

library(ggplot2)
ggplot(combinedData, aes(x = I_Rho, y = -C_Rho)) + 
  geom_point()+
  geom_smooth(method = "lm")


plot(combinedData$I_Rho, combinedData$V_Rho)
plot(-combinedData$C_Rho, combinedData$I_Rho)
plot(-combinedData$C_Rho, combinedData$V_Rho)




plot(combinedData$I_rank, combinedData$V_rank, xlim = c(0,3500), ylim = c(0,3500))
plot(combinedData$`I_rank`, combinedData$`C_rank`, xlim = c(0,4000), ylim = c(0,4000))
plot(combinedData$`V_rank`, combinedData$`C_rank`, xlim = c(0,4000), ylim = c(0,4000))
plot(combinedData$V_index, combinedData$V_p.adj)
hist(combinedData$V_p.adj)
length(combinedData$V_p.adj[which(combinedData$V_p.adj < 1)])
length(combinedData$I_p.adj[which(combinedData$I_p.adj < 1)])

?cbind
sigInsectGenes = InsectivoreGeneData[which(InsectivoreGeneData$p.adj <0.05),]
sigVertGenes = VertivoreGeneData[which(VertivoreGeneData$p.adj <0.05),]
sigCarnGenes = CarnivoreGeneData[which(CarnivoreGeneData$p.adj <0.05),]

geneVennData = list(
  Insectivore = rownames(sigInsectGenes),
  Vertivore = rownames(sigVertGenes),
  Carnivore = rownames(sigCarnGenes)
)
ggvenn(geneVennData, fill_color = c("blue", "red", "orange"))
require(venneuler)
v <- venneuler(c(Insectivore=194, Vertivore=145, Carnviore=598, "Insectivore&Vertivore"=0, "Insectivore&Carnviore"=294, "Vertivore&Carnviore"=125, "Carnviore&Vertivore&Insectivore"=75))
plot(v)

valuedInsectGenes = combinedData[which(combinedData$I_p.adj <1),]
valuedVertGenes = combinedData[which(combinedData$V_p.adj <1),]
valuedCarnGenes = combinedData[which(combinedData$C_p.adj <1),]

geneVennData = list(
  Insectivore = rownames(valuedInsectGenes),
  Vertivore = rownames(valuedVertGenes),
  Carnivore = rownames(valuedCarnGenes)
)

ggvenn(geneVennData, fill_color = c("blue", "red", "pink"))
# ---------------------------
library(RERconverge)
RERobject = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeRERFile.rds")
pathsObject = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalPathsFile.rds")
mainTrees = readRDS("Data/zoonomiaAllMammalsTrees.rds")
phenotypeSet = c("Herbivore", "Insectivore", "Omnivore", "Vertivore")
phenotypeSet = c(" Herbivore", " Insectivore", "Omnivore", " Vertivore")
colorset = c( "darkgreen", "darkblue", "red", "black")

plotRers(RERobject, "BPIFB1", pathsObject)

source("Src/Reu/rerViolinPlot.R")
rerViolinPlot(mainTrees, RERobject, pathsObject, phenotypeSet , geneOfInterest = "SLC14A2", colorScale = colorset)







#-------------------------------------------
manualAnnotsTrimmed = manualAnnots
which(manualAnnots$ZoonomiaTip %in% names(phenotypeVector))
manualAnnotsTrimmed = manualAnnotsTrimmed[which(manualAnnots$ZoonomiaTip %in% names(phenotypeVector)), ]

table(manualAnnotsTrimmed$MSWC_Family)

# ----------------------------
lowCategoryGeneDropper(mainTrees, phenotypeVector)

# ----------------------------
categoricalCorrelation = correlateWithCategoricalPhenotype(RERObject, pathsObject, min.sp = 400, min.pos = 2) #Calculate with categorical, min 2 species per category 
overalCategorical = categoricalCorrelation[[1]]                               #select the results relating to overall difference between all categories
correlation = overalCategorical                                               # and classify it as the main correlation file

#process the pairwise outputs
pairwiseCategorical = categoricalCorrelation[[2]]                             #select the group of pairwise comparisons

phenotypeVectorFilename = paste(outputFolderName, filePrefix, "CategoricalPhenotypeVector.rds",sep="") #select the phenotype vector based on prefix
phenotypeVector = readRDS(phenotypeVectorFilename)                            #load in the phenotype vector 
categories = map_to_state_space(phenotypeVector)                              #and use it to connect branch lengths to phenotype name
categoryNames = categories$name2index                                         #store the length-phenotype connection

pairwiseTableNames = names(pairwiseCategorical)                               #Prepare to repalce the number-number titles with phenotype-phenotype titles
for(i in 1:length(categoryNames)){                                            #for each phenotype
  pairwiseTableNames= gsub(i, names(categoryNames)[i], pairwiseTableNames)                        #replace the number with the phenotype name  
}
names(pairwiseCategorical) = pairwiseTableNames                               #update the dataframe titles

pairwiseCorrelationFileName = paste(outputFolderName, filePrefix, "PairwiseCorrelationFile", sep= "") #make a name for the pairwise comparisons based on prefix
write.csv(pairwiseCategorical, file= paste(pairwiseCorrelationFileName, ".csv", sep=""), row.names = T, quote = F) #save the correlations as a csv
saveRDS(pairwiseCategorical, paste(pairwiseCorrelationFileName, ".rds", sep="")) #and as an rds 

combinedCategoricalCorrelationFilename = pairwiseCorrelationFileName = paste(outputFolderName, filePrefix, "CombinedCategoricalCorrelationFile", sep= "") # make this file for later functions that want it in combo
saveRDS(categoricalCorrelation, paste(combinedCategoricalCorrelationFilename, ".rds", sep="")) #and as an rds 

#save the outputs to subdirectories 
outputSubdirectoryNoslash = paste(outputFolderName, "Overall", sep = "")
if(!dir.exists(outputSubdirectoryNoslash)){                       #create that directory if it does not exist
  dir.create(outputSubdirectoryNoslash)
}
outputSubdirectory = paste(outputSubdirectoryNoslash, "/", sep="")

correlationsOverallFilename = paste(outputSubdirectory, filePrefix, "OverallCorrelationFile.rds", sep= "")
saveRDS(categoricalCorrelation[[1]], correlationsOverallFilename)

for(i in 1:length(pairwiseTableNames)){
  pairwiseTableNames= gsub(" ", "", pairwiseTableNames)
  
  outputSubdirectoryNoslash = paste(outputFolderName, pairwiseTableNames[i], sep = "")
  if(!dir.exists(outputSubdirectoryNoslash)){                       #create that directory if it does not exist
    dir.create(outputSubdirectoryNoslash)
  }
  outputSubdirectory = paste(outputSubdirectoryNoslash, "/", sep="")
  
  correlationsPairFilename = paste(outputSubdirectory, filePrefix, pairwiseTableNames[i], "CorrelationFile",".rds", sep= "")
  saveRDS(categoricalCorrelation[[2]][[i]], correlationsPairFilename)
}





# -----------------------

report = mainTrees$report
view(report)

test = hist(rowSums(report))

length(which(rowSums(report)<400))

colnames(report) %in% names(phenotypeVector)
reportPruned = report[,colnames(report) %in% names(phenotypeVector)]

ncol(reportPruned)
colnames(reportPruned) %in% names(phenotypeVector)

test = hist(rowSums(reportPruned))
length(which(rowSums(reportPruned)<170))

?hist()



# -------------------------

mainTrees$masterTree$edge.length[1:length(mainTrees$masterTree$edge.length)] = 1

char2TreeCategoricalStates = function (tipvals, treesObj, useSpecies = NULL, model = "ER", 
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
    return(states)
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

commonStates = char2TreeCategoricalStates(commonPhenotypeVector, commonMainTrees, commonSpeciesFilter, model = modelType, anctrait = ancestralTrait, plot = T)
states = char2TreeCategoricalStates(phenotypeVector, mainTrees, speciesFilter, model = modelType, anctrait = ancestralTrait, plot = T) #use the phenotype vector to make a tree

categoricalTree
commonCategoricalTree

length(speciesFilter)
length(commonSpeciesFilter)

all.equal(categoricalTree$edge.length, commonCategoricalTree$edge.length)
all.equal(states, commonStates)

commonCategoricalTree$tip.label[which(duplicated(commonCategoricalTree$tip.label))]

manualAnnots$CommonName[which(duplicated(manualAnnots$CommonName))]

commonCategoricalTree = ZoonomTreeNameToCommon(categoricalTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)
stableMaintrees = mainTrees
mainTrees$masterTree$edge.length[1:length(mainTrees$masterTree$edge.length)] = 1
stableMaintrees = readRDS(mainTreesLocation)
stableCommonMainTrees = stableMaintrees
stableCommonMainTrees$masterTree = ZoonomTreeNameToCommon(stableCommonMainTrees$masterTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)

?plotTreeCategorical
plotTreeCategorical(commonCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = stableCommonMainTrees$masterTree)

plotTreeCategorical(categoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = stableMaintrees$masterTree)



plotTreeCategorical(categoricalTree, c("Carnivore", "Herbivore", "Omnivore"), master = stableMaintrees$masterTree)

plotTreeCategorical(commonCategoricalTree, c("Carnivore", "Herbivore", "Omnivore"), master = stableCommonMainTrees$masterTree)



# --------------------------

phenotypeVectorFilename = paste(outputFolderName, filePrefix, "CategoricalPhenotypeVector.rds",sep="") #make a filename based on the prefix
phenotypeVector = readRDS(phenotypeVectorFilename)

lowCategoryGeneDropper = function(mainTrees, phenotypeVector){
  genesToDrop = vector()
  for(i in 1:length(mainTrees$trees)){
    currentTree = mainTrees$trees[[i]]
    currentTreeName = names(mainTrees$trees[i])
    currentTips = currentTree$tip.label
    
    phenotypedTips = currentTips[which(currentTips %in% names(phenotypeVector))]
    
    phenotypeValues = phenotypeVector[match(phenotypedTips, names(phenotypeVector))]
    phenotypeNumbers = table(phenotypeValues)
    if(any(phenotypeNumbers <3)){
      message(currentTreeName)
      print(phenotypeNumbers)
      genesToDrop = append(genesToDrop, currentTreeName)
    }
  }
  return(genesToDrop)
}  

lowCategoryGeneDropper(mainTrees, phenotypeVector)

# -------------------
ZoonomTreeNameToCommon(commonMainTrees$masterTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)
commonRERs = RERObject

colnames(commonRERs) = ZonomNameConvertVectorCommon(colnames(commonRERs), annotationLocation = spreadSheetLocation, tipCol = nameColumn)

plotRers(commonRERs, "ECI2", pathsObject)

?plotRers



# -----------------
testTree = readRDS("Output/CategoricalMobivoreTree/CategoricalMobivoreTreeCategoricalTree.rds")
table(testTree$edge.length)


REROne = readRDS("Output/CategoricalMobivoreTree/CategoricalMobivoreTreeRERFile-Other.rds")
RERTwo = readRDS("Output/CategoricalMobivoreTree/CategoricalMobivoreTreeRERFileOLd.rds")

all.equal(REROne, RERTwo)

demoTree = readRDS
mainTrees$masterTree$tip.label[mainTrees$masterTree$tip.label %in% "vs_OrnAna3"]

mainTrees$masterTree$edge.length[1:length(mainTrees$masterTree$edge.length)] = 1


testMain = readRDS("Data/zoonomiaAllMammalsTrees.rds")
testMain$master

grep("Ana", mainTrees$masterTree$tip.label)

phenotypeVector = readRDS(phenotypeVectorFilename)

categoricalTree$tip.label[!categoricalTree$tip.label %in% testTree$tip.label]
testTree$tip.label[!testTree$tip.label %in% categoricalTree$tip.label]
length(commonSpeciesFilter)

commonPhenotypeVector[names(commonPhenotypeVector) %in% "Platypus"]

prunedTree$edge.length

mainTrees2 = readRDS("Data/zoonomiaAllMammalsTrees.rds")
length(mainTrees2$masterTree$tip.label)

all.equal(mainTrees$masterTree, mainTrees2$masterTree)

length(mainTrees$masterTree$tip.label)
grep("ornAna", mainTrees$masterTree$tip.label)
mainTrees$masterTree$edge.length


testTrees = readRDS("data/RemadeTreesAllZoonomiaSpecies.rds")
testTrees$masterTree$edge.length

which(names(phenotypeVector) == "vs_HLlniGeo1")
phenotypeVector[94]

plotTreeCategorical2 = function (tree, category_names = NULL, master = NULL, node_states = NULL) 
{
  n = length(unique(tree$edge.length))
  if (n > length(palette())) {
    colors = colorRampPalette(palette())(n)
  }
  else {
    colors = c("yellowgreen", "lightgray", "yellow", "darkgreen", "darkblue", "lightblue", "gold", "black", "pink", "red")
  }
  edge_colors = tree$edge.length
  edge_colors = sapply(edge_colors, function(x) {
    colors[x]
  })
  par(mar = c(5, 4, 4, 10), xpd = TRUE)
  if (!is.null(master)) {
    cm = intersect(master$tip.label, tree$tip.label)
    master = pruneTree(master, cm)
    if (!is.null(node_states)) {
      node_colors = node_states
      node_colors = sapply(node_colors, function(x) {
        colors[x]
      })
      plot(master, cex = 0.25, edge.color = edge_colors, 
           node.color = node_colors)
    }
    else {
      plot(master, cex = 0.25, edge.color = edge_colors)
    }
  }
  else {
    if (!is.null(node_states)) {
      node_colors = node_states
      node_colors = sapply(node_colors, function(x) {
        colors[x]
      })
      plot(tree, cex = 0.25, edge.color = edge_colors, 
           use.edge.length = FALSE, node.depth = 2, node.color = node_colors)
    }
    else {
      plot(tree, cex = 0.25, edge.color = edge_colors, 
           use.edge.length = FALSE, node.depth = 2)
    }
  }
  if (!is.null(category_names)) {
    legend(x = "bottomright", inset = c(-0.25, 0), cex = 0.5, 
           legend = category_names, col = colors, lwd = 2)
  }
}

char2TreeCategorical2 = function (tipvals, treesObj, useSpecies = NULL, model = "ER", 
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
      plotTreeCategorical2(tree, category_names = intlabels$state_names, 
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



treeImageFilename = paste(outputFolderName, filePrefix, "CategoricalTree.pdf", sep="") #make a filename based on the prefix
pdf(treeImageFilename, height = length(phenotypeVector)/18)                     #make a pdf to store the plot, sized based on tree size
char2TreeCategorical2(commonPhenotypeVector, commonMainTrees, commonSpeciesFilter, model = modelType, anctrait = ancestralTrait, plot = T)

categoricalTree = char2TreeCategorical2(phenotypeVector, mainTrees, speciesFilter, model = modelType, anctrait = ancestralTrait, plot = T) #use the phenotype vector to make a tree
dev.off()   



mainTrees3 = read.tree("Results/NewZoonomiaMasterTreePrunedToAlignmentSpecies.nwk")







# ------------------------------

fullTree = readRDS("Data/zoonomiaAllMammalsTrees.rds")
plotTree(commonMainTrees$masterTree)

fullTreeTrees = fullTree[[1]]

fullTree$masterTree

tipNumberList = numeric()
for(i in 1:length(fullTreeTrees)){
  currentTipNumber = length(fullTreeTrees[[i]]$tip.label)
  message(currentTipNumber)
  tipNumberList = append(tipNumberList, currentTipNumber)
}
length(tipNumberList)

max(tipNumberList)

tipList = list()
for(i in 1:length(fullTreeTrees)){
  currentTips = fullTreeTrees[[i]]$tip.label
  tipList = append(tipList, list(currentTips))
}

tipNumberList[order(tipNumberList, decreasing = T)]

plotTree(mainTrees$masterTree)
mainTrees$masterTree$tip.label
grep("ana", mainTrees$masterTree$tip.label)


biggestTree = fullTreeTrees[[8866]]$tip.label

singleMissing = fullTree$masterTree$tip.label[which(!fullTree$masterTree$tip.label %in% biggestTree)]

fullTreeTips = fullTree$masterTree$tip.label


tipTreeNumber = numeric()
for(i in 1:length(fullTreeTips)){
  currentTip = fullTreeTips[i]
  currentTipTreeNumber = length(which(sapply(tipList, function(x) currentTip %in% x)))
  names(currentTipTreeNumber) = currentTip
  message(currentTipTreeNumber)
  tipTreeNumber = append(tipTreeNumber, currentTipTreeNumber)
}

lowTipTrees = tipTreeNumber[order(tipTreeNumber)]



which(sapply(tipList, length) > 467)
highTipGenes = tipList[which(sapply(tipList, length) > 467)]

n=1
TestTip1 = names(lowTipTrees[1])
TestTip2 = names(lowTipTrees[2])


which(sapply(highTipGenes, function(x) !TestTip1 %in% x))
which(sapply(highTipGenes, function(x) !TestTip2 %in% x))

tipListTestDrop = tipList

lowTipTreesDropping = lowTipTrees[lowTipTrees < 10000]
tipsToDrop = names(lowTipTreesDropping) 
length(tipsToDrop)

tipListTestDrop = lapply(tipListTestDrop, function(x) Filter(function(y) !(y %in% tipsToDrop), x))


dropedLengths = sapply(tipListTestDrop, length)
dropedLengths[order(dropedLengths, decreasing = T)]

tressWith458 = which(dropedLengths == 458)

all10kspeciesTrees = fullTreeTrees[tressWith458]


all10kspeciesTreesSameTest = all10kspeciesTrees

all10kspeciesTreesSameTest$KAT7$edge.length = rep(1, length(all10kspeciesTreesSameTest$KAT7$edge.length))
all10kspeciesTreesSameTest$METTL1$edge.length = rep(1, length(all10kspeciesTreesSameTest$KAT7$edge.length))
all10kspeciesTreesSameTest$WNT2B$edge.length = rep(1, length(all10kspeciesTreesSameTest$KAT7$edge.length))
all10kspeciesTreesSameTest$TGFBI$edge.length = rep(1, length(all10kspeciesTreesSameTest$KAT7$edge.length))
all10kspeciesTreesSameTest$CFAP97D1$edge.length = rep(1, length(all10kspeciesTreesSameTest$KAT7$edge.length))


all.equal(all10kspeciesTreesSameTest$KAT7, all10kspeciesTreesSameTest$METTL1)
all.equal(all10kspeciesTreesSameTest$KAT7, all10kspeciesTreesSameTest$WNT2B)
all.equal(all10kspeciesTreesSameTest$KAT7, all10kspeciesTreesSameTest$TGFBI)
all.equal(all10kspeciesTreesSameTest$KAT7, all10kspeciesTreesSameTest$CFAP97D1)


report = fullTree$report
#write.csv(report, file= "Results/geneTreesReport.csv")


sum(report$vs_HLthyCyn1)

# ---

report = read.csv("Results/geneTreesReport.csv")
rownames(report) = report$X
report = report[,-1]

numSpecies = rowSums(report)[order(rowSums(report), decreasing = T)]
test = table(numSpecies)

topGeneNames = names(numSpecies[1:length(which(numSpecies >453))])
topGeneNames = names(numSpecies)

topGenes = report[which(row.names(report)%in% topGeneNames),]

speciesInTopTrees = colSums(topGenes)[order(colSums(topGenes))]
speciesMissingFromTopTrees = speciesInTopTrees[speciesInTopTrees < length(topGeneNames)]
length(speciesMissingFromTopTrees)


speciesToKeep = c("vs_HLornAna3", "vs_HLtacAcu1", "vs_HLgymLea1", "vs_HLpseCup1", "vs_ptePar1", "vs_HLpseCor1")


testDrop = topGenes
i=1
tipsToDrop = character()
while(T){
  currentLowestSpecies = names(speciesMissingFromTopTrees[i])
  
  message(" -------------- ")
  message(" i = ", i )
  message(currentLowestSpecies)
  if(!currentLowestSpecies %in% speciesToKeep){
  #if(T){
    tipsToDrop = append(tipsToDrop, currentLowestSpecies)
    dropCol = which(colnames(testDrop) == currentLowestSpecies)
    
    testDrop = testDrop[,-dropCol]
  }

  
  #check
  ncol(testDrop)
  print(rowSums(testDrop)[order(rowSums(testDrop))])
  
  numberOFFullTrees = length(which(rowSums(testDrop) == ncol(testDrop)))
  message(paste("number of matching trees =",numberOFFullTrees))
  
  if(numberOFFullTrees >9){
    message("Tips to drop:")
    print(tipsToDrop)
    break()
  }else(
    i = i+1
  )
}

tipsToDrop1 = tipsToDrop
tipsToDrop2 = tipsToDrop
tipsToDrop3 = tipsToDrop
tipsToDrop4 = tipsToDrop 

tipsToDropFinal = tipsToDrop


ZonomNameConvertVectorCommon(tipsToDrop, tipColumn = "Zoonomia")


tipsToDrop1 %in% tipsToDrop2
tipsToDrop4 %in% tipsToDrop1



tipsInNewMaster3 = colnames(testDrop)


tipsInNewMaster2 = colnames(testDrop)

tipsInNewMaster = colnames(testDrop)

all.equal(tipsInNewMaster3, tipsInNewMaster2)

saveRDS(tipsInNewMaster2, "Results/newZoMasterTips.rds")

tipsInUpdatedMaster = colnames(testDrop)


length(tipsInUpdatedMaster)
length(tipsInNewMaster)


togaTree = read.tree("Data/togaTree.nwk")
plot.phylo(togaTree)

tipsToDrop = togaTree$tip.label[!togaTree$tip.label %in% tipsInNewMaster2]

togaPruned = drop.tip(togaTree, tipsToDrop)
plot.phylo(togaPruned)

write.tree(togaPruned, "Results/NewZoonomiaMasterTreePrunedToAlignmentSpecies.nwk")

tipsToDrop = togaTree$tip.label[!togaTree$tip.label %in% tipsInUpdatedMaster]
togaPruned = drop.tip(togaTree, tipsToDrop)
plot.phylo(togaPruned)
write.tree(togaPruned, "Results/NewZoonomiaMasterTreePrunedToAlignmentSpecies.nwk")


#

test
?hist

dev.off()

which(is.na(commonMainTrees$masterTree$tip.label))
testPlot = plotTree(mainTrees$masterTree)

pdf(file = "test.pdf", width = 1000, height = 1000)
dev.off()
commonMainTrees$masterTree = commonMainTrees$masterTree

length(mainTrees$masterTree$edge.length)
which(is.na(mainTrees$masterTree$edge.length))

plot.phylo(mainTrees$masterTree)

mainTrees$masterTree$edge[which(is.na(mainTrees$masterTree$edge.length)),]

length(mainTrees$masterTree$tip.label)

NAFindTree = mainTrees$masterTree

NAFindTree$edge.length[which(is.na(NAFindTree$edge.length))] = 0.123456

plotTree(NAFindTree)

pdf("output/TestTree.pdf", height = length(NAFindTree$tip.label)/18)  

plotTreeHighlightBranches(NAFindTree, hlspecies = which(NAFindTree$edge.length == 0.123456), hlcols = "blue")

dev.off()

?plotTree


hillerMain = readRDS("data/NewHillerMainTrees.rds")
oldZoMain = readRDS("data/RemadeTreesAllZoonomiaSpecies.rds")

which(is.na(hillerMain$masterTree$edge.length))
which(is.na(oldZoMain$masterTree$edge.length))



#------ 

mainTrees = readRDS("Data/zoonomiaAllMammalsTrees.rds")
masterTree = mainTrees$masterTree

masterTree$tip.label
