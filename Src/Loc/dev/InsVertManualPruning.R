

mainTrees = readRDS("data/zoonomiaAllMammalsTrees.rds")
phenotypeVector = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalPhenotypeVector.rds")
categoricalTreeFilename = "Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalTree.rds"
length(phenotypeVector)




# ----- Rodents to drop -----
manualAnnotsTrimmed = manualAnnots
which(manualAnnots$ZoonomiaTip %in% names(phenotypeVector))
manualAnnotsTrimmed = manualAnnotsTrimmed[which(manualAnnots$ZoonomiaTip %in% names(phenotypeVector)), ]
#write.csv(manualAnnotsTrimmed, "Results/InsVertPruningTemp.csv")
testTree = readRDS(categoricalTreeFilename)
length(testTree$tip.label)
table(manualAnnotsTrimmed$MSWC_Order)



rodentTree = mainTrees$masterTree
length(rodentTree$tip.label)
prunedTips = rodentTree$tip.label[which(!rodentTree$tip.label %in% testTree$tip.label)]

rodentTree = drop.tip(rodentTree, prunedTips)

testRodent = rodentTree
testRodent$edge.length[1:length(testRodent$edge.length)] = 1

testPruned = testTree
testPruned$edge.length[1:length(testPruned$edge.length)] = 1

all.equal(testPruned, testRodent)

rodentTips = manualAnnotsTrimmed$ZoonomiaTip[manualAnnotsTrimmed$MSWC_Order == "Rodentia"]
nonRodentTips = rodentTree$tip.label[which(!rodentTree$tip.label %in% rodentTips)]

rodentTree = drop.tip(rodentTree, nonRodentTips)


autopruner(rodentTree, dropPercent = 0.5, nameConversionColumn = nameColumn, nameConversionData = spreadSheetLocation, tipsToKeep = manualPruningProtections)
autopruner(rodentTree, dropPercent = 1, nameConversionColumn = nameColumn, nameConversionData = spreadSheetLocation)

rodentsToKeep = c(
  "vs_HLmicAgr2", "voleClade",
  "vs_HLmusAve1", "vs_HLgraMur1", "vs_HLgliGli1", "doorMouseTransition",
  "vs_HLsciCar1", "vs_HLxerIna1", "vs_HLmarFla1", "vs_HLspeDau1", "vs_HLcynGun1", "SquirrelPlusTransitions",
  "vs_mm10", "mouseClade", "commonMouse",
  "vs_HLcteGun1", "vs_HLcavTsc1", "gundiGuineaPigClade"
  
)
rodentsToDrop = c(
  "vs_HLellTal1", "vs_HLellLut1", "vs_HLarvAmp1", "voleClade",
  "vs_HLmusSpi1", "vs_HLmusCar1", "vs_HLmasCou1", "vs_HLmusPah1", "mouseClade",
  "vs_HLhysCri1", "vs_HLthrSwi1", "vs_HLpetTyp1", "vs_hetGla2", "vs_chiLan1", "vs_HLdinBra1", "vs_HLcteSoc1", "vs_octDeg1", "vs_HLcoePre1", "vs_HLdasPun1", "vs_HLdolPat1", "gundiGuineaPigClade"
)

length(rodentsToDrop)
# ---------

# ---- Ruminants to drop ----

prunedTree = testTree
ruminantsPhenTree = extract.clade(prunedTree, which(prunedTree$node.label == "Ruminantia")+length(prunedTree$tip.label))

ruminantsTree = mainTrees$masterTree
irrelevantTips = ruminantsTree$tip.label[which(!ruminantsTree$tip.label %in% ruminantsPhenTree$tip.label)]
ruminantsTree = drop.tip(ruminantsTree, irrelevantTips)

ruminantsTestTree = ruminantsTree
ruminantsTestTree$edge.length[1:length(ruminantsTestTree$edge.length)] = 1
prunedRumTestTree = ruminantsPhenTree
prunedRumTestTree$edge.length[1:length(prunedRumTestTree$edge.length)] = 1
all.equal(ruminantsTestTree, prunedRumTestTree)

plot(ruminantsTree, show.node.label = T)
table(manualAnnotsTrimmed$MSWC_Family)

comRuminantsTree = ZoonomTreeNameToCommon(ruminantsTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)

bovidTree = extract.clade(ruminantsTree, which(ruminantsTree$node.label == "Bovidae")+length(ruminantsTree$tip.label))
comBovidTree = ZoonomTreeNameToCommon(bovidTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)
plot(comBovidTree, show.node.label = T)
plot(ruminantsTree, show.node.label = T)

autopruner(ruminantsTree, dropPercent = 1, nameConversionColumn = nameColumn, nameConversionData = spreadSheetLocation)
autopruner(bovidTree, dropPercent = 1, nameConversionColumn = nameColumn, nameConversionData = spreadSheetLocation)

ruminantsToKeep = c(
  "vs_HLoviNivLyd1", "Caprinae", 
  "vs_HLproPrz1", "Antelopinae",
  "vs_bisBis1", "Bovinae",
  "vs_HLodoHem1", "Cervidae"
)
ruminantsToDrop = c(
  "vs_HLoryGaz1", "vs_HLbeaHun1", "vs_HLkobLecLec1", "vs_HLkobLecLec1", "vs_HLmadKir1", "vs_HLneoPyg1", "vs_HLphiMax1", "vs_HLoreOre1", "vs_HLneoMos1", "vs_HLaepMel1", "vs_HLtraImb1", "Bovidae",
  "vs_HLhydIne1", "vs_HLmunMun1", "Cervidae",
  "vs_HLtraKan1", "mouseDeerOtherIsKept"
)

# -----------------------------  Figuring out how/where to start the ruminant tree ---
plot(testTree, show.node.label = T)
?plotTree()

prunedTree = testTree

prunedTree$node.label
which(prunedTree$node.label == "Boreoeutheria")

trimmedTree = prunedTree

trimmedTree = extract.clade(trimmedTree, 4+length(trimmedTree$tip.label))

comTrimmedTree = ZoonomTreeNameToCommon(trimmedTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)

plot(comTrimmedTree, show.node.label = T)
which(prunedTree$node.label == "Laurasiatheria")

comTrimmedTree = extract.clade(comTrimmedTree, 93+length(trimmedTree$tip.label))
plot(comTrimmedTree, show.node.label = T)

which(comTrimmedTree$node.label == "Cetartiodactyla_-even-toed_ungulates-")
comTrimmedTree = extract.clade(comTrimmedTree, 2+length(comTrimmedTree$tip.label))
plot(comTrimmedTree, show.node.label = T)

comTrimmedTree = extract.clade(comTrimmedTree, 3+length(comTrimmedTree$tip.label))
plot(comTrimmedTree, show.node.label = T)

comTrimmedTree = extract.clade(comTrimmedTree, 15+length(comTrimmedTree$tip.label))
plot(comTrimmedTree, show.node.label = T)



# ------- Other families that need pruning ------

manualAnnotsTrimmed = manualAnnots
which(manualAnnots$ZoonomiaTip %in% names(phenotypeVector))
manualAnnotsTrimmed = manualAnnotsTrimmed[which(manualAnnots$ZoonomiaTip %in% names(phenotypeVector)), ]

familySize = table(manualAnnotsTrimmed$MSWC_Family)
familySize[order(familySize, decreasing = T)]

# --- Checking each outlier family. ----

# - Cricetidae  - 
#ten species, some transitions, any removable species?
#maybe the canyon/cactus mouse. 

#Vespertilionidae   
#NEEDS WORK, no transitions. 

#Cercopithecidae        
#Transitions, is clean

#Mustelidae 
#Transitions, is good. 

#Phyllostomidae           
#Transitions, is good. 

#Muridae      
#Transitions, is good. 

#Pteropodidae       
#NEEDS WORK, no transitions. 

#Delphinidae 
#Single transition, may be excessive, but saying is fine for now. 

#Sciuridae         
#Transitions, is good. 

# ------- Vespertilionidae ------

prunedTree = testTree

prunedTree$node.label
plot(prunedTree, show.node.label = T)
grep("Chirop", prunedTree$node.label)
prunedTree$node.label[142]
batTree = extract.clade(prunedTree, which(prunedTree$node.label == "Chiroptera_-bats-")+length(prunedTree$tip.label))
plot(batTree, show.node.label = T)
grep("Vespert", prunedTree$node.label)
prunedTree$node.label[145]

insectTree = extract.clade(prunedTree, which(prunedTree$node.label == "Vespertilionidae_sensu_lato")+length(prunedTree$tip.label))
plot(insectTree, show.node.label = T)
autopruner(insectTree, dropPercent = 1, nameConversionColumn = nameColumn, nameConversionData = spreadSheetLocation)
autopruner(insectTree, dropPercent = 1, nameConversionColumn = nameColumn, nameConversionData = spreadSheetLocation, tipsToKeep = c("vs_HLmyoLuc1", "vs_eptFus1"))

insectBatsToKeep = c(
  "vs_HLminSch1", "outerVespert", 
  "vs_HLmyoSep1", "Nearctic",
  "vs_HLmyoMyo6", "Myotis",
  "vs_eptFus1", "Vespertilioninae",
)
insectBatsToDrop = c(
  "vs_HLmurAurFea1", "outerVespert",
  "vs_HLmyoLuc1", "Nearctic",
  "vs_myoDav1", "Myotis",
  "vs_HLpipPip1", "vs_HLlasBor1", "vs_HLnycHum2", "Vespertilioninae",
)


# ------- Pteropodidae ------

#Pteropodidae       
#NEEDS WORK, no transitions. 

prunedTree$node.label
plot(prunedTree, show.node.label = T)
fruitTree = prunedTree

grep("Chirop", prunedTree$node.label)
prunedTree$node.label[142]
batTree = extract.clade(prunedTree, which(prunedTree$node.label == "Chiroptera_-bats-")+length(prunedTree$tip.label))
plot(batTree, show.node.label = T)
grep("Yin", prunedTree$node.label)
prunedTree$node.label[164]



fruitTree = extract.clade(prunedTree, which(prunedTree$node.label == "Yinpterochiroptera_-fruit_bats+rhino+hippo+mega-")+length(prunedTree$tip.label))
fruitTree = extract.clade(prunedTree, which(prunedTree$node.label == "Pteropodidae_-fruit_bats-")+length(prunedTree$tip.label))

plot(fruitTree, show.node.label = T)
autopruner(fruitTree, dropPercent = 1, nameConversionColumn = nameColumn, nameConversionData = spreadSheetLocation)


fruitBatsToKeep = c(
    "vs_HLpteGig1", "FoxLongTounge", 
    "vs_HLcynBra1", "outerPeropodidae",
    "vs_HLrouLes1", "Roussetinae",
  )
fruitBatsToDrop = c(
  "vs_HLmacSob1", "FoxLongTounge",
  "vs_HLeidHel2", "outerPeropodidae",
  "vs_HLeonSpe1", "Roussetinae",
)