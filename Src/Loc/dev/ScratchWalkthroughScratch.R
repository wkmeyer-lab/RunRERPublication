mainMergedData = read.csv("Data/MergedData.csv")

mainMergedData$ScientificNameFull[mainMergedData$CommonName %in% demoTrees$masterTree$tip.label]

demoTrees$masterTree$tip.label[which(!demoTrees$masterTree$tip.label %in% mainMergedData$CommonName)]


testTipLabel = demoTrees$masterTree$tip.label
testTipLabel = sub("_", " ", testTipLabel)
testTipLabel = sub("_", " ", testTipLabel)
testTipLabel[which(!testTipLabel %in% mainMergedData$CommonName)]
testTipLabel[which(testTipLabel == "Star nosed mole")] = "Star-nosed mole"
testTipLabel[which(testTipLabel == "Elephant shrew")] = "Cape elephant shrew"
testTipLabel[which(testTipLabel == "Seal")] = "Grey seal"
testTipLabel[which(testTipLabel == "Walrus")] = "Pacific walrus"
testTipLabel[which(testTipLabel == "Flying fox")] = "Indian flying fox"
testTipLabel[which(testTipLabel == "Seal")] = "Grey seal"
testTipLabel[which(testTipLabel == "Brown bat")] = "Big brown bat"
testTipLabel[which(testTipLabel == "Myotis bat")] = "David's myotis bat"
testTipLabel[which(testTipLabel == "Rhinoceros")] = "Black Rhinoceros"
testTipLabel[which(testTipLabel == "Bactrian camel")] = "Bactrian Camel"
testTipLabel[which(testTipLabel == "Sheep")] = "Domestic sheep"
testTipLabel[which(testTipLabel == "Goat")] = "Domestic goat"
testTipLabel[which(testTipLabel == "Jerboa")] = "Gobi jerboa"
testTipLabel[which(testTipLabel == "Mouse")] = "Common mouse"
testTipLabel[which(testTipLabel == "Vole")] = "Common vole"
testTipLabel[which(testTipLabel == "Naked mole rat")] = "Naked mole-rat"
testTipLabel[which(testTipLabel == "Guinea pig")] = "Domestic guinea pig"
testTipLabel[which(testTipLabel == "Brush tailed rat")] = "Brush-tailed rat"
testTipLabel[which(testTipLabel == "Baboon")] = "Olive Baboon"
testTipLabel[which(testTipLabel == "Rhesus macaque")] = "Rhesus Macaca"
testTipLabel[which(testTipLabel == "Crab eating macaque")] = "Crab-eating macaque"


testTipLabel[which(!testTipLabel %in% mainMergedData$CommonName)]


scientificNames = mainMergedData$ScientificNameFull[match(testTipLabel, mainMergedData$CommonName)]


dput(scientificNames)


commonNames = mainMergedData$CommonName[match(testTipLabel, mainMergedData$CommonName)]
dput(commonNames)

length(demoTrees$masterTree$tip.label)


dietPhen = mainMergedData$MeyerTrophicLevel[match(testTipLabel, mainMergedData$CommonName)]

dput(dietPhen)

demoTrees = toyTrees
