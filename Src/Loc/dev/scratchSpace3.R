library(RERconverge)

testData = readRDS("Data/newHillerMainTrees.rds")

plotTree(testData$masterTree)

DomesticationTree = readRDS("Output/Domestication/DomesticationBinaryForegroundTree.rds")
plotTree(DomesticationTree)

source("Src/Reu/ZoonomTreeNameToCommon.R")
ZoonomTreeNameToCommon(DomesticationTree)


hiller4PhenTree = readRDS("Output/NewHiller4Phen/NewHiller4PhenCategoricalTree.rds")
plotTree(hiller4PhenTree)

hillerPhenotypeTable = read.csv("Data/HillerZoonomiaPhenotypeTable.csv")
supraPrimates = hillerPhenotypeTable[which(hillerPhenotypeTable$Supraprimates ==1),]
supraPrimates$FaName


plotTreeHighlightBranches(hiller4PhenTree,hlspecies = supraPrimates$FaName, hlcols = "blue")


all.equal(hillerPhenotypeTable$phenotype, hillerPhenotypeTable$phenotypeSimplified)
hillerPhenotypeTable$phenotype[which(!hillerPhenotypeTable$phenotype == hillerPhenotypeTable$phenotypeSimplified)]

GOValues = enrichmentResultSets[4][[1]]

GOTable = GOValues[order(GOValues$pval),]

GOTable[c(15,35),]   

GSEAValues = enrichmentResultSets[5][[1]]
GSEATable = GSEAValues[order(GSEAValues$pval),]
      
GSEATable[c(8,17,18,23,37),]

AAMetabloismTable = rbind(GOTable[c(15,35),], GSEATable[c(8,17,18,23,37),])

dput(rownames(AAMetabloismTable))
pathwayNames = rownames(AAMetabloismTable)

pathwayNames = gsub("_", " ", pathwayNames)
pathwayNames = str_to_title(pathwayNames)
pathwayNames = gsub("Gobp ", "", pathwayNames)
pathwayNames = gsub(" \\(.*", "", pathwayNames)


rownames(AAMetabloismTable) = pathwayNames

names(AAMetabloismTable) = c("Stat", "P Value", "", "Genes in pathway", "Gene Enrichment Values")

AAMetabloismTableTrimmed = AAMetabloismTable[,c(1,2,4,5)]

AAMetabloismTableTrimmed = AAMetabloismTableTrimmed[order(AAMetabloismTableTrimmed$stat),]

AAMEtabolismTableList = list(AAMetabloismTableTrimmed)
testPlot = makeGeListPlot(AAMEtabolismTableList, enrichmentSortColumn, 40, "Top pathways by permulation", enrichmentOrderDecrease)
testPlot

demoPDF = "Output/CVHRemake/AminoAcidPathways.pdf"
pdf(demoPDF)
pdf(file = demoPDF, width = 15, height = pdfLength)
testPlot
dev.off()



mainTrees = readRDS("Data/RemadeTreesAllZoonomiaSpecies.rds")
RERObject = CVHRERs = readRDS("Output/CVHRemake/CVHRemakeRERFile.rds")
phenotypeTree = readRDS("Output/CVHRemake/CVHRemakeBinaryForegroundTree.rds")
foregroundSpecies = readRDS("Output/CVHRemake/CVHRemakeBinaryTreeForegroundSpecies.rds")
CVHCorrelations = readRDS("Output/CVHRemake/CVHRemakeCorrelationFile.rds")
CVHPermulationP = readRDS("Output/CVHRemake/CVHRemakeCombinedPrunedFastAllPermulationsPValue.rds")
CorrelationData = CVHCorrelations
CorrelationData$permPValue = CVHPermulationP
source("Src/Reu/rerViolinPlot.R")
quickViolin = function(geneInterest){
  plot = rerViolinPlot(mainTrees, CVHRERs, phenotypeTree, foregroundSpecies, geneInterest, "Carnivore", "Herbivore", "salmon", "darkgreen", CorrelationData)
  print(plot)
}

quickViolin("SDS")


mainTrees$masterTree


newDemoEnrichment = 
  


RERObject = readRDS("Output/NewHil")

plotRers(rermat = RERObject, index = "OR10J5", phenv = pathsObject)

phenVec = readRDS("OUtput/NewHiller4Phen/NewHIller4PhenCategoricalPhenotypeVector.rds")

solPhenVec = c(hg38 = "_Omnivore", panTro5 = "_Omnivore", panPan2 = "_Omnivore",
               gorGor5 = "Herbivore", ponAbe2 = "_Omnivore", nomLeu3 = "Herbivore",
               rheMac8 = "_Omnivore", macFas5 = "_Omnivore", macNem1 = "_Omnivore",
               papAnu3 = "_Omnivore", manLeu1 = "_Omnivore", cerAty1 = "_Omnivore",
               chlSab2 = "_Omnivore", nasLar1 = "Herbivore", rhiRox1 = "Herbivore",
               rhiBie1 = "Herbivore", colAng1 = "Herbivore", HLpilTep1 = "Herbivore",
               calJac3 = "_Omnivore", aotNan1 = "_Omnivore", saiBol1 = "_Omnivore",
               cebCap1 = "_Omnivore", tarSyr2 = "Insectivore", otoGar3 = "Herbivore",
               micMur3 = "_Omnivore", proCoq1 = "Herbivore", galVar1 = "Herbivore",
               jacJac1 = "Herbivore", micOch1 = "Herbivore", criGri1 = "Herbivore",
               mesAur1 = "_Omnivore", perManBai1 = "_Omnivore", HLmusCar1 = "_Omnivore",
               HLmusPah1 = "_Omnivore", rn6 = "_Omnivore", HLmerUng1 = "Herbivore",
               nanGal1 = "Herbivore", HLcasCan1 = "Herbivore", HLdipOrd2 = "_Omnivore",
               hetGla2 = "Herbivore", HLfukDam1 = "_Omnivore", cavPor3 = "Herbivore",
               chiLan1 = "Herbivore", octDeg1 = "Herbivore", speTri2 = "_Omnivore",
               HLmarMar1 = "_Omnivore", oryCun2 = "Herbivore", ochPri3 = "Herbivore",
               vicPac2 = "Herbivore", HLcamFer2 = "_Omnivore", HLcamBac1 = "_Omnivore",
               HLcamDro1 = "_Omnivore", HLturTru3 = "Carnivore", HLorcOrc1 = "Carnivore",
               HLdelLeu1 = "Carnivore", lipVex1 = "Carnivore", phyCat1 = "Carnivore",
               bosTau8 = "Herbivore", HLbosInd1 = "Herbivore", bisBis1 = "Herbivore",
               bosMut1 = "Herbivore", bubBub1 = "Herbivore", HLoviAri4 = "Herbivore",
               HLoviCan1 = "Herbivore", HLcapHir2 = "Herbivore", panHod1 = "Herbivore",
               HLodoVir1 = "Herbivore", HLcerEla1 = "Herbivore", susScr11 = "_Omnivore",
               HLequCab3 = "Herbivore", equPrz1 = "Herbivore", HLequAsi1 = "Herbivore",
               cerSim1 = "Herbivore", felCat8 = "Carnivore", HLaciJub1 = "Carnivore",
               panTig1 = "Carnivore", HLpanPar1 = "Carnivore", canFam3 = "Carnivore",
               HLlycPic1 = "Carnivore", musFur1 = "Carnivore", HLailFul1 = "_Omnivore",
               ailMel1 = "_Omnivore", ursMar1 = "Carnivore", odoRosDiv1 = "Carnivore",
               lepWed1 = "Carnivore", neoSch1 = "Carnivore", manPen1 = "Insectivore",
               HLmanJav1 = "Insectivore", pteAle1 = "Herbivore", HLpteVam2 = "Herbivore",
               rouAeg1 = "Herbivore", HLrhiSin1 = "Insectivore", HLhipArm1 = "Insectivore",
               eptFus1 = "Insectivore", myoDav1 = "Insectivore", myoBra1 = "Insectivore",
               myoLuc2 = "Insectivore", HLminNat1 = "Insectivore", eriEur2 = "_Omnivore",
               sorAra2 = "_Omnivore", conCri1 = "Carnivore", loxAfr3 = "Herbivore",
               triMan1 = "Herbivore", HLproCap2 = "Herbivore", chrAsi1 = "Insectivore",
               eleEdw1 = "Insectivore", oryAfe1 = "Insectivore", dasNov3 = "Insectivore",
               HLchoHof2 = "Herbivore", monDom5 = "_Omnivore", sarHar1 = "_Omnivore",
               HLphaCin1 = "Herbivore", ornAna2 = "Insectivore")


all.equal(phenVec, solPhenVec)
which(!phenVec == solPhenVec)
phenVec[82]
solPhenVec[82]
source("Src/Reu/ZonomNameConvertVector.R")
commonPhenVec = ZonomNameConvertVectorCommon(solPhenVec, manualAnnotLocation = "Data/HillerZoonomiaPhenotypeTable.csv")
names(commonPhenVec) = ZonomNameConvertVectorCommon(names(solPhenVec), manualAnnotLocation = "Data/HillerZoonomiaPhenotypeTable.csv")
commonPhenVec[82]


rownames(enrichmentReorder)

phenotypeVector = readRDS("Output/IPCRelaxTest/HillerIPCPhenotype.rds")
speciesFilter = names(phenotypeVector) 
speciesFilter = speciesFilter[-which(speciesFilter == "ornAna2")]


phenotypeVector = readRDS("Output/HMGRelaxTest/HillerHGMPhenotype.rds")
speciesFilter = names(phenotypeVector) 


library(RERconverge)
mainTrees = readRDS("Data/RemadeTreesAllZoonomiaSpecies.rds")
masterTree = mainTrees$masterTree
write.tree(masterTree, "Results/ZoonomiaMasterTree.tree")


list.files("data")

masterTree2 =  ZoonomTreeNameToCommon(masterTree, scientific = T)
write.tree(masterTree2, "Results/zoonomiaMasterTreeScientific.tree")


source("src/reu/ZoonomTreeNameToCommon.R")

Phen5Tree = readRDS("Output/CategoricalDiet5Phen/CategoricalDiet5PhenCategoricalTree.rds")

ZoonomTreeNameToCommon(Phen5Tree)

hillerConversionTableLocation = "Data/HillerZoonomiaPhenotypeTable.csv"
hillerConversionTable = read.csv("Data/HIllerZoonomPhenotypeTable.csv")

relevantSpecies

hillerNames = match(speciesNames, hillerConversionTable$Zoonomia)
?match()

newHillerMasterTree = read.tree("Data/newHillerMasterTree.txt")
plotTree(newHillerMasterTree)
nodelabels(col = "red", adj = c(0, -0), frame = "none")
tiplabels(col = "blue", frame = "none")
edgelabels(col = "darkgreen", frame = "none")


is.binary(newHillerMasterTree)

newHillerTrees = readTrees('Data/newHillerMainTrees2.txt')


testTree = newHillerTrees$trees[[1]]

is.binary(newHillerMasterTree)

is.binary(testTree)
plotTree(testTree)

demoTree= testTree
demoTree$edge.length = rep(1, length = length(demoTree$edge))

plotTree(demoTree)
ZoonomTreeNameToCommon(demoTree, manualAnnotLocation = hillerConversionTableLocation)


lapply(newHillerTrees$trees[[]], is.binary())


oldHilTrees = oldHillerTrees$trees


all(sapply(oldHilTrees, is.binary))


?readTrees

saveRDS(newHillerTrees, "Data/newHillerMainTrees.rds")

oldHillerTrees = readRDS("Data/mam120aa_trees.rds")

newHillerTrees = names(newHillerTrees$trees)
oldHillerTrees = names(oldHillerTrees$trees)

missingGenes = oldHillerTrees[which(oldHillerTrees %in% newHillerTrees)]


testSubfolderData = readRDS("Output/CategoricalDiet3Phen/Overall/CategoricalDiet3PhenOverallPermulationsCorrelations.rds")


?plotRers

plotRers(RERObject, "ALOX15", phenv = pathsObject)

# --------------------------

length(phenotypeVector)
length(mainTrees$masterTree$tip.label)

removedSpecies = mainTrees$masterTree$tip.label[which(!mainTrees$masterTree$tip.label %in% speciesNames)]


droppedspecies = manualAnnots[[annotColumn]][which(manualAnnots$FaName %in% removedSpecies)]
names(droppedspecies) = manualAnnots$FaName[which(manualAnnots$FaName %in% removedSpecies)]
droppedspecies

which(manualAnnots$FaName %in% removedSpecies)



correlation[order(correlation$p.adj),]


?getPermsBinary
getPermsBinary
function (trees, mastertree, root_sp, fg_vec, sisters_list = NULL, 
          pathvec, plotTreeBool = F) 
{
  tip.labels = mastertree$tip.label
  res = getForegroundInfoClades(fg_vec, sisters_list, trees, 
                                plotTree = F, useSpecies = tip.labels)
  fg_tree = res$tree
  fg.table = res$fg.sisters.table
  t = root.phylo(trees$masterTree, root_sp, resolve.root = T)
  rm = ratematrix(t, pathvec)
  if (!is.null(sisters_list)) {
    fg_tree_info = getBinaryPermulationInputsFromTree(fg_tree)
    num_tip_sisters_true = unlist(fg_tree_info$sisters_list)
    num_tip_sisters_true = num_tip_sisters_true[which(num_tip_sisters_true %in% 
                                                        tip.labels)]
    num_tip_sisters_true = length(num_tip_sisters_true)
    fg_tree_depth_order = getDepthOrder(fg_tree)
  }
  else {
    fg_tree_depth_order = NULL
  }
  fgnum = length(which(fg_tree$edge.length == 1))
  if (!is.null(sisters_list)) {
    internal = nrow(fg.table)
  }
  else {
    internal = 0
  }
  tips = fgnum - internal
  testcondition = FALSE
  while (!testcondition) {
    blsum = 0
    while (blsum != fgnum) {
      sims = sim.char(t, rm, nsim = 1)
      nam = rownames(sims)
      s = as.data.frame(sims)
      simulatedvec = s[, 1]
      names(simulatedvec) = nam
      top = names(sort(simulatedvec, decreasing = TRUE))[1:tips]
      t_iter = foreground2Tree(top, trees, clade = "all", 
                               plotTree = F)
      blsum = sum(t_iter$edge.length)
    }
    t_info = getBinaryPermulationInputsFromTree(t_iter)
    if (!is.null(sisters_list)) {
      num_tip_sisters_fake = unlist(t_info$sisters_list)
      num_tip_sisters_fake = num_tip_sisters_fake[which(num_tip_sisters_fake %in% 
                                                          tip.labels)]
      num_tip_sisters_fake = length(num_tip_sisters_fake)
      t_depth_order = getDepthOrder(t_iter)
      testcondition = setequal(sort(t_depth_order), sort(fg_tree_depth_order)) && 
        (num_tip_sisters_fake == num_tip_sisters_true)
    }
    else {
      t_depth_order = getDepthOrder(t_iter)
      testcondition = setequal(sort(t_depth_order), sort(fg_tree_depth_order))
    }
  }
  if (plotTreeBool) {
    plot(t_iter)
  }
  return(t_iter)
}

function (tree, numperms, trees, root_sp, fg_vec, sisters_list, 
          pathvec, permmode = "cc") 
{
  if (permmode == "cc") {
    tree_rep = lapply(1:numperms, rep_tree, tree = trees)
    permulated.binphens = lapply(tree_rep, simBinPhenoCC, 
                                 mastertree = trees$masterTree, root_sp = root_sp, 
                                 fg_vec = fg_vec, sisters_list = sisters_list, pathvec = pathvec, 
                                 plotTreeBool = F)
  }
  else if (permmode == "ssm") {
    tree_rep = lapply(1:numperms, rep_tree, tree = tree)
    permulated.binphens = lapply(tree_rep, simBinPhenoSSM, 
                                 trees = trees, root_sp = root_sp, fg_vec = fg_vec, 
                                 sisters_list = sisters_list, pathvec = pathvec)
  }
  else {
    stop("Invalid binary permulation mode.")
  }
  output.list <- list()
  output.list[[1]] <- permulated.binphens
  return(output.list)
}


function (numperms, fg_vec, sisters_list, root_sp, RERmat, trees, 
          mastertree, permmode = "cc", method = "k", min.pos = 2, trees_list = NULL, 
          calculateenrich = F, annotlist = NULL) 
{
  pathvec = foreground2Paths(fg_vec, trees, clade = "all", 
                             plotTree = F)
  col_labels = colnames(trees$paths)
  names(pathvec) = col_labels
  if (permmode == "cc") {
    print("Running CC permulation")
    print("Generating permulated trees")
    permulated.binphens = generatePermulatedBinPhen(trees$masterTree, 
                                                    numperms, trees, root_sp, fg_vec, sisters_list, pathvec, 
                                                    permmode = "cc")
    permulated.fg = mapply(getForegroundsFromBinaryTree, 
                           permulated.binphens[[1]])
    permulated.fg.list = as.list(data.frame(permulated.fg))
    phenvec.table = mapply(foreground2Paths, permulated.fg.list, 
                           MoreArgs = list(treesObj = trees, clade = "all"))
    phenvec.list = lapply(seq_len(ncol(phenvec.table)), function(i) phenvec.table[, 
                                                                                  i])
    print("Calculating correlations")
    corMatList = lapply(phenvec.list, correlateWithBinaryPhenotype, 
                        RERmat = RERmat)
    permPvals = data.frame(matrix(ncol = numperms, nrow = nrow(RERmat)))
    rownames(permPvals) = rownames(RERmat)
    permRhovals = data.frame(matrix(ncol = numperms, nrow = nrow(RERmat)))
    rownames(permRhovals) = rownames(RERmat)
    permStatvals = data.frame(matrix(ncol = numperms, nrow = nrow(RERmat)))
    rownames(permStatvals) = rownames(RERmat)
    for (i in 1:length(corMatList)) {
      permPvals[, i] = corMatList[[i]]$P
      permRhovals[, i] = corMatList[[i]]$Rho
      permStatvals[, i] = sign(corMatList[[i]]$Rho) * -log10(corMatList[[i]]$P)
    }
  }
  else if (permmode == "ssm") {
    print("Running SSM permulation")
    if (is.null(trees_list)) {
      trees_list = trees$trees
    }
    RERmat = RERmat[match(names(trees_list), rownames(RERmat)), 
    ]
    print("Generating permulated trees")
    permulated.binphens = generatePermulatedBinPhenSSMBatched(trees_list, 
                                                              numperms, trees, root_sp, fg_vec, sisters_list, pathvec)
    df.list = lapply(trees_list, getSpeciesMembershipStats, 
                     masterTree = mastertree, foregrounds = fg_vec)
    df.converted = data.frame(matrix(unlist(df.list), nrow = length(df.list), 
                                     byrow = T), stringsAsFactors = FALSE)
    attr = attributes(df.list[[1]])
    col_names = attr$names
    attr2 = attributes(df.list)
    row_names = attr2$names
    colnames(df.converted) = col_names
    rownames(df.converted) = row_names
    df.converted$num.fg = as.integer(df.converted$num.fg)
    df.converted$num.spec = as.integer(df.converted$num.spec)
    spec.members = df.converted$spec.members
    grouped.trees = groupTrees(spec.members)
    ind.unique.trees = grouped.trees$ind.unique.trees
    ind.unique.trees = unlist(ind.unique.trees)
    ind.tree.groups = grouped.trees$ind.tree.groups
    unique.trees = trees_list[ind.unique.trees]
    unique.map.list = mapply(matchAllNodesClades, unique.trees, 
                             MoreArgs = list(treesObj = trees))
    unique.permulated.binphens = permulated.binphens[ind.unique.trees]
    unique.permulated.paths = calculatePermulatedPaths_apply(unique.permulated.binphens, 
                                                             unique.map.list, trees)
    permulated.paths = vector("list", length = length(trees_list))
    for (j in 1:length(permulated.paths)) {
      permulated.paths[[j]] = vector("list", length = numperms)
    }
    for (i in 1:length(unique.permulated.paths)) {
      ind.unique.tree = ind.unique.trees[i]
      ind.tree.group = ind.tree.groups[[i]]
      unique.path = unique.permulated.paths[[i]]
      for (k in 1:length(ind.tree.group)) {
        permulated.paths[[ind.tree.group[k]]] = unique.path
      }
    }
    attributes(permulated.paths)$names = row_names
    print("Calculating correlations")
    RERmat.list = lapply(seq_len(nrow(RERmat[])), function(i) RERmat[i, 
    ])
    corMatList = mapply(calculateCorPermuted, permulated.paths, 
                        RERmat.list)
    permPvals = extractCorResults(corMatList, numperms, mode = "P")
    rownames(permPvals) = names(trees_list)
    permRhovals = extractCorResults(corMatList, numperms, 
                                    mode = "Rho")
    rownames(permRhovals) = names(trees_list)
    permStatvals = sign(permRhovals) * -log10(permPvals)
    rownames(permStatvals) = names(trees_list)
  }
  else {
    stop("Invalid binary permulation mode.")
  }
  if (calculateenrich) {
    realFgtree = foreground2TreeClades(fg_vec, sisters_list, 
                                       trees, plotTree = F)
    realpaths = tree2PathsClades(realFgtree, trees)
    realresults = getAllCor(RERmat, realpaths, method = method, 
                            min.pos = min.pos)
    realstat = sign(realresults$Rho) * -log10(realresults$P)
    names(realstat) = rownames(RERmat)
    realenrich = fastwilcoxGMTall(na.omit(realstat), annotlist, 
                                  outputGeneVals = F)
    groups = length(realenrich)
    c = 1
    while (c <= groups) {
      current = realenrich[[c]]
      realenrich[[c]] = current[order(rownames(current)), 
      ]
      c = c + 1
    }
    permenrichP = vector("list", length(realenrich))
    permenrichStat = vector("list", length(realenrich))
    c = 1
    while (c <= length(realenrich)) {
      newdf = data.frame(matrix(ncol = numperms, nrow = nrow(realenrich[[c]])))
      rownames(newdf) = rownames(realenrich[[c]])
      permenrichP[[c]] = newdf
      permenrichStat[[c]] = newdf
      c = c + 1
    }
    counter = 1
    while (counter <= numperms) {
      stat = permStatvals[, counter]
      names(stat) = rownames(RERmat)
      enrich = fastwilcoxGMTall(na.omit(stat), annotlist, 
                                outputGeneVals = F)
      groups = length(enrich)
      c = 1
      while (c <= groups) {
        current = enrich[[c]]
        enrich[[c]] = current[order(rownames(current)), 
        ]
        enrich[[c]] = enrich[[c]][match(rownames(permenrichP[[c]]), 
                                        rownames(enrich[[c]])), ]
        permenrichP[[c]][, counter] = enrich[[c]]$pval
        permenrichStat[[c]][, counter] = enrich[[c]]$stat
        c = c + 1
      }
      counter = counter + 1
    }
  }
  if (calculateenrich) {
    data = vector("list", 5)
    data[[1]] = permPvals
    data[[2]] = permRhovals
    data[[3]] = permStatvals
    data[[4]] = permenrichP
    data[[5]] = permenrichStat
    names(data) = c("corP", "corRho", "corStat", "enrichP", 
                    "enrichStat")
  }
  else {
    data = vector("list", 3)
    data[[1]] = permPvals
    data[[2]] = permRhovals
    data[[3]] = permStatvals
    names(data) = c("corP", "corRho", "corStat")
  }
  data
}
<bytecode: 0x00000174ff39f338>
  <environment: namespace:RERconverge>
  > generatePermulatedBinPhen
function (tree, numperms, trees, root_sp, fg_vec, sisters_list, 
          pathvec, permmode = "cc") 
{
  if (permmode == "cc") {
    tree_rep = lapply(1:numperms, rep_tree, tree = trees)
    permulated.binphens = lapply(tree_rep, simBinPhenoCC, 
                                 mastertree = trees$masterTree, root_sp = root_sp, 
                                 fg_vec = fg_vec, sisters_list = sisters_list, pathvec = pathvec, 
                                 plotTreeBool = F)
  }
  else if (permmode == "ssm") {
    tree_rep = lapply(1:numperms, rep_tree, tree = tree)
    permulated.binphens = lapply(tree_rep, simBinPhenoSSM, 
                                 trees = trees, root_sp = root_sp, fg_vec = fg_vec, 
                                 sisters_list = sisters_list, pathvec = pathvec)
  }
  else {
    stop("Invalid binary permulation mode.")
  }
  output.list <- list()
  output.list[[1]] <- permulated.binphens
  return(output.list)
}


function (numperms, fg_vec, sisters_list, root_sp, RERmat, trees, 
          mastertree, permmode = "cc", method = "k", min.pos = 2, trees_list = NULL, 
          calculateenrich = F, annotlist = NULL) 
{
  pathvec = foreground2Paths(fg_vec, trees, clade = "all", 
                             plotTree = F)
  col_labels = colnames(trees$paths)
  names(pathvec) = col_labels
  if (permmode == "cc") {
    print("Running CC permulation")
    print("Generating permulated trees")
    permulated.binphens = generatePermulatedBinPhen(trees$masterTree, 
                                                    numperms, trees, root_sp, fg_vec, sisters_list, pathvec, 
                                                    permmode = "cc")
    permulated.fg = mapply(getForegroundsFromBinaryTree, 
                           permulated.binphens[[1]])
    permulated.fg.list = as.list(data.frame(permulated.fg))
    phenvec.table = mapply(foreground2Paths, permulated.fg.list, 
                           MoreArgs = list(treesObj = trees, clade = "all"))
    phenvec.list = lapply(seq_len(ncol(phenvec.table)), function(i) phenvec.table[, 
                                                                                  i])
    print("Calculating correlations")
    corMatList = lapply(phenvec.list, correlateWithBinaryPhenotype, 
                        RERmat = RERmat)
    permPvals = data.frame(matrix(ncol = numperms, nrow = nrow(RERmat)))
    rownames(permPvals) = rownames(RERmat)
    permRhovals = data.frame(matrix(ncol = numperms, nrow = nrow(RERmat)))
    rownames(permRhovals) = rownames(RERmat)
    permStatvals = data.frame(matrix(ncol = numperms, nrow = nrow(RERmat)))
    rownames(permStatvals) = rownames(RERmat)
    for (i in 1:length(corMatList)) {
      permPvals[, i] = corMatList[[i]]$P
      permRhovals[, i] = corMatList[[i]]$Rho
      permStatvals[, i] = sign(corMatList[[i]]$Rho) * -log10(corMatList[[i]]$P)
    }
  }
  else if (permmode == "ssm") {
    print("Running SSM permulation")
    if (is.null(trees_list)) {
      trees_list = trees$trees
    }
    RERmat = RERmat[match(names(trees_list), rownames(RERmat)), 
    ]
    print("Generating permulated trees")
    permulated.binphens = generatePermulatedBinPhenSSMBatched(trees_list, 
                                                              numperms, trees, root_sp, fg_vec, sisters_list, pathvec)
    df.list = lapply(trees_list, getSpeciesMembershipStats, 
                     masterTree = mastertree, foregrounds = fg_vec)
    df.converted = data.frame(matrix(unlist(df.list), nrow = length(df.list), 
                                     byrow = T), stringsAsFactors = FALSE)
    attr = attributes(df.list[[1]])
    col_names = attr$names
    attr2 = attributes(df.list)
    row_names = attr2$names
    colnames(df.converted) = col_names
    rownames(df.converted) = row_names
    df.converted$num.fg = as.integer(df.converted$num.fg)
    df.converted$num.spec = as.integer(df.converted$num.spec)
    spec.members = df.converted$spec.members
    grouped.trees = groupTrees(spec.members)
    ind.unique.trees = grouped.trees$ind.unique.trees
    ind.unique.trees = unlist(ind.unique.trees)
    ind.tree.groups = grouped.trees$ind.tree.groups
    unique.trees = trees_list[ind.unique.trees]
    unique.map.list = mapply(matchAllNodesClades, unique.trees, 
                             MoreArgs = list(treesObj = trees))
    unique.permulated.binphens = permulated.binphens[ind.unique.trees]
    unique.permulated.paths = calculatePermulatedPaths_apply(unique.permulated.binphens, 
                                                             unique.map.list, trees)
    permulated.paths = vector("list", length = length(trees_list))
    for (j in 1:length(permulated.paths)) {
      permulated.paths[[j]] = vector("list", length = numperms)
    }
    for (i in 1:length(unique.permulated.paths)) {
      ind.unique.tree = ind.unique.trees[i]
      ind.tree.group = ind.tree.groups[[i]]
      unique.path = unique.permulated.paths[[i]]
      for (k in 1:length(ind.tree.group)) {
        permulated.paths[[ind.tree.group[k]]] = unique.path
      }
    }
    attributes(permulated.paths)$names = row_names
    print("Calculating correlations")
    RERmat.list = lapply(seq_len(nrow(RERmat[])), function(i) RERmat[i, 
    ])
    corMatList = mapply(calculateCorPermuted, permulated.paths, 
                        RERmat.list)
    permPvals = extractCorResults(corMatList, numperms, mode = "P")
    rownames(permPvals) = names(trees_list)
    permRhovals = extractCorResults(corMatList, numperms, 
                                    mode = "Rho")
    rownames(permRhovals) = names(trees_list)
    permStatvals = sign(permRhovals) * -log10(permPvals)
    rownames(permStatvals) = names(trees_list)
  }
  else {
    stop("Invalid binary permulation mode.")
  }
  if (calculateenrich) {
    realFgtree = foreground2TreeClades(fg_vec, sisters_list, 
                                       trees, plotTree = F)
    realpaths = tree2PathsClades(realFgtree, trees)
    realresults = getAllCor(RERmat, realpaths, method = method, 
                            min.pos = min.pos)
    realstat = sign(realresults$Rho) * -log10(realresults$P)
    names(realstat) = rownames(RERmat)
    realenrich = fastwilcoxGMTall(na.omit(realstat), annotlist, 
                                  outputGeneVals = F)
    groups = length(realenrich)
    c = 1
    while (c <= groups) {
      current = realenrich[[c]]
      realenrich[[c]] = current[order(rownames(current)), 
      ]
      c = c + 1
    }
    permenrichP = vector("list", length(realenrich))
    permenrichStat = vector("list", length(realenrich))
    c = 1
    while (c <= length(realenrich)) {
      newdf = data.frame(matrix(ncol = numperms, nrow = nrow(realenrich[[c]])))
      rownames(newdf) = rownames(realenrich[[c]])
      permenrichP[[c]] = newdf
      permenrichStat[[c]] = newdf
      c = c + 1
    }
    counter = 1
    while (counter <= numperms) {
      stat = permStatvals[, counter]
      names(stat) = rownames(RERmat)
      enrich = fastwilcoxGMTall(na.omit(stat), annotlist, 
                                outputGeneVals = F)
      groups = length(enrich)
      c = 1
      while (c <= groups) {
        current = enrich[[c]]
        enrich[[c]] = current[order(rownames(current)), 
        ]
        enrich[[c]] = enrich[[c]][match(rownames(permenrichP[[c]]), 
                                        rownames(enrich[[c]])), ]
        permenrichP[[c]][, counter] = enrich[[c]]$pval
        permenrichStat[[c]][, counter] = enrich[[c]]$stat
        c = c + 1
      }
      counter = counter + 1
    }
  }
  if (calculateenrich) {
    data = vector("list", 5)
    data[[1]] = permPvals
    data[[2]] = permRhovals
    data[[3]] = permStatvals
    data[[4]] = permenrichP
    data[[5]] = permenrichStat
    names(data) = c("corP", "corRho", "corStat", "enrichP", 
                    "enrichStat")
  }
  else {
    data = vector("list", 3)
    data[[1]] = permPvals
    data[[2]] = permRhovals
    data[[3]] = permStatvals
    names(data) = c("corP", "corRho", "corStat")
  }
  data
}


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
