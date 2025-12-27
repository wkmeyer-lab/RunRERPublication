a = b #prevent full runs
library(RERconverge)
library(tools)

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

collapsedClades$Platypus = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(1)])
collapsedClades$Opossums = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(3,4,5)])
collapsedClades$Koala = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(8,9)])
collapsedClades$Kangaroos = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(10,11,12,13)])
collapsedClades$Anteaters = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(23,24)])
collapsedClades$Sloths = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(25,26)])
collapsedClades$Elephant= MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(19,20,21)])
collapsedClades$Aardvark = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(14,15,16,17,18)])
collapsedClades$Strepsirrhini = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(27,28,29,30,31,32,33,34,35,36,37)])
collapsedClades$Atelidae = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(38,39,40,41,42,43)])
collapsedClades$Chimpanze = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(44:55)])
collapsedClades$Hares = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(56,57)])
collapsedClades$Squirrels = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(58:66)])
collapsedClades$Capybara = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(67:70)])
collapsedClades$Beaver = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(72:74)])
collapsedClades$Jerboa = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(75:77)])
collapsedClades$Deomyinae = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(90:97)])
collapsedClades$Vole = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(86:89)])
collapsedClades$Neotominae = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(80:85)])
collapsedClades$`African Hedgehogs` = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(99,100)])
collapsedClades$`Talpa europaea` = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(101:104)])
collapsedClades$`Flying Fox` = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(105:107)])
collapsedClades$Rhinolophidae = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(108:113)])
collapsedClades$`Big Brown Bat` = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(125:129)])
collapsedClades$Phyllostomidae = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(121:124)])
collapsedClades$Noctilio = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(114)])
collapsedClades$Horse = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(168:170)])
collapsedClades$Pig = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(172:173)])
collapsedClades$`bos bison` = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(194:196)])
collapsedClades$`Humpback Whale` = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(175:178)])
collapsedClades$`Dolphins` = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(182:186)])
collapsedClades$`Pangolin` = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(130:131)])
collapsedClades$`Lion` = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(137:139)])
collapsedClades$`Meerkat` = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(135:136)])
collapsedClades$`Dog` = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(140:141)])
collapsedClades$`Brown Bear` = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(142:144)])
collapsedClades$`Odobenus rosmarus` = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(146:149)])
collapsedClades$`Phocidae` = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(150:153)])
collapsedClades$`Procyon lotor` = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(156:158)])
collapsedClades$`Lontra provocax` = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(162:166)])
collapsedClades$`Tasmanian Devil` = MRCA(commonCategoricalTree, liamNonliamTipConversionIndex[c(6:7)])
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
