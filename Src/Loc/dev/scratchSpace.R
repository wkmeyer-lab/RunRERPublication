library(RERconverge)
library(insight)
library(stringr)

a = b #This line prevents accental full runs
function (realenrich, permvals) 
{
  groups = length(realenrich)
  c = 1
  while (c <= groups) {
    current = realenrich[[c]]
    realenrich[[c]] = current[match(rownames(permvals$enrichStat[[c]]), 
                                    rownames(current)), ]
    c = c + 1
  }
  permenrich = permvals$enrichStat
  enrichpvals = vector("list", length(realenrich))
  groups = length(realenrich)
  count = 1
  while (count <= groups) {
    currreal = realenrich[[count]]
    currenrich = permenrich[[count]]
    rowlen = nrow(currenrich)
    rowcount = 1
    pvallist = c()
    while (rowcount <= rowlen) {
      if (is.na(currreal[rowcount, ]$stat)) {
        pval = lessnum/denom
      }
      else {
        lessnum = sum(abs(currenrich[rowcount, ]) > abs(currreal[rowcount, 
        ]$stat), na.rm = T)
        denom = sum(!is.na(currenrich[rowcount, ]))
        pval = lessnum/denom
      }
      pvallist = c(pvallist, pval)
      rowcount = rowcount + 1
    }
    names(pvallist) = rownames(currreal)
    enrichpvals[[count]] = pvallist
    count = count + 1
  }
  enrichpvals
}



domesticationRERs = readRDS("Output/Domestication/DomesticationRERFile.rds")
domesticationCorrelations = readRDS("Output/Domestication/DomesticationCorrelationFile.rds")
domesticationPaths = readRDS("Output/Domestication/DomesticationPathsFile.rds")

makeTextPlot = function(data, collumn){
  ValueHead = head(data[order(collumn),], n=40)
  ValueHead$N = as.character(ValueHead$N)
  ValueHead$Rho = round(ValueHead$Rho, digits = 5)
  ValueHead$Rho = as.character(ValueHead$Rho)
  ValueHead = format_table(ValueHead, pretty_names = F, digits = "scientific5")
  ValueHead
}

domesticationTopVals = makeTextPlot(domesticationCorrelations, domesticationCorrelations$p.adj)

genesRankedPermP = rownames(domesticationCorrelations[order(domesticationCorrelations$p.adj),])


RERFile = readRDS("Data/CVHRemakeRERFile.rds")
PathsFile = readRDS("Data/CVHRemakePathsFile.rds")
plotRers(domesticationRERs, "AC010255", PathsFile)
for(i in 1:100){
  plotRers(domesticationRERs, genesRankedPermP[i], domesticationPaths)
  message(i)
  Sys.sleep(4)
}


domesticationRERs = readRDS("Output/Domestication/DomesticationRERFile.rds")
domesticationCorrelations = readRDS("Output/Domestication/DomesticationCorrelationFile.rds")
domesticationPaths = readRDS("Output/Domestication/DomesticationPathsFile.rds")

source("Src/Reu/showRERPlots.R")
showRERPlots(domesticationRERs, domesticationPaths, domesticationCorrelations)




library(RERconverge)
source("Src/Reu/ZoonomTreeNameToCommon.R")
currentDomTree = readRDS("Output/Domestication/DomesticationBinaryForegroundTree.rds")
newickTree = read.tree("Output/Domestication/IgnoreGit/DomesticationBinaryTreeNewickFurtherTrim.txt")
plotBinaryTreeViewable(newickTree)

ZoonomTreeNameToCommon(currentDomTree)

ZoonomTreeNameToCommon(newickTree)

plotBinaryTreeViewable(newickTree)
saveRDS(newickTree, "Output/Domestication/DomesticationBinaryForegroundTree.rds")
library(RERconverge)

fixPseudoroot(currentDomTree)

write.tree(currentDomTree, "Output/Domestication/IgnoreGit/DomesticationBinaryTreeNewickPriorWorkingRun.txt")


nobatTree = read.tree("Output/Domestication/IgnoreGit/DomesticationBinaryTreeNewickRemoveBat.txt")
plotBinaryTreeViewable(nobatTree)

saveRDS(nobatTree, "Output/Domestication/DomesticationBinaryForegroundTree.rds")

binaryPhenotypeTree = readRDS(binaryPhenotypeTreeLocation)
pathsObject = tree2Paths(binaryPhenotypeTree, mainTrees, binarize=T, useSpecies = speciesFilter)


fixPseudoroot(nobatTree, mainTrees$masterTree)


domesticationRERs = readRDS("Output/Domestication/DomesticationRERFile.rds")
domesticationCorrelations = readRDS("Output/Domestication/DomesticationCorrelationFile.rds")
domesticationPaths = readRDS("Output/Domestication/DomesticationPathsFile.rds")

source("Src/Reu/showRERPlots.R")
showRERPlots(domesticationRERs, domesticationPaths, domesticationCorrelations)


odlfilter = readRDS("Output/Domestication/DomesticationSpeciesFilter.rds")
domesticationTree = readRDS("Output/Domestication/DomesticationBinaryForegroundTree.rds")
newSpecies = domesticationTree$tip.label

odlfilter %in% newSpecies
newSpecies %in% odlfilter

newfilter = newSpecies
saveRDS(newfilter, "Output/Domestication/DomesticationSpeciesFilter.rds")

enrichSetNames = rownames(enrichmentResult2)

enrichSetNames[grep("GO", enrichSetNames)]

geneSetsInAnnotFile = annotationsList$enrichmentGmtFile$geneset.names

all.equal(enrichSetNames[grep("GO", annotationsList$enrichmentGmtFile$geneset.names)], enrichSetNames[grep("GO", enrichSetNames)])

genesetsInAnnotLIst = enrichSetNames[grep("GO", annotationsList$enrichmentGmtFile$geneset.names)]
genesetInOutput = enrichSetNames[grep("GO", enrichSetNames)]

genesetInOutput %in% geneSetsInAnnotFile

goSetsInAnnot = geneSetsInAnnotFile[grep("GO", geneSetsInAnnotFile)]
goSetsInEnrich = enrichSetNames[grep("GO", enrichSetNames)]

goSetsInEnrich %in% goSetsInAnnot

grep("GO_B", geneSetsInAnnotFile)

correlateWithBinaryPhenotype
getAllCor



corellateTest = function (RERmat, charP, method = "auto", min.sp = 10, min.pos = 2, 
          winsorizeRER = NULL, winsorizetrait = NULL, weighted = F) 
{
  RERna = (apply(is.na(RERmat), 2, all))
  iicharPna = which(is.na(charP))
  method = "p"
  winsorizeRER = 3
  winsorizetrait = 3

  win = function(x, w) {  # this is the function that does the windsorizing; which is to say, capping the top values to the lowest of
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
    names(tables) = rep(NA, n)
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
          yfacts = as.factor(y)
          df = data.frame(x, yfacts)
          colnames(df) = c("RER", "category")
          kres = kruskal.test(RER ~ category, data = df)
          kres_Hval = kres$statistic
          kres_pval = kres$p.value
          effect_size = kres_Hval/(nb - 1)
          corout[i, 1:3] = c(effect_size, nb, kres_pval)
          dunn = dunnTest(RER ~ category, data = df, 
                          method = "bonferroni")
          groups = dunn$res$Comparison
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
            tables[[name]][i, "Rho"] = dunn$res$Z[k]
            tables[[name]][i, "P"] = dunn$res$P.adj[k]
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

correlateWithContinuousPhenotype()

source("Src/Reu/ZoonomTreeNameToCommon.R")

pdf(height = 40, width = 40)
ZoonomTreeNameToCommon(mainTrees$masterTree)
dev.off()

mainTrees = readRDS("Data/RemadeTreesAllZoonomiaSpecies.rds")
RERObject = CVHRERs = readRDS("Output/Old/CVHRemake/CVHRemakeRERFile.rds")
phenotypeTree = readRDS("Output/CVHRemake/CVHRemakeBinaryForegroundTree.rds")
foregroundSpecies = readRDS("Output/Old/CVHRemake/CVHRemakeBinaryTreeForegroundSpecies.rds")
geneOfInterest = "FNDC11"
foregroundName = "Carnivore"
BackgroundName = "Herbivore"
CVHCorrelations = readRDS("Output/CVHRemake/CVHRemakeCorrelationFile.rds")
CVHPermulationP = readRDS("Output/CVHRemake/CVHRemakeCombinedPrunedFastAllPermulationsPValue.rds")
CorrelationData = CVHCorrelations
CorrelationData$permPValue = CVHPermulationP
CVHPaths = readRDS("Output/CVHRemake/CVHRemakePathsFile.rds")

sortedCorrelations = CorrelationData[order(abs(CorrelationData$Rho), decreasing = T),]

library(ggplot)
library(RERconverge)

rerViolinPlot = function(mainTrees, RERObject, phenotypeTree, foregroundSpecies, geneOfInterest, foregroundName = "Foreground", backgroundName = "Background", correlationFile = NULL){
  source("Src/Reu/RERConvergeFunctions.R")
  rerTree = returnRersAsTree(mainTrees, RERObject, geneOfInterest, foregroundSpecies, plot = F)
  relativeRate = rerTree$edge.length
  geneTree = mainTrees$trees[[geneOfInterest]]
  phenotypePath = tree2Paths(phenotypeTree, mainTrees)
  
  edgeIndexMap = edgeIndexRelativeMaster(geneTree, mainTrees$masterTree)
  relevantBranches = mainTrees$matIndex[edgeIndexMap[,c(2,1)]]
  relevantPath = phenotypePath[relevantBranches]
  
  phenotypeVector = c(rep(NA, length(relativeRate)))
  phenotypeVector[relevantPath == 1] = foregroundName
  phenotypeVector[relevantPath == 0] = BackgroundName
  
  #Adding Variables with different caps for better Aesthetics 
  Phenotype = phenotypeVector
  RelativeRate = relativeRate
  
  rateWithPhenotype = data.frame(RelativeRate, Phenotype)
  rateWithPhenotype = rateWithPhenotype[!is.na(rateWithPhenotype$Phenotype),]

  plot = ggplot(rateWithPhenotype, aes(x=Phenotype, y=RelativeRate, col=Phenotype)) + 
    geom_violin(adjust=1/3) +
    geom_jitter(position=position_jitter(0.2)) +
    theme_classic() +
    scale_color_manual(values=c("lightsalmon","darkgreen")) +
    theme(text = element_text(size = 20))+
    ggtitle(geneOfInterest)
  
  if(!is.null(correlationFile)){
    rowNumber = grep(geneOfInterest, rownames(correlationFile))
    if("permPValue" %in% colnames(correlationFile)){
      plot = plot + ggtitle(geneOfInterest, subtitle = paste("p.adj = ", signif(correlationFile$p.adj[rowNumber], 4), "    Permulated P = ", signif(correlationFile$permPValue[rowNumber],4))) 
    }else{
      plot = plot + ggtitle(geneOfInterest, subtitle = paste("p.adj = ", signif(correlationFile$p.adj[rowNumber], 4))) 
    }
  }
  plot
  return(plot)
}

source("Src/Reu/rerViolinPlot.R")

newPlot = rerViolinPlot(mainTrees, CVHRERs, phenotypeTree, foregroundSpecies, "SLC14A2", "Carnivore", "Herbivore", )
newPlot

quickViolin = function(geneInterest){
  plot = rerViolinPlot(mainTrees, CVHRERs, phenotypeTree, foregroundSpecies, geneInterest, "Carnivore", "Herbivore", "blue", "black", CorrelationData)
  print(plot)
}



geneNames = rownames(sortedCorrelations)
i=1

grep("SPECC1L", rownames(sortedCorrelations))

for(i in 23:100){
  quickViolin(geneNames[i])
  print(i)
}

quickViolin("IQANK1")

rerViolinPlot(mainTrees, CVHRERs, phenotypeTree, foregroundSpecies, "PCDH11X", "Carnivore", "Herbivore", "lightsalmon", "darkgreen", CorrelationData)

#----
#----
source("Src/Reu/makeMasterAndGeneTreePlots.R")

for(i in 1:100){
  makeMasterVsGeneTreePlots(mainTrees, CVHRERs, geneNames[i], foregroundSpecies)
  message(i)
  Sys.sleep(2)
}

makeMasterAndGeneTreePlots(mainTrees,"DNAH1", CVHRERs,  foregroundSpecies, correlationPlot = F, tipColumn = "ZoonomiaTip")
  
?ggtitle()

mainTrees$trees[["DNAH1"]]

newPlot = rerViolinPlot(mainTrees, CVHRERs, phenotypeTree, foregroundSpecies, "BCAT2", "Carnivore", "Herbivore", correlData)
newPlot

source("Src/Reu/ZonomNameConvertMatrixCommon.R")

commonRERs = ZonomNameConvertMatrixCommon(CVHRERs)

plotRersNew(commonRERs, "IQANK1", CVHPaths, sort = F)
plotRersNew(commonRERs, "CAPN8", CVHPaths, sort = F)

geneOfInterest= "BCAT2"
correlationFile = sortedCorrelations










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


gmtAnnotations$geneset.names[grep("S_22Q11_DELETION_DN", gmtAnnotations$geneset.names)]
gmtAnnotations$geneset.descriptions[679]
grep("METHIONINE_M", rownames(enrichmentResult2))
enrichmentResult2[2902,]

sortedCorrelations


makeTextPlotSlim = function(data, collumn){
  ValueHead = head(data[order(collumn),], n=40)
  ValueHead$N = as.character(ValueHead$N)
  ValueHead$Rho = round(ValueHead$Rho, digits = 3)
  ValueHead$Rho = as.character(ValueHead$Rho)
  ValueHead = format_table(ValueHead, pretty_names = F, digits = "scientific2")
  ValueHead
}
slimHead = makeTextPlotSlim(sortedCorrelations, sortedCorrelations$p.adj)
textplot(slimHead[1:4], mar = c(0,0,2,0), cmar = 1)
title(main = paste("Top phenotype genes"))

tes = read.delim("../../QuickGO-annotations-1682738614603-20230429.tsv")
tes$SYMBOL
aaCatabolismGenes = unique(tes$SYMBOL)
aacGenesInSet = aaCatabolismGenes[which(aaCatabolismGenes %in% names(mainTrees$trees))]

sortedCorrelations$Rank = 1:nrow(sortedCorrelations)

aaGeneCorels = sortedCorrelations[which(rownames(sortedCorrelations) %in% aacGenesInSet),] 

quickViolin(rownames(aaGeneCorels)[4])



# ---------
rerpath = find.package('RERconverge')
toytreefile = "subsetMammalGeneTrees.txt"
toyTrees=readTrees(paste(rerpath,"/extdata/",toytreefile,sep=""), max.read = 200)
data("sleepPattern")
names(sleepPattern)

cvhNames = readRDS("Output/CVHRemake/CVHRemakeBinaryTreeForegroundSpecies.rds")


cvhEnrich = readRDS("Output/CVHRemake/CVHRemakeEnrichmentFile.rds")


(correlation[order(correlation$p.adj),])


catTree1 = readRDS("Output/CategoricalDiet/CategoricalDietCategoricalTree.rds")
catTree2 = readRDS("Output/CategoricalDiet/CategoricalDietCategoricalTree--possibleTerminal.rds")

all.equal(catTree1, catTree2)

write.nexus(mainTrees, "Results/Maintrees.txt")

class(mainTrees)

mainTrees2 = mainTrees

mainTrees2$trees = mainTrees2$trees[-c(201:16209)]

phenotypeVector2 = phenotypeVector[1:20]

phenotypeHOP = phenotypeVector
phenotypeHGP = phenotypeVector
phenotypeGPI = phenotypeVector
phenotypeCPG = phenotypeVector
phenotypeGPO = phenotypeVector

phenotypeCHO = phenotypeVector


phenotypeHGPI = phenotypeVector
phenotypeCGPI = phenotypeVector

phenotypeHCGP = phenotypeVector
phenotypeCGPO = phenotypeVector

phenotypeGPCO = phenotypeCGPO
phenotypeGPCH = phenotypeHCGP


phenotypeCHOG = phenotypeVector
phenotypeCH_OG = gsub("Generalist", "Omnivore", phenotypeCHOG)

phenotypeAllDataGPCHOI = phenotypeVector


phenotypeVector = phenotypeHGP
phenotypeVector = phenotypeGPI
phenotypeVector = phenotypeCHO_G
phenotypeVector = phenotypeCHOG

permsStartTime = Sys.time()                                                     #get the time before start of permulations
permulationData = categoricalPermulations(mainTrees, phenotypeVector, rm = "SYM", rp = "auto", ntrees = 5)
permsEndTime = Sys.time()                                                       #get time at end of permulations


timeCH_OG = permsEndTime-permsStartTime

timeCH_OG = timeCHO_G

timeGPH = timeHGP
timePHO = timeHOP
phenotypeGPC = phenotypeCPG
phenotypeGPH = phenotypeHGP
phenotypePHO = phenotypeHOP
phenotypeGPCI = phenotypeCGPI

length(phenotypeGPC)
timeGPC
length(phenotypeCHO)
timeCHO
length(phenotypeCH_OG)
timeCH_OG


timePHO
timeGPH
timeGPO
timeGPC
timeGPI


timeGPC
timeGPO
timeGPCO

timeCHO
timeCH_OG
timeCHOG

phenotypeVectors = list(phenotypeGPC, phenotypeGPO, phenotypeGPI, phenotypeGPH, phenotypePHO, phenotypeCHO, phenotypeGPCO, phenotypeGPCH, phenotypeGPCI, phenotypeCHOG, phenotypeCH_OG, phenotypeAllDataGPCHOI)

names(phenotypeVectors) = c("phenotypeGPC", "phenotypeGPO", "phenotypeGPI", "phenotypeGPH", "phenotypePHO", "phenotypeCHO", "phenotypeGPCO", "phenotypeGPCH", "phenotypeGPCI", "phenotypeCHOG", "phenotypeCH_OG", "phenotypeAllDataGPCHOI")

#saveRDS(phenotypeVectors, file = "Results/CategoricalPermuationsTimingPhenotypes.rds")
TimeLists = list(timePHO, timeGPH, timeGPO, timeGPC, timeGPI, timeGPCO, timeCHO, timeCHOG, timeCH_OG)
names(TimeLists) = c("timePHO", "timeGPH", "timeGPO", "timeGPC", "timeGPI", "timeGPCO", "timeCHO", "timeCHOG", "timeCH_OG")
#saveRDS(TimeLists, file = "Results/CategoricalPermulationsTimes.rds")
#saveRDS(mainTrees, "Results/CategoricalPermulationsTimingTrees.rds")

names(phenotypeVectors)

?install.packages()


