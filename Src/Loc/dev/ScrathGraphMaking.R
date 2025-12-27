# -----------------
library(ggplot2)
library(dplyr)
library(RERconverge)
library(data.table)
library(gridExtra)
source("Src/Reu/ZoonomTreeNameToCommon.R")
a=b

# -------- 
# Figure 2 for presentation at carnivroy group 

# Load data
InsVertivoreData = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCorrelationFile.rds")
InsVertivoreDataPairwise = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreePairwiseCorrelationFile.rds")

InsVertivoreData = InsVertivoreData
InsVertvioreHvI = InsVertivoreDataPairwise$`Herbivore - Insectivore`
InsVertivoreHvO = InsVertivoreDataPairwise$`Herbivore - Omnivore`
InsVertivoreIvO = InsVertivoreDataPairwise$`Insectivore - Omnivore`
InsVertivoreHvV = InsVertivoreDataPairwise$`Herbivore - Vertivore`
InsVertivoreIvV = InsVertivoreDataPairwise$`Insectivore - Vertivore`
InsVertivoreOvV = InsVertivoreDataPairwise$`Omnivore - Vertivore`

InsVertivoreDataOrdered = InsVertivoreData[order(InsVertivoreData$p.adj),]
InsVertvioreHvIOrdered = InsVertivoreDataPairwise$`Herbivore - Insectivore`[order(InsVertivoreDataPairwise$`Herbivore - Insectivore`$p.adj),]
InsVertivoreHvOOrdered = InsVertivoreDataPairwise$`Herbivore - Omnivore`[order(InsVertivoreDataPairwise$`Herbivore - Omnivore`$p.adj),]
InsVertivoreIvOOrdered = InsVertivoreDataPairwise$`Insectivore - Omnivore`[order(InsVertivoreDataPairwise$`Insectivore - Omnivore`$p.adj),]
InsVertivoreHvVOrdered = InsVertivoreDataPairwise$`Herbivore - Vertivore`[order(InsVertivoreDataPairwise$`Herbivore - Vertivore`$p.adj),]
InsVertivoreIvVOrdered = InsVertivoreDataPairwise$`Insectivore - Vertivore`[order(InsVertivoreDataPairwise$`Insectivore - Vertivore`$p.adj),]
InsVertivoreOvVOrdered = InsVertivoreDataPairwise$`Omnivore - Vertivore`[order(InsVertivoreDataPairwise$`Omnivore - Vertivore`$p.adj),]


pvalDF = data.frame(InsVertivoreData$p.adj, InsVertvioreHvI$p.adj, InsVertivoreHvO$p.adj, InsVertivoreIvO$p.adj, InsVertivoreHvV$p.adj, InsVertivoreIvV$p.adj, InsVertivoreOvV$p.adj)
pvalDF[pvalDF == 1] = NA
names(pvalDF) = c("Overall", "HvI", "HvO", "IvO", "HvV", "IvV", "OvV")
rownames(pvalDF) = rownames(InsVertivoreData)
pvalDF

pvalOrder = pvalDF[order(pvalDF$Overall),]
pvalOrderTrim = pvalOrder[1:100,]

matrixTest = as.matrix(pvalOrderTrim)

heatmap(matrixTest)


# --------
#Figure 3 look ups
InsectivoryBinaryCorrelations = readRDS("Output/CategoricalBinaryInsectivoreTree/CategoricalBinaryInsectivoreTreeCombinedCategoricalCorrelationFile.rds")[[2]][[1]]
HerbivoryBinaryCorrelations = readRDS("Output/CategoricalBinaryHerbivoreTree/CategoricalBinaryHerbivoreTreeCombinedCategoricalCorrelationFile.rds")[[2]][[1]]
categoricalRER = readRDS("Output/CategoricalInsvertivoreTree/CategoricalInsVertivoreTreeRERFile.rds")
categoricalPaths = readRDS("Output/CategoricalInsvertivoreTree/CategoricalInsVertivoreTreeCategoricalPathsFile.rds")
categoricalPhenv = readRDS("Output/CategoricalInsvertivoreTree/CategoricalInsVertivoreTreeCategoricalPhenotypeVector.rds")
categoricalBigPaths = readRDS("Output/CategoricalMobivoreTree/CategoricalMobivoreTreeCategoricalPathsFile.rds")
InsectBinaryRER = readRDS("Output/CategoricalBinaryInsectivoreTree/CategoricalBinaryInsectivoreTreeRERFile.rds")
InsectBinaryPath = readRDS("Output/CategoricalBinaryInsectivoreTree/CategoricalBinaryInsectivoreTreeCategoricalPathsFile.rds")
HerbivoreBinaryRER = readRDS("Output/CategoricalBinaryHerbivoreTree/CategoricalBinaryHerbivoreTreeRERFile.rds")
HerbivoreBinaryPaths = readRDS("Output/CategoricalBinaryHerbivoreTree/CategoricalBinaryHerbivoreTreeCategoricalPathsFile.rds")

which(row.names(pvalDF) == "EHHADH")


plotRers2 = function (rermat = NULL, index = NULL, phenv = NULL, rers = NULL, 
                      method = "k", xlims = NULL, plot = 1, xextend = 0.2, sortrers = F, categoricalForce = F) 
{
  if (!is.null(phenv) && length(unique(phenv[!is.na(phenv)])) > 
      2 | categoricalForce) {
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
    g <- ggplot(df, aes(x = rer, y = factor(species, levels = unique(ifelse(rep(sortrers, 
                                                                                nrow(df)), species[order(rer)], sort(unique(species))))), 
                        col = mole, label = species)) + scale_size_manual(values = c(1, 
                                                                                     1, 1, 1)) + geom_point(aes(size = mole)) + scale_color_manual(values = c("deepskyblue3", 
                                                                                                                                                              "brown1")) + scale_x_continuous(limits = ll) + geom_text(hjust = 1, 
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
  if (plot) {
    print(g)
  }
  else {
    g
  }
}



which(row.names(InsectivoryBinaryCorrelations) == "EHHADH")
pvalDF[868,]
InsVertvioreHvI[868,]
InsectivoryBinaryCorrelations[868,]
HerbivoryBinaryCorrelations[868,]

plotRers(CategoricalRER, "EHHADH", phenv = categoricalBigPaths)
plotRers(CategoricalRER, "EHHADH", phenv = categoricalPaths)
plotRers2(InsectBinaryRER, "EHHADH", phenv = InsectBinaryPath, categoricalForce = T)
plotRers2(HerbivoreBinaryRER, "EHHADH", phenv = HerbivoreBinaryPaths, categoricalForce = T)


which(row.names(InsectivoryBinaryCorrelations) == "ACAA2")
InsVertvioreHvI[3434,]
InsectivoryBinaryCorrelations[3434,]
HerbivoryBinaryCorrelations[3434,]

plotRers(CategoricalRER, "ACAA2", phenv = categoricalBigPaths)
plotRers(CategoricalRER, "ACAA2", phenv = categoricalPaths)
plotRers2(InsectBinaryRER, "ACAA2", phenv = InsectBinaryPath, categoricalForce = T)
plotRers2(HerbivoreBinaryRER, "ACAA2", phenv = HerbivoreBinaryPaths, categoricalForce = T)


which(row.names(InsectivoryBinaryCorrelations) == "ACADM")
InsVertvioreHvI[9077,]
InsectivoryBinaryCorrelations[9077,]
HerbivoryBinaryCorrelations[9077,]


# --- Chemical synaptyic ----- 

which(row.names(InsectivoryBinaryCorrelations) == "GABRR1")
InsVertvioreHvI[13000,]
InsectivoryBinaryCorrelations[13000,]
HerbivoryBinaryCorrelations[13000,]

which(row.names(InsectivoryBinaryCorrelations) == "STXBP1")
InsVertvioreHvI[9498,]
InsectivoryBinaryCorrelations[9498,]
HerbivoryBinaryCorrelations[9498,]

plotRers(categoricalRER, "STXBP1", phenv = categoricalBigPaths)
plotRers(categoricalRER, "STXBP1", phenv = categoricalPaths)
plotRers2(InsectBinaryRER, "STXBP1", phenv = InsectBinaryPath, categoricalForce = T)
plotRers2(HerbivoreBinaryRER, "STXBP1", phenv = HerbivoreBinaryPaths, categoricalForce = T)

which(row.names(InsectivoryBinaryCorrelations) == "GABRB3")
InsVertvioreHvI[14052,]
InsectivoryBinaryCorrelations[14052,]
HerbivoryBinaryCorrelations[14052,]

plotRers(CategoricalRER, "GABRB3", phenv = categoricalBigPaths)
plotRers(CategoricalRER, "GABRB3", phenv = categoricalPaths)
plotRers2(InsectBinaryRER, "GABRB3", phenv = InsectBinaryPath, categoricalForce = T)
plotRers2(HerbivoreBinaryRER, "GABRB3", phenv = HerbivoreBinaryPaths, categoricalForce = T)



which(row.names(InsectivoryBinaryCorrelations) == "ALDH4A1")
InsVertvioreHvI[9821,]
InsectivoryBinaryCorrelations[9821,]
HerbivoryBinaryCorrelations[9821,]

which(row.names(InsectivoryBinaryCorrelations) == "TAT")
InsVertvioreHvI[7620,]
InsectivoryBinaryCorrelations[7620,]
HerbivoryBinaryCorrelations[7620,]

which(row.names(InsectivoryBinaryCorrelations) == "HIBCH")
InsVertvioreHvI[8011,]
InsectivoryBinaryCorrelations[8011,]
HerbivoryBinaryCorrelations[8011,]

# ------------------
#Proportional venn diagram. 
require(venneuler)
library(venneuler)
library(eulerr)

v <- venneuler(c(
  Insectivore=385, Vertivore=191, Carnivore = 707, 
  "Insectivore&Vertivore"=6, "Insectivore&Carnivore"= 168, "Vertivore&Carnivore"= 49,
  "Vertivore&Carnivore&Insectivore"= 42))
plot(v, , col = c("slateblue", "firebrick", "palevioletred"))

v <- venneuler(c(
  Insectivore=104, Vertivore=35, Carnivore = 1092, 
  "Insectivore&Vertivore"=0, "Insectivore&Carnivore"= 78, "Vertivore&Carnivore"= 25,
  "Vertivore&Carnivore&Insectivore"= 5))
plot(v, , col = c("slateblue", "firebrick", "palevioletred"))


v <- euler(c(
  Insectivore=385, Vertivore=191, Carnivore = 707, 
  "Insectivore&Vertivore"=6, "Insectivore&Carnivore"= 168, "Vertivore&Carnivore"= 49,
  "Vertivore&Carnivore&Insectivore"= 42))
plot(v, fills = c("slateblue", "firebrick", "palevioletred"))

v <- euler(c(
  Insectivore=104, Vertivore=35, Carnivore = 1092, 
  "Insectivore&Vertivore"=0, "Insectivore&Carnivore"= 78, "Vertivore&Carnivore"= 25,
  "Vertivore&Carnivore&Insectivore"= 5))
plot(v, fills = c("slateblue", "firebrick", "palevioletred"))

?venneuler
?plot.euler

# ------ Checking RERplots of specific genes ----- 
categoricalCommonRER = categoricalRER
colnames(categoricalCommonRER)
colnames(categoricalCommonRER) = ZonomNameConvertVectorCommon(colnames(categoricalCommonRER), tipColumn = "ZoonomiaTip")
palette(c( "darkgreen", "darkblue", "black", "red"))
categoricalPaths


plotRers2(categoricalCommonRER, "POLE", phenv = categoricalBigPaths, categoricalForce = T, sortrers = T)
plotRers2(categoricalCommonRER, "POLE", phenv = categoricalPaths, categoricalForce = T, sortrers = T)
plotRers2(categoricalCommonRER, "POLE", phenv = categoricalBigPaths, categoricalForce = T, sortrers = F)

plotRers2(categoricalCommonRER, "MAPK8IP2", phenv = categoricalBigPaths, categoricalForce = T, sortrers = T)
plotRers2(categoricalCommonRER, "MAPK8IP2", phenv = categoricalPaths, categoricalForce = T, sortrers = T)
plotRers2(categoricalCommonRER, "MAPK8IP2", phenv = categoricalBigPaths, categoricalForce = T, sortrers = F)


# --- Violin plots --- 
source("Src/Reu/rerViolinPlot.R")
phenotypeSet = c("Herbivore", "Insectivore", "Omnivore", "Vertivore")
colorScale = c("darkgreen", "darkblue", "black", "red")
rerViolinPlotCategorical(mainTrees, categoricalRER, categoricalBigPaths, phenotypeSet, "POLE", colorScale)

quickViolin = function(geneOfInterest, legend = T){rerViolinPlot(mainTrees, categoricalRER, categoricalBigPaths, phenotypeSet, geneOfInterest, colorScale, NULL, legend)}
quickPrunedViolin = function(geneOfInterest, legend = T){rerViolinPlot(mainTrees, categoricalRER, categoricalPaths, phenotypeSet, geneOfInterest, colorScale, NULL, legend)}
compareViolin = function(geneOfInterest, sideBySide = F){
  if(!sideBySide){
    print(quickViolin(geneOfInterest))
    print(quickPrunedViolin(geneOfInterest))
  }else{
    p1 = quickViolin(geneOfInterest, F)
    p2 = quickPrunedViolin(geneOfInterest, F)
    grid.arrange(p1, p2, ncol=2)
  }
}

compareViolin("POLE")
compareViolin("EHHADH")
compareViolin("ACAA2")
compareViolin("MAPK8IP2")

compareViolin("EHHADH", T)
compareViolin("ACAA2", T)

geneOfInterest = "POLE"

