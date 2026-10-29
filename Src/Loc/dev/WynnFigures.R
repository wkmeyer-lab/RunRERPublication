library(RERconverge)
library(tools)
library(scales)
library(data.table)

source("Src/Reu/ZoonomTreeNameToCommon.R")




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

#---
#Make the RER Plot 
#---
png("Output/Misc/WynnCategoricalPlot.png", width = 2000, height = 2000)
plotRers(wynnCategoricalRERCommon, difGene, wynnCategoricalPath, sortrers = T)
dev.off()


#---
# Make the Scatter plot 
#---
phenotypeTree = readRDS("Output/HarshalCategoricalRER/HarshalCategoricalRERCategoricalTree.rds")
speciesFilter = readRDS("Output/HarshalContinousRERMod/HarshalContinousRERModSpeciesFilter.rds")

RelativeEvolutionaryRate = wynnContinousRERCommon[which(rownames(wynnContinousRERCommon) == difGene),]
ChangeinTMM = wynnContinousPath
matchedPathsObject = tree2Paths(phenotypeTree, mainTrees, useSpecies = speciesFilter, categorical = TRUE) #do not binarize; the categorical data is already contained in the phenotype tree.

scatterPlotData = data.frame(RelativeEvolutionaryRate, ChangeinTMM, RERNames = names(RelativeEvolutionaryRate), categoricalPath = matchedPathsObject)


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


