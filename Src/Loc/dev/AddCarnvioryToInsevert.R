pairwiseCorrelationFileName = paste(outputFolderName, filePrefix, "PairwiseCorrelationFile.rds", sep= "") #make a name for the pairwise comparisons based on prefix

correlationResults = readRDS(pairwiseCorrelationFileName)


CarnivoreGeneData = readRDS("Output/CategoricalPrunedCarnivoreTree/Carnivore-Herbivore/CategoricalPrunedCarnivoreTreeCarnivore-HerbivoreCorrelationFile.rds")

correlationResults$'Carnivore - Herbivore' = CarnivoreGeneData
saveRDS(correlationResults, pairwiseCorrelationFileName)

# also coppy across the output folder 