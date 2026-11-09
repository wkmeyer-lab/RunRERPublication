#Check association between RER for a gene of interest and a quantitative trait
#of interest, using overlapping branches from the diet analysis.
#Additionally/alternatively, check distribution of inferred quantitative phenotype or
#inferred change in quantitative phenotype across branches with different
#categorical diets in the analyses.

library(RERconverge) #check version

#Need: Diet phenotype tree, RER matrix, quantitative phenotype values, gene name,
# trees object/master tree, (maybe) tip label conversion table
#Read in quantitative phenotype values
qphenofile <- "Data/TEContentOsmanski2023TS8.csv"
qphenobase <- "TEContent" #base filename for saving outputs
myqpheno <- read.csv(qphenofile)
#Find a way to specify the species name column and phenotype column
qspcol <- "Species"
qpcol <- "Young_DNA_proportion"

#Read in diet phenotype tree
dcattreefile <- "Output/ComplexDietCentralAnalysisSimplify2CategoricalScientificTree.rds"
mdcattree <- readRDS(dcattreefile) #tip labels will need conversion

#Find overlap between species
print(paste("Number of shared species between categorical",
            "diet phenotype tree and quantitative phenotype",
            "vector:",
            length(intersect(tolower(myqpheno[,qspcol]),
                       mdcattree$tip.label))))

#Read in trees object (any way to use just master tree?)
atreesfile <- "Data/zoonomiaAllMammalsTrees.rds"
myatrees <- readRDS(atreesfile)

#This has species names as vs_HL.... So will need to convert first.
convfile <- "Data/mergedData_try2.csv"
myconv <- read.csv(convfile)
#Want to convert from 'Scientific_Binomial' to 'ZoonomiaTip'
relcolsconv <- c("Scientific_Binomial", "ZoonomiaTip")
myconv2 <- myconv[,relcolsconv]
myqpheno$lowSpecies <- tolower(myqpheno[,qspcol]) #make lowercase
myqpheno2 <- merge(myqpheno, myconv2, by.x="lowSpecies", 
                   by.y = relcolsconv[1], all.x=T)
qspcol2 <- relcolsconv[2] #new species column name

#Use char2Paths to create a paths vector for the quantitative phenotype
#Does this work if the species are a subset?
#Make the phenotype a named vector
myqpvec <- as.vector(myqpheno2[,qpcol])
names(myqpvec) <- myqpheno2[,qspcol2]

#Try first: just do char2Paths
myqpathsdiff <- char2Paths(tip.vals = myqpvec, treesObj = myatrees,
                           metric = "diff") #this does run....

#The paths object has diffs for all possible branches in trees containing
#the species in the input dataset. We want to subset to the maximum branch set
#with overlap between the quantitative phenotype data and the categorical diet
#data.

#We could prune the categorical tree and run tree2Paths to get the paths for
#just the species in the overlap set. Then we would need to use some of the code
#from getAllCor to find the right indices within the paths to represent the
#branches for the tree containing all overlapping species.

#Or this may not be necessary. We could create an RER matrix with just a subset
#of genes, using just the overlapping species as useSpecies and no filtering.
#Then we determine which of the row values are not NA; these are the elements of
#the paths that we want to use.
#The 'doOnly' option for 'getAllResiduals' allows computing RERs for a single
#gene (supposedly). If this doesn't work, we could use 'maxT' (maximum number
#of trees).
#'doOnly' is numeric, so it does not take a gene name.
#Find a gene with data for all species in the trees file
gallsp <- which(rowSums(myatrees$report)==dim(myatrees$report)[2])
myasgene <- names(gallsp)[1]

#Read in categorical diet tree file with 'ZoonomiaTip' names (exists?)
ztcattreefile <- "Output/ComplexDietCentralAnalysisSimplify2CategoricalTree.rds"
myztcattree <- readRDS(ztcattreefile) #tip labels will need conversion

#Find overlap between species
myospecies <- intersect(names(myqpvec), myztcattree$tip.label)
print(paste("Number of shared species between categorical",
            "diet phenotype tree and quantitative phenotype",
            "vector, using tip labels:",
            length(unique(myospecies))))

#Test getting all residuals for just one gene, just overlapping species
mytestgar <- getAllResiduals(treesObj = myatrees, maxT = 5)
mytestrer <- getAllResiduals(treesObj = myatrees, cutoff = 0,
                             useSpecies = myospecies, doOnly = myasgene)
#This runs for rather a long time.... If it works, save this as an RDS object.
saveRDS(mytestrer, file=paste("Output/RERFor",myasgene,"SpeciesOverlapWith",
                              qphenobase,".rds",sep=""))

#If it works:

#To do: Subset the paths to the branches in the diet phenotype
# tree. Parse this to compare the distribution of (1) diffs, and 
# (2) mean values of the quantitative phenotype across diets.

#The relevant functions for subsetting paths are likely in the
#RERconverge correlation and/or plotting functions.
#I may have used them in the past for violin plots of RERs across
#categories.

#Do some exploratory plotting.
#The file is saved as a list, with names for each table.
#The values are stored in "edgeTable" within each element of the list.
#We will want to do a test of whether the distribution of continuous values 
#differs across branch phenotypes.
#We may also want to visualize this as a violin plot.
#For summary outputs, we will want the sample size per phenotype, and the 
#test statistics and P-values for diffs and means.
mydres <- readRDS("DnaReapirREsultsForWynn.rds")

ageDiffTab <- mydres[[1]]$edgeTable
ageMeanTab <- mydres[[2]]$edgeTable
TeDiffTab <- mydres[[3]]$edgeTable
TeMeanTab <- mydres[[4]]$edgeTable
for (i in c(1:6)) {
  mytab <- mydres[[i]]$edgeTable
  #Also include log transformation
  myconst <- min(abs(mytab$`Continuous Value`)[
    which(mytab$`Continuous Value` != 0)])/2
  mytab$`Log Continuous Value` <- sign(mytab$`Continuous Value`)*
      log(abs(mytab$`Continuous Value` + myconst))
  ggplot(data=mytab, aes(x=`Continuous Value`)) +
    geom_histogram() +
    facet_wrap(~Phenotype) +
    #geom_text(x=0,y=0,label=names(mydres)[i])
    labs(x=names(mydres)[i])
  ggsave(paste("Plots/HistogramsFor",names(mydres)[i],
               ".jpeg",sep=""))
  ggplot(data=mytab, aes(x=`Log Continuous Value`)) +
    geom_histogram() +
    facet_wrap(~Phenotype) +
    #geom_text(x=0,y=0,label=paste("log",names(mydres)[i]))
    labs(x=paste("log",names(mydres)[i]))
  ggsave(paste("Plots/HistogramsForLog",names(mydres)[i],
               ".jpeg",sep=""))
  print(names(mydres)[i])
  #kruskal.test(`Continuous Value` ~ Phenotype, data=mytab)
  dunn.test(mytab$`Continuous Value`, g = mytab$Phenotype,
            method="bh")
  print(paste("log",names(mydres)[i]))
  #kruskal.test(`Log Continuous Value` ~ Phenotype, data=mytab)
  dunn.test(mytab$`Log Continuous Value`, g=mytab$Phenotype,
            method="bh")
}



