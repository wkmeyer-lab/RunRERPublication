library(RERconverge)
source("Src/Loc/dev/MergeDataCreation.R")




# ------------------------------- Making your input Data ----------------------------------------

# The first thing we need is a set of genetrees with a mastertree based on the average of the gene trees. This set of trees is called a Trees Object, or a MainTrees. 
# In this case, we will be using the demo shiped with RERConerge. 

rerPath = find.package("RERconverge")
demoTreeFile = paste(rerPath, "/extdata/", "SubsetMammalGeneTrees.txt", sep="")

demoTrees = readTrees(demoTreeFile)

demoTrees$masterTree$tip.label


# The second thing we needs is a phenotype, and a way to connect the phenotype data with the tips on the tree. This is called a MergeData. 
# We will be creating one now. 
# MergeData have a few critical components:
  #The species scientific name
  # The species name as it appears on the tips in the tree
  # Phenotype data for species. 
#In this walkthrough, you have been provided diet phenotype data for the species in the tree. 

#MergeData is made using the "CombineDatasets" script. 
#CombineDatasets requires two things: a column with scientific names of the species, and a column of data you want to add. 
#Your final mergeData requires at least the following: 
  #A column of scientific names
  #A column of common names
  #A column of names as they appear in the Maintrees
  #A column containing your phenotype data. 


#The first table to add into CombineDatasets is the scientific names of species. 
ScientificName = c("Ornithorhynchus anatinus", "Monodelphis domestica", "Sarcophilus harrisii", 
                   "Macropus eugenii", "Dasypus novemcinctus", "Orycteropus afer afer", 
                   "Echinops telfairi", "Chrysochloris asiatica", "Elephantulus edwardii", 
                   "Trichechus manatus latirostris", "Loxodonta africana", "Condylura cristata", 
                   "Erinaceus europaeus", "Sorex araneus", "Felis catus", "Canis lupus familiaris", 
                   "Ailuropoda melanoleuca", "Mustela putorius", "Halichoerus grypus", 
                   "Odobenus rosmarus", "Pteropus vampyrus", "Pteropus giganteus", 
                   "Eptesicus fuscus", "Myotis davidii", "Myotis lucifugus", "Diceros bicornis", 
                   "Equus caballus", "Vicugna pacos", "Camelus bactrianus", "Sus scrofa", 
                   "Delphinus delphis", "Orcinus orca", "Bos taurus", "Pantholops hodgsonii", 
                   "Ovis aries", "Capra hircus", "Tupaia chinensis", "Ochotona princeps", 
                   "Oryctolagus cuniculus", "Allactaga bullata", "Rattus norvegicus", 
                   "Mus musculus", "Microtus arvalis", "Mesocricetus auratus", "Cricetulus griseus", 
                   "Spermophilus tridecemlineatus", "Heterocephalus glaber", "Cavia porcellus", 
                   "Chinchilla lanigera", "Octodon degus", "Otolemur garnettii", 
                   "Saimiri boliviensis", "Callithrix jacchus", "Chlorocebus sabaeus", 
                   "Papio anubis", "Macaca mulatta", "Macaca fascicularis", "Nomascus leucogenys", 
                   "Pongo abelii", "Gorilla gorilla", "Pan troglodytes", "Homo sapiens")
SpeciesScientificNamesDataframe = data.frame(ScientificName)


#We then use the createSortedTable function to create a table with a column with a consistently formatted scientific name that can be used to connect across datasets (explain space, sub spcies, ect)

starterData = createSortedTable(SpeciesScientificNamesDataframe, "ScientificName")

#Next, we add the common names 
CommonName = c("Platypus", "Opossum", "Tasmanian devil", "Wallaby", "Armadillo", 
               "Aardvark", "Tenrec", "Cape golden mole", "Cape elephant shrew", 
               "Manatee", "Elephant", "Star-nosed mole", "Hedgehog", "Shrew", 
               "Cat", "Dog", "Panda", "Ferret", "Grey seal", "Pacific walrus", 
               "Megabat", "Indian flying fox", "Big brown bat", "David's myotis bat", 
               "Microbat", "Black Rhinoceros", "Horse", "Alpaca", "Bactrian Camel", 
               "Pig", "Dolphin", "Killer whale", "Cow", "Tibetan antelope", 
               "Domestic sheep", "Domestic goat", "Chinese tree shrew", "Pika", 
               "Rabbit", "Gobi jerboa", "Rat", "Common mouse", "Common vole", 
               "Golden hamster", "Chinese hamster", "Squirrel", "Naked mole-rat", 
               "Domestic guinea pig", "Chinchilla", "Brush-tailed rat", "Bushbaby", 
               "Squirrel monkey", "Marmoset", "Green monkey", "Olive Baboon", 
               "Rhesus Macaca", "Crab-eating macaque", "Gibbon", "Orangutan", 
               "Gorilla", "Chimp", "Human")

SpeciesCommonNamesDataframe = data.frame(CommonName, ScientificName) #Note that all added data must have scientific name included in the added data


#### DEV NOTE: Figure out why a global combinedData has to exist for this to work. #######
#We then combine these into a single dataset using the CombineDatsets Function. By listing CommonName as a Name column, it is put before the divider between name columns and data columns. 
combinedData = starterData
combinedData = CombineDatasets(combinedDataInput = starterData, newDatasetInput = SpeciesCommonNamesDataframe, newDataScientificNameColumn = "ScientificName", addNewSpeciesValue = T, nameColumns = "CommonName")

#We then need to add phenotype data for our species to the mergedData 

DemoDietPhenotype = c("Carnivore", "Omnivore", "Omnivore", "Herbivore", "Carnivore", 
                      "Carnivore", NA, "Carnivore", "Carnivore", "Herbivore", "Herbivore", 
                      "Carnivore", "Omnivore", "Omnivore", "Carnivore", "Carnivore", 
                      "Omnivore", "Carnivore", NA, "Carnivore", "Herbivore", NA, "Carnivore", 
                      "Carnivore", "Carnivore", "Herbivore", "Herbivore", "Herbivore", 
                      "Omnivore", "Omnivore", NA, "Carnivore", "Herbivore", "Herbivore", 
                      "Herbivore", "Herbivore", NA, "Herbivore", "Herbivore", "Omnivore", 
                      "Omnivore", "Omnivore", NA, "Omnivore", "Herbivore", "Omnivore", 
                      "Herbivore", "Herbivore", "Herbivore", "Herbivore", "Herbivore", 
                      "Omnivore", "Omnivore", "Omnivore", "Omnivore", "Omnivore", "Omnivore", 
                      "Herbivore", "Omnivore", "Herbivore", "Omnivore", "Omnivore")

DietDataframe = data.frame(DemoDietPhenotype, ScientificName)
newDataset = DietDataframe
combinedData = CombineDatasets(combinedDataInput = combinedData, newDatasetInput = DietDataframe, newDataScientificNameColumn = "ScientificName", addNewSpeciesValue = T)

#Finally, we need to add the column of the tip names as they appear in MainTrees. 

demoTreeTipName = demoTrees$masterTree$tip.label
tipNameDataframe = data.frame(demoTreeTipName, ScientificName)

combinedData = CombineDatasets(combinedDataInput = combinedData, newDatasetInput = tipNameDataframe, newDataScientificNameColumn = "ScientificName", addNewSpeciesValue = T, nameColumns = "demoTreeTipName")


MergedData = combinedData 

#Here, you would save your mergedData as a result. 
write.csv(MergedData, "Results/DemoMergedData.csv")



#With your dataset and merged data complete, we can now begin using the functions of runRER. 



# ------------------ Making your Phenotype Tree ---------------------------------------

#The first script to use is MakeCategoricalPhenotypeTree. If you are not using a categorical phenotype, you would use the appropriate phenotype tree creator. 


#These scripts are operated by providing them a list of arguments at the start of the script. This is for compatibility with running on a cluster, which will be covered later. 
#When running a script, the first step is to read the list of args, and what they do. 

#There are a few args which are common to many scripts: 
# r   This inidcates the prefix being used. Each prefix has its own output folder, and allows you to run the script on mulitple analyses without overwriting the other analyses. 
  #In this case, out prefix will be "demo"
# v   This prefix will force the script to de-generate all of its input files, instead of using ones generated in the past. 
  #Useful if you changed something, and need to make updates. If you don't, leave v to False (F), because regenerating files wastes a lot of time. 
#m    This is the filename of your maintrees file. 
  #In this case, this will be the what is stored in the <DemoTreeFile> variaible.

#There are many arguments in Make Categorical Phenotype Tree. other than the three above, here are the ones we will use: 
#d    This is the location of your MergedData Spreadsheet. 
#a    This is the name of the column with the phenotype data. In this case, "DemoDietPhenotype". 
#c    This is a list of the categories you want to include in the analysis. Any categories not listed here will be ignored. This allows the script to ignore species either without phenoytpe data or with phenotypes not relevant to the analysis.
  #In this case, "Carnivore", "Herbivore", "Omnivore". 
#n    This is the column in your spreadsheet with the tip names in MainTrees. 
  #In this case, "demoTreeTipName". 

#The other arguments are involved in remaning phenotype values (u and o), using only a subset of the species in the mainTrees with the relevant phenotypes (s), details of ancestral state reconsturction (t and g) and tree pruning (z, x, y, and p).
  #We can ignore all of these for now -- if you do not include the in your arg string, they will use their default values, which are fine.  


##NOTE: WHEN RUNNING THE SCRIPTs, MAKE SURE TO SKIP THE LINE clusterRun = T AS THIS SETS THE CODE TO WORK ON THE CLUSTER INSTEAD OF A LOCAL MACHINE 

## TEST YOURSELF: What should the args you use for this script be? ##
#####DEV Note: In full version, better explain why args are a string (need to use c), why use ' quotes, why args with mulitple values need another c and use " quotes. ########
##Spoiler this##
args = c('r=Demo', 'v=F', 'm=demoTreeFile', 'd=Results/DemoMergedData.csv', 'a=DemoDietPhenotype', 'c=c("Carnivore", "Herbivore", "Omnivore")', 'n=demoTreeTipName')
##end spoiler##

#If you run that script, you should find a copy of your tree in "Ouput/Demo/DemoCategoricalTree.pdf
#If you look at it, and the branches are colored, good job, it worked! 
#However, the colors might seem a little mis-aligned. Run the following command, and then run the code that generates the pdf (and ONLY that part of the code), and see if which phenotype has which color has changed. 
palette(c("red", "darkgreen", "black"))


# ----------------------------------- Running RERConverge and getting your RER Results --------------------------------------

# Now that you have a phenotype tree (phenoytpe data) and a MainTrees (gene data), we can compare the two, to see if they are linked. 

#The first step is to compare the number of differences in each gene across species, to determine the relative evolutionary rate of each gene. 
  #This is done by the getAllResiduals function (contained in runRER), and is stored in the PrefixRERFile.rds file in the output folder. 
#Once we know which genes are evolving faster and slower, we can then see if a gene's relative evolutionary rate is correlated with our phenotype. 
  #This is done using the various Correlation functions within runRER. These are then output to PrefixCorrelationFile.csv (and .rds) in the output folder. 

#To use RunRERandCorrelation, we do something very similar to Make CategoricalPhenotypeTree; it is also operated using arugments. 

#The arguments we needs are: 
#The shared arguments from before (r, v, and m)
#s    This sets what type of phenotype you are using.
  #In this case, "g" or "categorical" (they do the same thing)

#The other arguments are for overriding the normal functions by manually providing files (p and f), used for continuous phenotypes (c), or change the minimum species required to run correlation on a gene (l)

## TEST YOURSELF: What should the args be for this script? ##
##Spoiler this##
args = c('r=Demo', 'v=F', 'm=demoTreeFile', 's=g')
##End spoiler## 

#Use those arguments, and run the script. This may take a while. 
#Once the script finishes, take a look at the output correlation files. 
####DEV NOTE: In full version explain p, p.adj, MTH correction, and what RHO means ####
###DEV NOTE: May want to switch out phenotype to one that produces significant results from this dataset. ####

#Now that you have run this script, you can see folders in your output folder. These are for each of the pairwise comparisons (For example Herbivore-Carnivore, Herbivore-Omnivore, etc.)
#Now that you have your results, they can still be somewhat hard to interpret. There is a script that can help with that to a degree. 

#This script is VisualizeResultsNew, and is operated like the others. There are two new arguments that new need to change for this script: 
#s  This tells the script which subdirectorys to make a visualization of. We want to do it for all of them. 
  #In this case, it is "Carnivore-Herbivore", "Carnivore-Omnivore", "Herbivore-Omnivore", amd "Overall
#g  This sets is gene enrichment should be included in the visualization. For now, we need to set this to False (F)\
#p  This sets if the 

#The others are involved in visualizing enrichment (l, o, and u) which we will be using later. The others are used for permulations (p and f), which we will be using much later. 

##Test yourself: What should the args for this script be? ##
##Spoiler this##
args = c('r=Demo', 'v=F', 'g=F', 's=c("Carnivore-Herbivore", "Carnivore-Omnivore", "Herbivore-Omnivore", "Overall")')

#ONce that complete, look at the new visualize output files in the OUptu subdirectories. 


# -------------------------- Running Enrichment Analysis -----------------------------------

#Now that we have connected genes with out pehnotype, we can try to see if any groups of genes (pathways, similar functions, etc.) are linked with it. 
#This is done using Gene Enrichment Categories, which are datasets of groups of genes invilved in a shared pathways or process. 
#Various Gene Enrichment sets can be found online, in this case we will use Kegg-Reactome and DisGeNET as examples. 

####Dev Note: Figure out the obtaining of these files when not shipped with repo #####

#To Do this, we will use the EnrichmentAnalysisNew script. It is operated with args, just like the others. 

#For this, there are a few arguments that we have seen before: 
#r, and v 
#s  The subdirectory argument 

#There is a new argument for this script as well, though it shares a letter with other scripts. 
#m    In this case, m is for the gmt files, instead of the MainTrees file. provide it a vector of the gmt files you want to you. 
  #In this case, "Data/KeggReactome.gmt" and "Data/DisGeNET.gmt". 

##Test yourself: What should the args for this script be? ##
##Spoiler this##
args = c('r=Demo', 'v=F', 'm=c("Data/KeggReactome.gmt", "Data/DisGeNET.gmt")', 's=c("Carnivore-Herbivore", "Carnivore-Omnivore", "Herbivore-Omnivore", "Overall")')

#ONce you have the enrichments completed, we can go back and use the visualize script again, but this time with the enrichment turned on. 
##Test yourself: What should the args for the new VisualizeResultsNew be? ##
##Spoiler this##
args = c('r=Demo', 'v=F', 'g=T', 's=c("Carnivore-Herbivore", "Carnivore-Omnivore", "Herbivore-Omnivore", "Overall")')

####DEV NOTE: FIx the need to manually change the number of enrichment plots in visualizeResults####
#####DEV NOTE: Fix the visualizer breaking paritailly when only using 2 GO categories#####

#With that, you have completed a basic RERConverge analysis of this demo phenotype! 
#Later walkthroughs will cover using permulations and use on the cluster. 