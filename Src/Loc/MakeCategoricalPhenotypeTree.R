# -- Libraries 
clusterRun = F
clusterRun = T
if(clusterRun){.libPaths("/share/ceph/wym219group/shared/libraries/R4")} #add path to custom libraries to searched locations
library(RERconverge)
library(tools)
source("Src/Reu/cmdArgImport.R")
source("Src/Reu/ZoonomTreeNameToCommon.R")
# -- Usage:
# This script creates a categorical tree of a phenotype which has been annotated in the Manual Annotations spreadsheet of the meyer lab. 
# In theory, this script could be used on any spreadsheet, so long as the column containing the tip.labels is specified using the n argument, and the column with common names is named "CommonName". 

# -- Command arguments list
# r = filePrefix                                         This is a prefix used to organize and separate files by analysis run. Always required. 
# v = <T or F>                                           This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# m = mainTreeFilename.txt or .rds                       This sets the location of the maintrees file
# d = spreadSheetFilename.csv                            This sets the spreadsheet to read the data from 
# a = "annotCollumn"                                     This is the column in the manual annotations spreadsheet to use
# c = <c("nameOfCategory1,"nameOfCategory2")>            This is the list of category names 
# u = list(c("replace1", "with1"),c("replace2, with2"))
# o = list(c("phenotype1", "intophen1"), c("2", "i2"))   This causes combination phenotypes to be merged into the second phenotype, but does not replace standalone phenotypes
# s = "screenCollumn"                                    This is a collumn which must have a value of 1 for the species to be included. 
# t = <ER or SYM or ARD>                                 This sets the model type used to estimate ancestral branches 
# g = "ancestralTrait"                                   This can be used to set all non-terminal branches to this category. Use be one of the categories in the list. 
# n = "nameColumn"                                       This sets the column with the tip names as they appear in the maintrees file. 
# z = <minimum branch length>                            This sets the minimum branch length for terminal branches in the master tree. Branches shorter than this will be removed. 
# x = "pruningPrefrenceColumn"                           This sets a column, where if the value is 1, the tip will be preferentially kept. If the value is TRUE, the tip will never be pruned.
# y = "c('unprunedtip1', 'unprunedtip2')"                This allows you to add a list of specific tips to not be dropped during pruning. Must use the tip name, not common name. 
# p = "c('prunedtip1', 'prunedtip2')"                    This allows you to manually specify additional branches to be pruned

#----------------
{
args = c('r=CategoricalInsectRoot4Phen', 'a=Meyer.Lab.Classification', 'c=c("Carnivore", "Omnivore", "Herbivore", "Insectivore")', 'u=list(c("Generalist","_Omnivore"),c("Omnivore","_Omnivore"), c("Piscivore", "Carnivore"))',   'm=data/RemadeTreesAllZoonomiaSpecies.rds', 'v=T', 't=ER', "n=Insectivore")
args = c('r=BinaryCVHApplesToApples', 'a=Meyer.Lab.Classification', 'c=c("Carnivore", "Herbivore")', 'u=list(c("Carnivore","_Carnivore"), c("Piscivore", "_Carnivore"))',   'm=data/RemadeTreesAllZoonomiaSpecies.rds', 'v=T', 't=ER')
args = c('r=CategoricalER5Phen', 'a=Meyer.Lab.Classification', 'c=c("Carnivore", "Omnivore", "Herbivore", "Insectivore", "Piscivore")', 'u=list(c("Generalist","_Omnivore"),c("Omnivore","_Omnivore"))',   'm=data/RemadeTreesAllZoonomiaSpecies.rds', 'v=T', 't=ER')
args = c('r=CategoricalARD5Phen', 'a=Meyer.Lab.Classification', 'c=c("Carnivore", "Omnivore", "Herbivore", "Insectivore", "Piscivore")', 'u=list(c("Generalist","_Omnivore"),c("Omnivore","_Omnivore"))',   'm=data/RemadeTreesAllZoonomiaSpecies.rds', 'v=T', 't=ARD')
args = c('r=CategoricalARD3Phen', 'a=Meyer.Lab.Classification', 'c=c("Carnivore", "Omnivore", "Herbivore")', 'u=list(c("Generalist","_Omnivore"),c("Omnivore","_Omnivore")), c("Piscivore", "Carnivore")',   'm=data/RemadeTreesAllZoonomiaSpecies.rds', 'v=T', 't=ARD')
args = c('r=RubyRegenARD',   'm=data/mam120aa_trees.rds', 'v=F', 't=ARD')
args = c('r=RubyRegenER',   'm=data/mam120aa_trees.rds', 'v=F', 't=ER', 'a=Meyer.Lab.Classification')
args = c('r=Categorical3PhenARDTest', 'a=Meyer.Lab.Classification', 'c=c("Carnivore", "Herbivore", "Omnivore")', 'u=list(c("Omnivore","_Omnivore"), c("Piscivore", "Carnivore"))',   'm=data/RemadeTreesAllZoonomiaSpecies.rds', 'v=T', 't=rm')
args = c('r=Categorical5PhenDemo', 'a=Meyer.Lab.Classification', 'c=c("Carnivore", "Omnivore", "Herbivore", "Insectivore", "Piscivore")', 'u=list(c("Generalist","_Omnivore"),c("Omnivore","_Omnivore"))',   'm=data/RemadeTreesAllZoonomiaSpecies.rds', 'v=T', 't=ER')



args = c('r=OnetwentyWay6Phen', 'm=data/mam120aa_trees.rds', 'a=Meyer.Lab.Classification', 'c=c("Carnivore", "Omnivore", "Herbivore", "Insectivore", "Piscivore", "Generalist")', 'u=list(c("Generalist","Anthropivore"), c("Generalist", "Omnivore"))', 'o=list(c("Carnivore", "Piscivore"))', 'v=F', 't=SYM')
args = c('r=ThreePhenLikeihoodTest', 'm=data/mam120aa_trees.rds', 'a=Meyer.Lab.Classification', 'c=c("Carnivore", "Omnivore", "Herbivore", "Piscivore", "Generalist", "Insectivore")', 'u=list(c("Anthropivore","_Omninivore"), c("Omnivore", "_Omninivore"), c("Piscivore", "Carnivore"), c("Insectivore", "Carnivore"))', 'o=list(c("Carnivore", "Piscivore"))', 'v=F', 't=ARD')
args = c('r=HMGRelaxTest', 'm=data/mam120aa_trees.rds', 'a=DEBUG_using_preexisting_phenotypeVector', 'v=F', 't=ER')
args = c('r=HMGUnRelaxTest', 'm=data/mam120aa_trees.rds', 'a=DEBUG_using_preexisting_phenotypeVector', 'v=F', 't=ER')
args = c('r=IPCRelaxTest', 'm=data/mam120aa_trees.rds', 'a=DEBUG_using_preexisting_phenotypeVector', 'v=F', 't=ER')
args = c('r=IPCUnRelaxTest', 'm=data/mam120aa_trees.rds', 'a=DEBUG_using_preexisting_phenotypeVector', 'v=F', 't=ER')


args = c('r=NewHiller4Phen', 'm=data/newHillerMainTrees.rds', 'd=Data/HillerZoonomiaPhenotypeTable.csv', 'a=phenotype', 'c=c("Carnivore", "Omnivore", "Herbivore", "Insectivore", "Piscivore", "Generalist")', 'u=list(c("Generalist","Omnivore"), c("Omnivore", "_Omnivore"), c("Piscivore", "Carnivore"))', 'v=F', 't=ER')
args = c('r=NewHiller2Phen', 'm=data/newHillerMainTrees.rds', 'd=Data/HillerZoonomiaPhenotypeTable.csv', 'a=phenotype', 'c=c("Carnivore", "Herbivore", "Piscivore")', 'u=list(c("Piscivore", "Carnivore"))', 'v=F', 't=ER')
args = c('r=PruningTest', 'm=data/RemadeTreesAllZoonomiaSpecies.rds', 'd=Data/mergedData.csv', 'a=Meyer.Lab.Classification', 'c=c("Carnivore", "Omnivore", "Herbivore", "Insectivore", "Piscivore", "Generalist")', 'u=list(c("Generalist", "Omnivore"), )', 'o=list(c("Carnivore", "Piscivore"))','v=T', 't=ER', 'n=Zoonomia', "z=0.01", "x=HillerZoonomiaOverlap")


args = c('r=NewHiller4Phen', 'm=data/newHillerMainTrees.rds', 'd=Data/HillerZoonomiaPhenotypeTable.csv', 'a=phenotypeSimplified', 'c=c("Carnivore", "Omnivore", "Herbivore", "Insectivore", "Piscivore", "Generalist")', 'u=list(c("Generalist","Omnivore"), c("Omnivore", "_Omnivore"), c("Piscivore", "Carnivore"))', 'o=list(c("Carnivore", "Insectivore"))','v=T', 't=ER')
args = c('r=NewHiller2Phen', 'm=data/newHillerMainTrees.rds', 'd=Data/HillerZoonomiaPhenotypeTable.csv', 'a=phenotypeSimplified', 'c=c("Carnivore", "Herbivore", "Piscivore")', 'u=list(c("Herbivore", "_Herbivore"), c("Piscivore", "Carnivore"))','v=T', 't=ER')
args = c('r=NewHillerTestSupraPrimates', 'm=data/newHillerMainTrees.rds', 'd=Data/HillerZoonomiaPhenotypeTable.csv', 'a=Supraprimates', 'c=c("1", "0")', 'u=list(c("Herbivore", "_Herbivore"), c("Piscivore", "Carnivore"))','v=T', 't=ER')

args = c('r=CIvHBinaryHiller', 'm=data/newHillerMainTrees.rds', 'd=Data/mergedData.csv', 'a=Meyer.Diet.Categorical', 'c=c("Carnivore", "Herbivore", "Insectivore", "Piscivore")', 'u=list(c("Herbivore", "_Herbivore"), c("Piscivore", "Carnivore"), c("Insectivore", "Carnivore"))','v=T', 't=ER', 'n=HillerName', 's=NoAoudad')
args = c('r=CIvHBinaryZoonomia', 'm=data/RemadeTreesAllZoonomiaSpecies.rds', 'd=Data/mergedData.csv', 'a=Meyer.Lab.Classification', 'c=c("Carnivore", "Herbivore", "Insectivore", "Piscivore")', 'u=list(c("Herbivore", "_Herbivore"), c("Piscivore", "Carnivore"), c("Insectivore", "Carnivore"))','v=T', 't=ER', 'n=Zoonomia', 's=NoAoudad')


args = c('r=ZoonomiaCategoricalRefrenceTree', 'm=data/RemadeTreesAllZoonomiaSpecies.rds', 'd=Data/mergedData.csv', 'a=Meyer.Lab.Classification.Clean', 'c=c("Carnivore", "Omnivore", "Herbivore", "Insectivore", "Piscivore", "Generalist", "Planktivore")', 'u=list(c("Generalist", "Omnivore"), c("Omnivore-IH", "Omnivore"), c("Omnivore", "_Omnivore"))', 'o=list(c("Piscivore", "Carnivore"), c("Planktivore", "Carnivore"), c("Insectivore", "Carnivore"), c("Piscivore", "Insectivore"))','v=T', 't=ER', 'n=Zoonomia')
args = c('r=HillerCategoricalRefrenceTree', 'm=data/NewHillerMainTrees.rds', 'd=Data/mergedData.csv', 'a=Meyer.Lab.Classification.Clean', 'c=c("Carnivore", "Omnivore", "Herbivore", "Insectivore", "Piscivore", "Generalist", "Planktivore")', 'u=list(c("Generalist", "Omnivore"), c("Omnivore-IH", "Omnivore"), c("Omnivore", "_Omnivore"))', 'o=list(c("Piscivore", "Carnivore"), c("Planktivore", "Carnivore"), c("Insectivore", "Carnivore"), c("Piscivore", "Insectivore"))','v=T', 't=ER', 'n=HillerName', 's=NoAoudad')

args = c('r=TrueCategoricalRefrenceTree', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=Meyer.Lab.Classification.Clean', 'c=c("Carnivore", "Omnivore", "Herbivore", "Insectivore", "Piscivore", "Generalist", "Planktivore")', 'u=list(c("Generalist", "Omnivore"), c("Omnivore-IH", "Omnivore"), c("Omnivore", "_Omnivore"))', 'o=list(c("Piscivore", "Carnivore"), c("Planktivore", "Carnivore"), c("Insectivore", "Carnivore"), c("Piscivore", "Insectivore"))','v=T', 't=ER', 'n=Zoonomia')
args = c('r=makeLalithaTree', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=AlliePheno', 'c=c("1", "2", "3")' ,'v=T', 't=ER', 'n=ZoonomiaTip')

}
{
args = c('r=TrueCategoricalRefrenceTreeComplex', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90', 
          'c=c("C-Invertebrate-eater", "C-Endotherm-Carnivore", "C-Herpetivore", "C-Piscivore", "C-Nonspecific-Vertebrate-eater", "C-Scavenger", 
              "O-For Examination", "O-Scavenger", 
              "H-Frugivore", "H-Nectarivore", "H-Granivore", "H-Nonspecific-Herbivore", 
              "C-Terrestrial-vertebrates-eater", "C-All-vertebrate-eater", "C-All-Animals-Eater", 
              "H-High-sugar-plants-Eater", "H-Low-sugar-plants-Eater", "H-All-plants-Eater", 
              "O-Generalist", 
              "Insectivore", "Herpetivore", "Piscivore", "Vertivore", "InsVertivore", "Omnivore", "Frugivore", "Nectarivore", "Glucivore", "Herbivore", "Generalist"
            )', 
          'u=list(
            c("C-Invertebrate-eater", "Insectivore"), 
            c("C-Herpetivore", "Herpetivore"),
            c("C-Piscivore", "Piscivore"),
            c("C-Endotherm-Carnivore", "Vertivore"), c("C-Scavenger", "Vertivore"), c("C-Nonspecific-Vertebrate-eater", "Vertivore"),
            c("C-Terrestrial-vertebrates-eater", "Vertivore"), c("C-All-vertebrate-eater", "Vertivore"), 
            c("C-All-Animals-Eater", "InsVertivore"), 
            c("O-For Examination", "Omnivore"), c("O-Scavenger", "Omnivore"),
            c("H-Frugivore", "Frugivore"), 
            c("H-Nectarivore", "Nectarivore"), 
            c("H-High-sugar-plants-Eater", "Glucivore"),
            c("H-Granivore", "Herbivore"), c("H-Nonspecific-Herbivore", "Herbivore"), 
            c("H-Low-sugar-plants-Eater", "Herbivore"), c("H-All-plants-Eater", "Herbivore"),
            c("O-Generalist", "Generalist")
          )', 
          'v=T', 't=ER', 'n=ZoonomiaTip')

args = c('r=TrueCategoricalRefrenceTree', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90', 
         'c=c("C-Invertebrate-eater", "C-Endotherm-Carnivore", "C-Herpetivore", "C-Piscivore", "C-Nonspecific-Vertebrate-eater", "C-Scavenger", 
              "O-For Examination", "O-Scavenger", 
              "H-Frugivore", "H-Nectarivore", "H-Granivore", "H-Nonspecific-Herbivore", 
              "C-Terrestrial-vertebrates-eater", "C-All-vertebrate-eater", "C-All-Animals-Eater", 
              "H-High-sugar-plants-Eater", "H-Low-sugar-plants-Eater", "H-All-plants-Eater", 
              "O-Generalist", 
              "Insectivore", "Herpetivore", "Piscivore", "Vertivore", "InsVertivore", "Omnivore", "Frugivore", "Nectarivore", "Glucivore", "Herbivore", "Generalist"
            )', 
         'u=list(
            c("C-Invertebrate-eater", "Insectivore"), 
            c("C-Herpetivore", "Herpetivore"),
            c("C-Piscivore", "Piscivore"),
            c("C-Endotherm-Carnivore", "Vertivore"), c("C-Scavenger", "Vertivore"), c("C-Nonspecific-Vertebrate-eater", "Vertivore"),
            c("C-Terrestrial-vertebrates-eater", "Vertivore"), c("C-All-vertebrate-eater", "Vertivore"), 
            c("C-All-Animals-Eater", "InsVertivore"), 
            c("O-For Examination", "Omnivore"), c("O-Scavenger", "Omnivore"),
            c("H-Frugivore", "Frugivore"), 
            c("H-Nectarivore", "Nectarivore"), 
            c("H-High-sugar-plants-Eater", "Glucivore"),
            c("H-Granivore", "Herbivore"), c("H-Nonspecific-Herbivore", "Herbivore"), 
            c("H-Low-sugar-plants-Eater", "Herbivore"), c("H-All-plants-Eater", "Herbivore"),
            c("O-Generalist", "Omnivore")
          )', 
         'v=T', 't=ER', 'n=ZoonomiaTip')

args = c('r=TrueCategoricalRefrenceTreeCompressed', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90', 
         'c=c("C-Invertebrate-eater", "C-Endotherm-Carnivore", "C-Herpetivore", "C-Piscivore", "C-Nonspecific-Vertebrate-eater", "C-Scavenger", 
              "O-For Examination", "O-Scavenger", 
              "H-Frugivore", "H-Nectarivore", "H-Granivore", "H-Nonspecific-Herbivore", 
              "C-Terrestrial-vertebrates-eater", "C-All-vertebrate-eater", "C-All-Animals-Eater", 
              "H-High-sugar-plants-Eater", "H-Low-sugar-plants-Eater", "H-All-plants-Eater", 
              "O-Generalist", 
              "Insectivore", "Herpetivore", "Piscivore", "Vertivore", "InsVertivore", "Omnivore", "Frugivore", "Nectarivore", "Glucivore", "Herbivore", "Generalist"
            )', 
         'u=list(
            c("C-Invertebrate-eater", "Insectivore"), 
            c("C-Herpetivore", "Vertivore"),
            c("C-Piscivore", "Vertivore"),
            c("C-Endotherm-Carnivore", "Vertivore"), c("C-Scavenger", "Vertivore"), c("C-Nonspecific-Vertebrate-eater", "Vertivore"),
            c("C-Terrestrial-vertebrates-eater", "Vertivore"), c("C-All-vertebrate-eater", "Vertivore"), 
            c("C-All-Animals-Eater", "InsVertivore"), 
            c("O-For Examination", "Omnivore"), c("O-Scavenger", "Omnivore"),
            c("H-Frugivore", "Herbivore"), 
            c("H-Nectarivore", "Herbivore"), 
            c("H-High-sugar-plants-Eater", "Herbivore"),
            c("H-Granivore", "Herbivore"), c("H-Nonspecific-Herbivore", "Herbivore"), 
            c("H-Low-sugar-plants-Eater", "Herbivore"), c("H-All-plants-Eater", "Herbivore"),
            c("O-Generalist", "Omnivore")
          )', 
         'v=T', 't=ER', 'n=ZoonomiaTip')

args = c('r=CategoricalWorkingTree', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90', 
         'c=c("C-Invertebrate-eater", "C-Endotherm-Carnivore", "C-Herpetivore", "C-Piscivore", "C-Nonspecific-Vertebrate-eater", "C-Scavenger", 
              "O-For Examination", "O-Scavenger", 
              "H-Frugivore", "H-Nectarivore", "H-Granivore", "H-Nonspecific-Herbivore", 
              "C-Terrestrial-vertebrates-eater", "C-All-vertebrate-eater", "C-All-Animals-Eater", 
              "H-High-sugar-plants-Eater", "H-Low-sugar-plants-Eater", "H-All-plants-Eater", 
              "O-Generalist", 
              "Insectivore", "Herpetivore", "Piscivore", "Vertivore", "InsVertivore", "Omnivore", "Frugivore", "Nectarivore", "Glucivore", "Herbivore", "Generalist"
            )', 
         'u=list(
            c("C-Invertebrate-eater", "Insectivore"), 
            c("C-Herpetivore", "Vertivore"),
            c("C-Piscivore", "Vertivore"),
            c("C-Endotherm-Carnivore", "Vertivore"), c("C-Scavenger", "Vertivore"), c("C-Nonspecific-Vertebrate-eater", "Vertivore"),
            c("C-Terrestrial-vertebrates-eater", "Vertivore"), c("C-All-vertebrate-eater", "Vertivore"), 
            c("C-All-Animals-Eater", "InsVertivore"), 
            c("O-For Examination", "Omnivore"), c("O-Scavenger", "Omnivore"),
            c("H-Frugivore", "Herbivore"), 
            c("H-Nectarivore", "Herbivore"), 
            c("H-High-sugar-plants-Eater", "Herbivore"),
            c("H-Granivore", "Herbivore"), c("H-Nonspecific-Herbivore", "Herbivore"), 
            c("H-Low-sugar-plants-Eater", "Herbivore"), c("H-All-plants-Eater", "Herbivore"),
            c("O-Generalist", "Omnivore")
          )', 
         'v=T', 't=ER', 'n=ZoonomiaTip')

args = c('r=CategoricalMobivoreTree', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 
         'c=c("C-Invertebrate-eater", "C-Endotherm-Carnivore", "C-Herpetivore", "C-Piscivore", "C-Nonspecific-Vertebrate-eater", "C-Scavenger", 
              "O-For Examination", "O-Scavenger", 
              "H-Frugivore", "H-Nectarivore", "H-Granivore", "H-Nonspecific-Herbivore", 
              "C-Terrestrial-vertebrates-eater", "C-All-vertebrate-eater", "C-All-Animals-Eater", 
              "H-High-sugar-plants-Eater", "H-Low-sugar-plants-Eater", "H-All-plants-Eater", 
              "O-Generalist", 
              "C-InsVertivore-Mixed", "C-InsVertivore-Piscivore", "C-InsVertivore-Insectivore","C-InsVertivore-Carnivore",
              "Insectivore", "Herpetivore", "Piscivore", "Vertivore", "InsVertivore", "Omnivore", "Frugivore", "Nectarivore", "Glucivore", "Herbivore", "Generalist"
            )', 
         'u=list(
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
          )', 
         'v=T', 't=ER', 'n=ZoonomiaTip')
args = c('r=RyanPermulationsDemo', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 
         'c=c("C-Invertebrate-eater", "C-Endotherm-Carnivore", "C-Herpetivore", "C-Piscivore", "C-Nonspecific-Vertebrate-eater", "C-Scavenger", 
              "O-For Examination", "O-Scavenger", 
              "H-Frugivore", "H-Nectarivore", "H-Granivore", "H-Nonspecific-Herbivore", 
              "C-Terrestrial-vertebrates-eater", "C-All-vertebrate-eater", "C-All-Animals-Eater", 
              "H-High-sugar-plants-Eater", "H-Low-sugar-plants-Eater", "H-All-plants-Eater", 
              "O-Generalist", 
              "C-InsVertivore-Mixed", "C-InsVertivore-Piscivore", "C-InsVertivore-Insectivore","C-InsVertivore-Carnivore",
              "Insectivore", "Herpetivore", "Piscivore", "Vertivore", "InsVertivore", "Omnivore", "Frugivore", "Nectarivore", "Glucivore", "Herbivore", "Generalist"
            )', 
         'u=list(
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
          )', 
         'v=T', 't=ER', 'n=ZoonomiaTip')

args = c('r=CategoricalCarnivoreTree', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 
         'c=c("C-Invertebrate-eater", "C-Endotherm-Carnivore", "C-Herpetivore", "C-Piscivore", "C-Nonspecific-Vertebrate-eater", "C-Scavenger", 
              "O-For Examination", "O-Scavenger", 
              "H-Frugivore", "H-Nectarivore", "H-Granivore", "H-Nonspecific-Herbivore", 
              "C-Terrestrial-vertebrates-eater", "C-All-vertebrate-eater", "C-All-Animals-Eater", 
              "H-High-sugar-plants-Eater", "H-Low-sugar-plants-Eater", "H-All-plants-Eater", 
              "O-Generalist", 
              "C-InsVertivore-Mixed", "C-InsVertivore-Piscivore", "C-InsVertivore-Insectivore","C-InsVertivore-Carnivore",
              "Insectivore", "Herpetivore", "Piscivore", "Vertivore", "InsVertivore", "Omnivore", "Frugivore", "Nectarivore", "Glucivore", "Herbivore", "Generalist"
            )', 
         'u=list(
            c("C-Invertebrate-eater", "Carnivore"), c("C-InsVertivore-Insectivore", "Carnivore"),
            c("C-Herpetivore", "Carnivore"),
            c("C-Piscivore", "Carnivore"), c("C-InsVertivore-Piscivore", "Carnivore"),
            c("C-Endotherm-Carnivore", "Carnivore"), c("C-Scavenger", "Carnivore"), c("C-Nonspecific-Vertebrate-eater", "Carnivore"),
            c("C-Terrestrial-vertebrates-eater", "Carnivore"), c("C-All-vertebrate-eater", "Carnivore"), c("C-InsVertivore-Carnivore", "Carnivore"),
            c("C-InsVertivore-Mixed", "Omnivore"), 
            c("O-For Examination", "Omnivore"), c("O-Scavenger", "Omnivore"),
            c("H-Frugivore", "Herbivore"), 
            c("H-Nectarivore", "Herbivore"), 
            c("H-High-sugar-plants-Eater", "Herbivore"),
            c("H-Granivore", "Herbivore"), c("H-Nonspecific-Herbivore", "Herbivore"), 
            c("H-Low-sugar-plants-Eater", "Herbivore"), c("H-All-plants-Eater", "Herbivore"),
            c("O-Generalist", "Omnivore")
          )', 
         'v=T', 't=ER', 'n=ZoonomiaTip')
}
args = c('r=CategoricalInsVertivoreTree', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01',
         'c=c(
              "C-Invertebrate-eater", "C-Endotherm-Carnivore", "C-Herpetivore", "C-Piscivore", "C-Nonspecific-Vertebrate-eater", "C-Scavenger", 
              "O-For Examination", "O-Scavenger", 
              "H-Frugivore", "H-Nectarivore", "H-Granivore", "H-Nonspecific-Herbivore", 
              "C-Terrestrial-vertebrates-eater", "C-All-vertebrate-eater", "C-All-Animals-Eater", 
              "H-High-sugar-plants-Eater", "H-Low-sugar-plants-Eater", "H-All-plants-Eater", 
              "O-Generalist", 
              "C-InsVertivore-Mixed", "C-InsVertivore-Piscivore", "C-InsVertivore-Insectivore","C-InsVertivore-Carnivore",
              "Insectivore", "Herpetivore", "Piscivore", "Vertivore", "InsVertivore", "Omnivore", "Frugivore", "Nectarivore", "Glucivore", "Herbivore", "Generalist"
            )', 
         'u=list(
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
          )', 
         'y=c(
          "vs_HLornAna3", "vs_HLtacAcu1", "vs_HLdidVir1", "vs_HLgymLea1", "vs_HLpseCup1", 
          "vs_HLmyrTri1", "vs_HLchoDid1", "vs_HLchoHof3", "vs_HLproCap3", 
          "vs_HLpanLeo1", "LionClade", "vs_HLpanOnc1", "vs_HLaciJub2", "CheetahClade", "vs_HLursThi1", 
          "vs_ursMar1", "vs_HLursArc1", "vs_lepWed1", "vs_HLailMel2", 
          "vs_HLmirAng2", "vs_HLphoVit1", "HarborSealGreySealClade", "vs_HLeriBar1", 
          "vs_HLodoRos1", "vs_HLcalUrs1", "vs_HLzalCal1", 
          "vs_HLmelCap1", "vs_HLgulGul1", "vs_HLneoVis1", 
          "vs_HLpteBra1", "vs_HLlutLut1", "vs_enhLutKen1", 
          "vs_HLlycPic2", "vs_HLgloMel1", "vs_HLpepEle1", 
          "vs_HLturAdu1", "DolphinClade", "vs_orcOrc1", "vs_HLescRob1", "vs_HLlniGeo1", "amazonRiverDolphinFromYeast",
          "vs_HLbalEde1", "vs_HLmegNov1", "vs_HLcynGun1", 
          "vs_HLmerUng1", "vs_HLeulMon1", "vs_HLeulFul1", 
          "vs_eulMac1", "vs_HLeulFla1", "vs_ponAbe3", 
          "vs_panTro6", "vs_HLrhiRox2", "LangurClade", "vs_HLallNig1", 
          "vs_HLeryPat1", "vs_chlSab2", "geunonClade", "vs_HLtheGel1", 
          "vs_HLpapAnu5", "vs_HLmanSph1", "DrillMandrillClade", "vs_cerAty1", 
          "vs_HLtraJav1", "vs_mm10", "vs_HLmarFla1", "marmotClade"
         )')
{
args = c('r=CategoricalPrunedCarnivoreTree', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 
         'c=c("C-Invertebrate-eater", "C-Endotherm-Carnivore", "C-Herpetivore", "C-Piscivore", "C-Nonspecific-Vertebrate-eater", "C-Scavenger", 
              "O-For Examination", "O-Scavenger", 
              "H-Frugivore", "H-Nectarivore", "H-Granivore", "H-Nonspecific-Herbivore", 
              "C-Terrestrial-vertebrates-eater", "C-All-vertebrate-eater", "C-All-Animals-Eater", 
              "H-High-sugar-plants-Eater", "H-Low-sugar-plants-Eater", "H-All-plants-Eater", 
              "O-Generalist", 
              "C-InsVertivore-Mixed", "C-InsVertivore-Piscivore", "C-InsVertivore-Insectivore","C-InsVertivore-Carnivore",
              "Insectivore", "Herpetivore", "Piscivore", "Vertivore", "InsVertivore", "Omnivore", "Frugivore", "Nectarivore", "Glucivore", "Herbivore", "Generalist"
            )', 
         'u=list(
            c("C-Invertebrate-eater", "Carnivore"), c("C-InsVertivore-Insectivore", "Carnivore"),
            c("C-Herpetivore", "Carnivore"),
            c("C-Piscivore", "Carnivore"), c("C-InsVertivore-Piscivore", "Carnivore"),
            c("C-Endotherm-Carnivore", "Carnivore"), c("C-Scavenger", "Carnivore"), c("C-Nonspecific-Vertebrate-eater", "Carnivore"),
            c("C-Terrestrial-vertebrates-eater", "Carnivore"), c("C-All-vertebrate-eater", "Carnivore"), c("C-InsVertivore-Carnivore", "Carnivore"),
            c("C-InsVertivore-Mixed", "Omnivore"), 
            c("O-For Examination", "Omnivore"), c("O-Scavenger", "Omnivore"),
            c("H-Frugivore", "Herbivore"), 
            c("H-Nectarivore", "Herbivore"), 
            c("H-High-sugar-plants-Eater", "Herbivore"),
            c("H-Granivore", "Herbivore"), c("H-Nonspecific-Herbivore", "Herbivore"), 
            c("H-Low-sugar-plants-Eater", "Herbivore"), c("H-All-plants-Eater", "Herbivore"),
            c("O-Generalist", "Omnivore")
          )', 
         'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01', 
         'y=c(
          "vs_HLornAna3", "vs_HLtacAcu1", "vs_HLdidVir1", "vs_HLgymLea1", "vs_HLpseCup1", 
          "vs_HLmyrTri1", "vs_HLchoDid1", "vs_HLchoHof3", "vs_HLproCap3", 
          "vs_HLpanLeo1", "LionClade", "vs_HLpanOnc1", "vs_HLaciJub2", "CheetahClade", "vs_HLursThi1", 
          "vs_ursMar1", "vs_HLursArc1", "vs_lepWed1", "vs_HLailMel2", 
          "vs_HLmirAng2", "vs_HLphoVit1", "HarborSealGreySealClade", "vs_HLeriBar1", 
          "vs_HLodoRos1", "vs_HLcalUrs1", "vs_HLzalCal1", 
          "vs_HLmelCap1", "vs_HLgulGul1", "vs_HLneoVis1", 
          "vs_HLpteBra1", "vs_HLlutLut1", "vs_enhLutKen1", 
          "vs_HLlycPic2", "vs_HLgloMel1", "vs_HLpepEle1", 
          "vs_HLturAdu1", "DolphinClade", "vs_orcOrc1", "vs_HLescRob1", "vs_HLlniGeo1", "amazonRiverDolphinFromYeast",
          "vs_HLbalEde1", "vs_HLmegNov1", "vs_HLcynGun1", 
          "vs_HLmerUng1", "vs_HLeulMon1", "vs_HLeulFul1", 
          "vs_eulMac1", "vs_HLeulFla1", "vs_ponAbe3", 
          "vs_panTro6", "vs_HLrhiRox2", "LangurClade", "vs_HLallNig1", 
          "vs_HLeryPat1", "vs_chlSab2", "geunonClade", "vs_HLtheGel1", 
          "vs_HLpapAnu5", "vs_HLmanSph1", "DrillMandrillClade", "vs_cerAty1", 
          "vs_HLtraJav1", "vs_mm10", "vs_HLmarFla1", "marmotClade"
         )',
         'p=c(
          "vs_HLellTal1", "vs_HLellLut1", "vs_HLarvAmp1", "voleClade",
          "vs_HLmusSpi1", "vs_HLmusCar1", "vs_HLmasCou1", "vs_HLmusPah1", "mouseClade",
          "vs_HLhysCri1", "vs_HLthrSwi1", "vs_HLpetTyp1", "vs_hetGla2", "vs_chiLan1", "vs_HLdinBra1", "vs_HLcteSoc1", "vs_octDeg1", "vs_HLcoePre1", "vs_HLdasPun1", "vs_HLdolPat1", "gundiGuineaPigClade",
          "vs_HLoryGaz1", "vs_HLbeaHun1", "vs_HLkobLecLec1", "vs_HLkobLecLec1", "vs_HLmadKir1", "vs_HLneoPyg1", "vs_HLphiMax1", "vs_HLoreOre1", "vs_HLneoMos1", "vs_HLaepMel1", "vs_HLtraImb1", "Bovidae",
          "vs_HLhydIne1", "vs_HLmunMun1", "Cervidae",
          "vs_HLtraKan1", "mouseDeerOtherIsKept",
          "vs_HLmurAurFea1", "outerVespert",
          "vs_HLmyoLuc1", "Nearctic",
          "vs_myoDav1", "Myotis",
          "vs_HLpipPip1", "vs_HLlasBor1", "vs_HLnycHum2", "Vespertilioninae",
          "vs_HLmacSob1", "FoxLongTounge",
          "vs_HLeidHel2", "outerPeropodidae",
          "vs_HLeonSpe1", "Roussetinae"
         )')
}
{
args = c('r=CategoricalBinaryInsectivoreTree', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 
         'c=c(
              "C-Invertebrate-eater", "C-Endotherm-Carnivore", "C-Herpetivore", "C-Piscivore", "C-Nonspecific-Vertebrate-eater", "C-Scavenger", 
              "O-For Examination", "O-Scavenger", 
              "H-Frugivore", "H-Nectarivore", "H-Granivore", "H-Nonspecific-Herbivore", 
              "C-Terrestrial-vertebrates-eater", "C-All-vertebrate-eater", "C-All-Animals-Eater", 
              "H-High-sugar-plants-Eater", "H-Low-sugar-plants-Eater", "H-All-plants-Eater", 
              "O-Generalist", 
              "C-InsVertivore-Mixed", "C-InsVertivore-Piscivore", "C-InsVertivore-Insectivore","C-InsVertivore-Carnivore",
              "Insectivore", "Herpetivore", "Piscivore", "Vertivore", "InsVertivore", "Omnivore", "Frugivore", "Nectarivore", "Glucivore", "Herbivore", "Generalist"
            )', 
         'u=list(
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
            c("O-Generalist", "Omnivore"),
            c("Omnivore", "Background"), c("Vertivore", "Background"), c("Herbivore", "Background")
          )', 
         'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01',
         'y=c(
          "vs_HLornAna3", "vs_HLtacAcu1", "vs_HLdidVir1", "vs_HLgymLea1", "vs_HLpseCup1", 
          "vs_HLmyrTri1", "vs_HLchoDid1", "vs_HLchoHof3", "vs_HLproCap3", 
          "vs_HLpanLeo1", "LionClade", "vs_HLpanOnc1", "vs_HLaciJub2", "CheetahClade", "vs_HLursThi1", 
          "vs_ursMar1", "vs_HLursArc1", "vs_lepWed1", "vs_HLailMel2", 
          "vs_HLmirAng2", "vs_HLphoVit1", "HarborSealGreySealClade", "vs_HLeriBar1", 
          "vs_HLodoRos1", "vs_HLcalUrs1", "vs_HLzalCal1", 
          "vs_HLmelCap1", "vs_HLgulGul1", "vs_HLneoVis1", 
          "vs_HLpteBra1", "vs_HLlutLut1", "vs_enhLutKen1", 
          "vs_HLlycPic2", "vs_HLgloMel1", "vs_HLpepEle1", 
          "vs_HLturAdu1", "DolphinClade", "vs_orcOrc1", "vs_HLescRob1", "vs_HLlniGeo1", "amazonRiverDolphinFromYeast",
          "vs_HLbalEde1", "vs_HLmegNov1", "vs_HLcynGun1", 
          "vs_HLmerUng1", "vs_HLeulMon1", "vs_HLeulFul1", 
          "vs_eulMac1", "vs_HLeulFla1", "vs_ponAbe3", 
          "vs_panTro6", "vs_HLrhiRox2", "LangurClade", "vs_HLallNig1", 
          "vs_HLeryPat1", "vs_chlSab2", "geunonClade", "vs_HLtheGel1", 
          "vs_HLpapAnu5", "vs_HLmanSph1", "DrillMandrillClade", "vs_cerAty1", 
          "vs_HLtraJav1", "vs_mm10", "vs_HLmarFla1", "marmotClade"
         )',
         'p=c(
          "vs_HLellTal1", "vs_HLellLut1", "vs_HLarvAmp1", "voleClade",
          "vs_HLmusSpi1", "vs_HLmusCar1", "vs_HLmasCou1", "vs_HLmusPah1", "mouseClade",
          "vs_HLhysCri1", "vs_HLthrSwi1", "vs_HLpetTyp1", "vs_hetGla2", "vs_chiLan1", "vs_HLdinBra1", "vs_HLcteSoc1", "vs_octDeg1", "vs_HLcoePre1", "vs_HLdasPun1", "vs_HLdolPat1", "gundiGuineaPigClade",
          "vs_HLoryGaz1", "vs_HLbeaHun1", "vs_HLkobLecLec1", "vs_HLkobLecLec1", "vs_HLmadKir1", "vs_HLneoPyg1", "vs_HLphiMax1", "vs_HLoreOre1", "vs_HLneoMos1", "vs_HLaepMel1", "vs_HLtraImb1", "Bovidae",
          "vs_HLhydIne1", "vs_HLmunMun1", "Cervidae",
          "vs_HLtraKan1", "mouseDeerOtherIsKept",
          "vs_HLmurAurFea1", "outerVespert",
          "vs_HLmyoLuc1", "Nearctic",
          "vs_myoDav1", "Myotis",
          "vs_HLpipPip1", "vs_HLlasBor1", "vs_HLnycHum2", "Vespertilioninae",
          "vs_HLmacSob1", "FoxLongTounge",
          "vs_HLeidHel2", "outerPeropodidae",
          "vs_HLeonSpe1", "Roussetinae"
         )')
  args = c('r=CategoricalBinaryVertivoreTree', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 
           'c=c(
              "C-Invertebrate-eater", "C-Endotherm-Carnivore", "C-Herpetivore", "C-Piscivore", "C-Nonspecific-Vertebrate-eater", "C-Scavenger", 
              "O-For Examination", "O-Scavenger", 
              "H-Frugivore", "H-Nectarivore", "H-Granivore", "H-Nonspecific-Herbivore", 
              "C-Terrestrial-vertebrates-eater", "C-All-vertebrate-eater", "C-All-Animals-Eater", 
              "H-High-sugar-plants-Eater", "H-Low-sugar-plants-Eater", "H-All-plants-Eater", 
              "O-Generalist", 
              "C-InsVertivore-Mixed", "C-InsVertivore-Piscivore", "C-InsVertivore-Insectivore","C-InsVertivore-Carnivore",
              "Insectivore", "Herpetivore", "Piscivore", "Vertivore", "InsVertivore", "Omnivore", "Frugivore", "Nectarivore", "Glucivore", "Herbivore", "Generalist"
            )', 
           'u=list(
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
            c("O-Generalist", "Omnivore"),
            c("Omnivore", "Background"), c("Insectivore", "Background"), c("Herbivore", "Background")
          )', 
          'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01',
          'y=c(
          "vs_HLornAna3", "vs_HLtacAcu1", "vs_HLdidVir1", "vs_HLgymLea1", "vs_HLpseCup1", 
          "vs_HLmyrTri1", "vs_HLchoDid1", "vs_HLchoHof3", "vs_HLproCap3", 
          "vs_HLpanLeo1", "LionClade", "vs_HLpanOnc1", "vs_HLaciJub2", "CheetahClade", "vs_HLursThi1", 
          "vs_ursMar1", "vs_HLursArc1", "vs_lepWed1", "vs_HLailMel2", 
          "vs_HLmirAng2", "vs_HLphoVit1", "HarborSealGreySealClade", "vs_HLeriBar1", 
          "vs_HLodoRos1", "vs_HLcalUrs1", "vs_HLzalCal1", 
          "vs_HLmelCap1", "vs_HLgulGul1", "vs_HLneoVis1", 
          "vs_HLpteBra1", "vs_HLlutLut1", "vs_enhLutKen1", 
          "vs_HLlycPic2", "vs_HLgloMel1", "vs_HLpepEle1", 
          "vs_HLturAdu1", "DolphinClade", "vs_orcOrc1", "vs_HLescRob1", "vs_HLlniGeo1", "amazonRiverDolphinFromYeast",
          "vs_HLbalEde1", "vs_HLmegNov1", "vs_HLcynGun1", 
          "vs_HLmerUng1", "vs_HLeulMon1", "vs_HLeulFul1", 
          "vs_eulMac1", "vs_HLeulFla1", "vs_ponAbe3", 
          "vs_panTro6", "vs_HLrhiRox2", "LangurClade", "vs_HLallNig1", 
          "vs_HLeryPat1", "vs_chlSab2", "geunonClade", "vs_HLtheGel1", 
          "vs_HLpapAnu5", "vs_HLmanSph1", "DrillMandrillClade", "vs_cerAty1", 
          "vs_HLtraJav1", "vs_mm10", "vs_HLmarFla1", "marmotClade"
         )',
         'p=c(
          "vs_HLellTal1", "vs_HLellLut1", "vs_HLarvAmp1", "voleClade",
          "vs_HLmusSpi1", "vs_HLmusCar1", "vs_HLmasCou1", "vs_HLmusPah1", "mouseClade",
          "vs_HLhysCri1", "vs_HLthrSwi1", "vs_HLpetTyp1", "vs_hetGla2", "vs_chiLan1", "vs_HLdinBra1", "vs_HLcteSoc1", "vs_octDeg1", "vs_HLcoePre1", "vs_HLdasPun1", "vs_HLdolPat1", "gundiGuineaPigClade",
          "vs_HLoryGaz1", "vs_HLbeaHun1", "vs_HLkobLecLec1", "vs_HLkobLecLec1", "vs_HLmadKir1", "vs_HLneoPyg1", "vs_HLphiMax1", "vs_HLoreOre1", "vs_HLneoMos1", "vs_HLaepMel1", "vs_HLtraImb1", "Bovidae",
          "vs_HLhydIne1", "vs_HLmunMun1", "Cervidae",
          "vs_HLtraKan1", "mouseDeerOtherIsKept",
          "vs_HLmurAurFea1", "outerVespert",
          "vs_HLmyoLuc1", "Nearctic",
          "vs_myoDav1", "Myotis",
          "vs_HLpipPip1", "vs_HLlasBor1", "vs_HLnycHum2", "Vespertilioninae",
          "vs_HLmacSob1", "FoxLongTounge",
          "vs_HLeidHel2", "outerPeropodidae",
          "vs_HLeonSpe1", "Roussetinae"
         )')
  args = c('r=CategoricalBinaryHerbivoreTree', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 
           'c=c(
              "C-Invertebrate-eater", "C-Endotherm-Carnivore", "C-Herpetivore", "C-Piscivore", "C-Nonspecific-Vertebrate-eater", "C-Scavenger", 
              "O-For Examination", "O-Scavenger", 
              "H-Frugivore", "H-Nectarivore", "H-Granivore", "H-Nonspecific-Herbivore", 
              "C-Terrestrial-vertebrates-eater", "C-All-vertebrate-eater", "C-All-Animals-Eater", 
              "H-High-sugar-plants-Eater", "H-Low-sugar-plants-Eater", "H-All-plants-Eater", 
              "O-Generalist", 
              "C-InsVertivore-Mixed", "C-InsVertivore-Piscivore", "C-InsVertivore-Insectivore","C-InsVertivore-Carnivore",
              "Insectivore", "Herpetivore", "Piscivore", "Vertivore", "InsVertivore", "Omnivore", "Frugivore", "Nectarivore", "Glucivore", "Herbivore", "Generalist"
            )', 
           'u=list(
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
            c("O-Generalist", "Omnivore"),
            c("Omnivore", "Background"), c("Insectivore", "Background"), c("Vertivore", "Background")
          )', 
          'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01',
          'y=c(
          "vs_HLornAna3", "vs_HLtacAcu1", "vs_HLdidVir1", "vs_HLgymLea1", "vs_HLpseCup1", 
          "vs_HLmyrTri1", "vs_HLchoDid1", "vs_HLchoHof3", "vs_HLproCap3", 
          "vs_HLpanLeo1", "LionClade", "vs_HLpanOnc1", "vs_HLaciJub2", "CheetahClade", "vs_HLursThi1", 
          "vs_ursMar1", "vs_HLursArc1", "vs_lepWed1", "vs_HLailMel2", 
          "vs_HLmirAng2", "vs_HLphoVit1", "HarborSealGreySealClade", "vs_HLeriBar1", 
          "vs_HLodoRos1", "vs_HLcalUrs1", "vs_HLzalCal1", 
          "vs_HLmelCap1", "vs_HLgulGul1", "vs_HLneoVis1", 
          "vs_HLpteBra1", "vs_HLlutLut1", "vs_enhLutKen1", 
          "vs_HLlycPic2", "vs_HLgloMel1", "vs_HLpepEle1", 
          "vs_HLturAdu1", "DolphinClade", "vs_orcOrc1", "vs_HLescRob1", "vs_HLlniGeo1", "amazonRiverDolphinFromYeast",
          "vs_HLbalEde1", "vs_HLmegNov1", "vs_HLcynGun1", 
          "vs_HLmerUng1", "vs_HLeulMon1", "vs_HLeulFul1", 
          "vs_eulMac1", "vs_HLeulFla1", "vs_ponAbe3", 
          "vs_panTro6", "vs_HLrhiRox2", "LangurClade", "vs_HLallNig1", 
          "vs_HLeryPat1", "vs_chlSab2", "geunonClade", "vs_HLtheGel1", 
          "vs_HLpapAnu5", "vs_HLmanSph1", "DrillMandrillClade", "vs_cerAty1", 
          "vs_HLtraJav1", "vs_mm10", "vs_HLmarFla1", "marmotClade"
         )',
         'p=c(
          "vs_HLellTal1", "vs_HLellLut1", "vs_HLarvAmp1", "voleClade",
          "vs_HLmusSpi1", "vs_HLmusCar1", "vs_HLmasCou1", "vs_HLmusPah1", "mouseClade",
          "vs_HLhysCri1", "vs_HLthrSwi1", "vs_HLpetTyp1", "vs_hetGla2", "vs_chiLan1", "vs_HLdinBra1", "vs_HLcteSoc1", "vs_octDeg1", "vs_HLcoePre1", "vs_HLdasPun1", "vs_HLdolPat1", "gundiGuineaPigClade",
          "vs_HLoryGaz1", "vs_HLbeaHun1", "vs_HLkobLecLec1", "vs_HLkobLecLec1", "vs_HLmadKir1", "vs_HLneoPyg1", "vs_HLphiMax1", "vs_HLoreOre1", "vs_HLneoMos1", "vs_HLaepMel1", "vs_HLtraImb1", "Bovidae",
          "vs_HLhydIne1", "vs_HLmunMun1", "Cervidae",
          "vs_HLtraKan1", "mouseDeerOtherIsKept",
          "vs_HLmurAurFea1", "outerVespert",
          "vs_HLmyoLuc1", "Nearctic",
          "vs_myoDav1", "Myotis",
          "vs_HLpipPip1", "vs_HLlasBor1", "vs_HLnycHum2", "Vespertilioninae",
          "vs_HLmacSob1", "FoxLongTounge",
          "vs_HLeidHel2", "outerPeropodidae",
          "vs_HLeonSpe1", "Roussetinae"
         )')
  args = c('r=CategoricalBinaryOmnivoreTree', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 
           'c=c(
              "C-Invertebrate-eater", "C-Endotherm-Carnivore", "C-Herpetivore", "C-Piscivore", "C-Nonspecific-Vertebrate-eater", "C-Scavenger", 
              "O-For Examination", "O-Scavenger", 
              "H-Frugivore", "H-Nectarivore", "H-Granivore", "H-Nonspecific-Herbivore", 
              "C-Terrestrial-vertebrates-eater", "C-All-vertebrate-eater", "C-All-Animals-Eater", 
              "H-High-sugar-plants-Eater", "H-Low-sugar-plants-Eater", "H-All-plants-Eater", 
              "O-Generalist", 
              "C-InsVertivore-Mixed", "C-InsVertivore-Piscivore", "C-InsVertivore-Insectivore","C-InsVertivore-Carnivore",
              "Insectivore", "Herpetivore", "Piscivore", "Vertivore", "InsVertivore", "Omnivore", "Frugivore", "Nectarivore", "Glucivore", "Herbivore", "Generalist"
            )', 
           'u=list(
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
            c("O-Generalist", "Omnivore"),
            c("Herbivore", "Background"), c("Insectivore", "Background"), c("Vertivore", "Background")
          )', 
          'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01',
          'y=c(
          "vs_HLornAna3", "vs_HLtacAcu1", "vs_HLdidVir1", "vs_HLgymLea1", "vs_HLpseCup1", 
          "vs_HLmyrTri1", "vs_HLchoDid1", "vs_HLchoHof3", "vs_HLproCap3", 
          "vs_HLpanLeo1", "LionClade", "vs_HLpanOnc1", "vs_HLaciJub2", "CheetahClade", "vs_HLursThi1", 
          "vs_ursMar1", "vs_HLursArc1", "vs_lepWed1", "vs_HLailMel2", 
          "vs_HLmirAng2", "vs_HLphoVit1", "HarborSealGreySealClade", "vs_HLeriBar1", 
          "vs_HLodoRos1", "vs_HLcalUrs1", "vs_HLzalCal1", 
          "vs_HLmelCap1", "vs_HLgulGul1", "vs_HLneoVis1", 
          "vs_HLpteBra1", "vs_HLlutLut1", "vs_enhLutKen1", 
          "vs_HLlycPic2", "vs_HLgloMel1", "vs_HLpepEle1", 
          "vs_HLturAdu1", "DolphinClade", "vs_orcOrc1", "vs_HLescRob1", "vs_HLlniGeo1", "amazonRiverDolphinFromYeast",
          "vs_HLbalEde1", "vs_HLmegNov1", "vs_HLcynGun1", 
          "vs_HLmerUng1", "vs_HLeulMon1", "vs_HLeulFul1", 
          "vs_eulMac1", "vs_HLeulFla1", "vs_ponAbe3", 
          "vs_panTro6", "vs_HLrhiRox2", "LangurClade", "vs_HLallNig1", 
          "vs_HLeryPat1", "vs_chlSab2", "geunonClade", "vs_HLtheGel1", 
          "vs_HLpapAnu5", "vs_HLmanSph1", "DrillMandrillClade", "vs_cerAty1", 
          "vs_HLtraJav1", "vs_mm10", "vs_HLmarFla1", "marmotClade"
         )',
         'p=c(
          "vs_HLellTal1", "vs_HLellLut1", "vs_HLarvAmp1", "voleClade",
          "vs_HLmusSpi1", "vs_HLmusCar1", "vs_HLmasCou1", "vs_HLmusPah1", "mouseClade",
          "vs_HLhysCri1", "vs_HLthrSwi1", "vs_HLpetTyp1", "vs_hetGla2", "vs_chiLan1", "vs_HLdinBra1", "vs_HLcteSoc1", "vs_octDeg1", "vs_HLcoePre1", "vs_HLdasPun1", "vs_HLdolPat1", "gundiGuineaPigClade",
          "vs_HLoryGaz1", "vs_HLbeaHun1", "vs_HLkobLecLec1", "vs_HLkobLecLec1", "vs_HLmadKir1", "vs_HLneoPyg1", "vs_HLphiMax1", "vs_HLoreOre1", "vs_HLneoMos1", "vs_HLaepMel1", "vs_HLtraImb1", "Bovidae",
          "vs_HLhydIne1", "vs_HLmunMun1", "Cervidae",
          "vs_HLtraKan1", "mouseDeerOtherIsKept",
          "vs_HLmurAurFea1", "outerVespert",
          "vs_HLmyoLuc1", "Nearctic",
          "vs_myoDav1", "Myotis",
          "vs_HLpipPip1", "vs_HLlasBor1", "vs_HLnycHum2", "Vespertilioninae",
          "vs_HLmacSob1", "FoxLongTounge",
          "vs_HLeidHel2", "outerPeropodidae",
          "vs_HLeonSpe1", "Roussetinae"
         )')
  args = c('r=CategoricalBinaryCarnivoreTree', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 
           'c=c(
              "C-Invertebrate-eater", "C-Endotherm-Carnivore", "C-Herpetivore", "C-Piscivore", "C-Nonspecific-Vertebrate-eater", "C-Scavenger", 
              "O-For Examination", "O-Scavenger", 
              "H-Frugivore", "H-Nectarivore", "H-Granivore", "H-Nonspecific-Herbivore", 
              "C-Terrestrial-vertebrates-eater", "C-All-vertebrate-eater", "C-All-Animals-Eater", 
              "H-High-sugar-plants-Eater", "H-Low-sugar-plants-Eater", "H-All-plants-Eater", 
              "O-Generalist", 
              "C-InsVertivore-Mixed", "C-InsVertivore-Piscivore", "C-InsVertivore-Insectivore","C-InsVertivore-Carnivore",
              "Insectivore", "Herpetivore", "Piscivore", "Vertivore", "InsVertivore", "Omnivore", "Frugivore", "Nectarivore", "Glucivore", "Herbivore", "Generalist", "Carnivore"
            )', 
           'u=list(
            c("C-Invertebrate-eater", "Insectivore"), c("C-InsVertivore-Insectivore", "Insectivore"),
            c("C-Herpetivore", "Vertivore"),
            c("C-Piscivore", "Vertivore"), c("C-InsVertivore-Piscivore", "Vertivore"),
            c("C-Endotherm-Carnivore", "Vertivore"), c("C-Scavenger", "Vertivore"), c("C-Nonspecific-Vertebrate-eater", "Vertivore"),
            c("C-Terrestrial-vertebrates-eater", "Vertivore"), c("C-All-vertebrate-eater", "Vertivore"), c("C-InsVertivore-Carnivore", "Vertivore"),
            c("C-InsVertivore-Mixed", "Mixed"), 
            c("O-For Examination", "Omnivore"), c("O-Scavenger", "Omnivore"),
            c("H-Frugivore", "Herbivore"), 
            c("H-Nectarivore", "Herbivore"), 
            c("H-High-sugar-plants-Eater", "Herbivore"),
            c("H-Granivore", "Herbivore"), c("H-Nonspecific-Herbivore", "Herbivore"), 
            c("H-Low-sugar-plants-Eater", "Herbivore"), c("H-All-plants-Eater", "Herbivore"),
            c("O-Generalist", "Omnivore"),
            c("Vertivore", "Carnivore"), c("Insectivore", "Carnivore"), c("Mixed", "Carnivore"), c("Omnivore", "Background"), c("Herbivore", "Background")
          )', 
          'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01',
          'y=c(
          "vs_HLornAna3", "vs_HLtacAcu1", "vs_HLdidVir1", "vs_HLgymLea1", "vs_HLpseCup1", 
          "vs_HLmyrTri1", "vs_HLchoDid1", "vs_HLchoHof3", "vs_HLproCap3", 
          "vs_HLpanLeo1", "LionClade", "vs_HLpanOnc1", "vs_HLaciJub2", "CheetahClade", "vs_HLursThi1", 
          "vs_ursMar1", "vs_HLursArc1", "vs_lepWed1", "vs_HLailMel2", 
          "vs_HLmirAng2", "vs_HLphoVit1", "HarborSealGreySealClade", "vs_HLeriBar1", 
          "vs_HLodoRos1", "vs_HLcalUrs1", "vs_HLzalCal1", 
          "vs_HLmelCap1", "vs_HLgulGul1", "vs_HLneoVis1", 
          "vs_HLpteBra1", "vs_HLlutLut1", "vs_enhLutKen1", 
          "vs_HLlycPic2", "vs_HLgloMel1", "vs_HLpepEle1", 
          "vs_HLturAdu1", "DolphinClade", "vs_orcOrc1", "vs_HLescRob1", "vs_HLlniGeo1", "amazonRiverDolphinFromYeast",
          "vs_HLbalEde1", "vs_HLmegNov1", "vs_HLcynGun1", 
          "vs_HLmerUng1", "vs_HLeulMon1", "vs_HLeulFul1", 
          "vs_eulMac1", "vs_HLeulFla1", "vs_ponAbe3", 
          "vs_panTro6", "vs_HLrhiRox2", "LangurClade", "vs_HLallNig1", 
          "vs_HLeryPat1", "vs_chlSab2", "geunonClade", "vs_HLtheGel1", 
          "vs_HLpapAnu5", "vs_HLmanSph1", "DrillMandrillClade", "vs_cerAty1", 
          "vs_HLtraJav1", "vs_mm10", "vs_HLmarFla1", "marmotClade"
         )',
         'p=c(
          "vs_HLellTal1", "vs_HLellLut1", "vs_HLarvAmp1", "voleClade",
          "vs_HLmusSpi1", "vs_HLmusCar1", "vs_HLmasCou1", "vs_HLmusPah1", "mouseClade",
          "vs_HLhysCri1", "vs_HLthrSwi1", "vs_HLpetTyp1", "vs_hetGla2", "vs_chiLan1", "vs_HLdinBra1", "vs_HLcteSoc1", "vs_octDeg1", "vs_HLcoePre1", "vs_HLdasPun1", "vs_HLdolPat1", "gundiGuineaPigClade",
          "vs_HLoryGaz1", "vs_HLbeaHun1", "vs_HLkobLecLec1", "vs_HLkobLecLec1", "vs_HLmadKir1", "vs_HLneoPyg1", "vs_HLphiMax1", "vs_HLoreOre1", "vs_HLneoMos1", "vs_HLaepMel1", "vs_HLtraImb1", "Bovidae",
          "vs_HLhydIne1", "vs_HLmunMun1", "Cervidae",
          "vs_HLtraKan1", "mouseDeerOtherIsKept",
          "vs_HLmurAurFea1", "outerVespert",
          "vs_HLmyoLuc1", "Nearctic",
          "vs_myoDav1", "Myotis",
          "vs_HLpipPip1", "vs_HLlasBor1", "vs_HLnycHum2", "Vespertilioninae",
          "vs_HLmacSob1", "FoxLongTounge",
          "vs_HLeidHel2", "outerPeropodidae",
          "vs_HLeonSpe1", "Roussetinae"
         )')
  
}
{
args = c('r=CategoricalNoMegabranchInsvertTree', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 
         'c=c(
              "C-Invertebrate-eater", "C-Endotherm-Carnivore", "C-Herpetivore", "C-Piscivore", "C-Nonspecific-Vertebrate-eater", "C-Scavenger", 
              "O-For Examination", "O-Scavenger", 
              "H-Frugivore", "H-Nectarivore", "H-Granivore", "H-Nonspecific-Herbivore", 
              "C-Terrestrial-vertebrates-eater", "C-All-vertebrate-eater", "C-All-Animals-Eater", 
              "H-High-sugar-plants-Eater", "H-Low-sugar-plants-Eater", "H-All-plants-Eater", 
              "O-Generalist", 
              "C-InsVertivore-Mixed", "C-InsVertivore-Piscivore", "C-InsVertivore-Insectivore","C-InsVertivore-Carnivore",
              "Insectivore", "Herpetivore", "Piscivore", "Vertivore", "InsVertivore", "Omnivore", "Frugivore", "Nectarivore", "Glucivore", "Herbivore", "Generalist"
            )', 
         'u=list(
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
          )', 
         'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01',
         'y=c(
          "vs_HLornAna3", "vs_HLtacAcu1", "vs_HLdidVir1", "vs_HLgymLea1", "vs_HLpseCup1", 
          "vs_HLmyrTri1", "vs_HLchoDid1", "vs_HLchoHof3", "vs_HLproCap3", 
          "vs_HLpanLeo1", "LionClade", "vs_HLpanOnc1", "vs_HLaciJub2", "CheetahClade", "vs_HLursThi1", 
          "vs_ursMar1", "vs_HLursArc1", "vs_lepWed1", "vs_HLailMel2", 
          "vs_HLmirAng2", "vs_HLphoVit1", "HarborSealGreySealClade", "vs_HLeriBar1", 
          "vs_HLodoRos1", "vs_HLcalUrs1", "vs_HLzalCal1", 
          "vs_HLmelCap1", "vs_HLgulGul1", "vs_HLneoVis1", 
          "vs_HLpteBra1", "vs_HLlutLut1", "vs_enhLutKen1", 
          "vs_HLlycPic2", "vs_HLgloMel1", "vs_HLpepEle1", 
          "vs_HLturAdu1", "DolphinClade", "vs_orcOrc1", "vs_HLescRob1", "vs_HLlniGeo1", "amazonRiverDolphinFromYeast",
          "vs_HLbalEde1", "vs_HLmegNov1", "vs_HLcynGun1", 
          "vs_HLmerUng1", "vs_HLeulMon1", "vs_HLeulFul1", 
          "vs_eulMac1", "vs_HLeulFla1", "vs_ponAbe3", 
          "vs_panTro6", "vs_HLrhiRox2", "LangurClade", "vs_HLallNig1", 
          "vs_HLeryPat1", "vs_chlSab2", "geunonClade", "vs_HLtheGel1", 
          "vs_HLpapAnu5", "vs_HLmanSph1", "DrillMandrillClade", "vs_cerAty1", 
          "vs_HLtraJav1", "vs_mm10", "vs_HLmarFla1", "marmotClade"
         )',
         'p=c(
          "vs_HLellTal1", "vs_HLellLut1", "vs_HLarvAmp1", "voleClade",
          "vs_HLmusSpi1", "vs_HLmusCar1", "vs_HLmasCou1", "vs_HLmusPah1", "mouseClade",
          "vs_HLhysCri1", "vs_HLthrSwi1", "vs_HLpetTyp1", "vs_hetGla2", "vs_chiLan1", "vs_HLdinBra1", "vs_HLcteSoc1", "vs_octDeg1", "vs_HLcoePre1", "vs_HLdasPun1", "vs_HLdolPat1", "gundiGuineaPigClade",
          "vs_HLoryGaz1", "vs_HLbeaHun1", "vs_HLkobLecLec1", "vs_HLkobLecLec1", "vs_HLmadKir1", "vs_HLneoPyg1", "vs_HLphiMax1", "vs_HLoreOre1", "vs_HLneoMos1", "vs_HLaepMel1", "vs_HLtraImb1", "Bovidae",
          "vs_HLhydIne1", "vs_HLmunMun1", "Cervidae",
          "vs_HLtraKan1", "mouseDeerOtherIsKept",
          "vs_HLmurAurFea1", "outerVespert",
          "vs_HLmyoLuc1", "Nearctic",
          "vs_myoDav1", "Myotis",
          "vs_HLpipPip1", "vs_HLlasBor1", "vs_HLnycHum2", "Vespertilioninae",
          "vs_HLmacSob1", "FoxLongTounge",
          "vs_HLeidHel2", "outerPeropodidae",
          "vs_HLeonSpe1", "Roussetinae"
         )')

args = c('r=CategoricalNoMegabranchInsvertTree', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 
         'c=c(
              "C-Invertebrate-eater", "C-Endotherm-Carnivore", "C-Herpetivore", "C-Piscivore", "C-Nonspecific-Vertebrate-eater", "C-Scavenger", 
              "O-For Examination", "O-Scavenger", 
              "H-Frugivore", "H-Nectarivore", "H-Granivore", "H-Nonspecific-Herbivore", 
              "C-Terrestrial-vertebrates-eater", "C-All-vertebrate-eater", "C-All-Animals-Eater", 
              "H-High-sugar-plants-Eater", "H-Low-sugar-plants-Eater", "H-All-plants-Eater", 
              "O-Generalist", 
              "C-InsVertivore-Mixed", "C-InsVertivore-Piscivore", "C-InsVertivore-Insectivore","C-InsVertivore-Carnivore",
              "Insectivore", "Herpetivore", "Piscivore", "Vertivore", "InsVertivore", "Omnivore", "Frugivore", "Nectarivore", "Glucivore", "Herbivore", "Generalist"
            )', 
         'u=list(
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
          )', 
         'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01',
         'y=c(
          "vs_HLornAna3", "vs_HLtacAcu1", "vs_HLdidVir1", "vs_HLgymLea1", "vs_HLpseCup1", 
          "vs_HLmyrTri1", "vs_HLchoDid1", "vs_HLchoHof3", "vs_HLproCap3", 
          "vs_HLpanLeo1", "LionClade", "vs_HLpanOnc1", "vs_HLaciJub2", "CheetahClade", "vs_HLursThi1", 
          "vs_ursMar1", "vs_HLursArc1", "vs_lepWed1", "vs_HLailMel2", 
          "vs_HLmirAng2", "vs_HLphoVit1", "HarborSealGreySealClade", "vs_HLeriBar1", 
          "vs_HLodoRos1", "vs_HLcalUrs1", "vs_HLzalCal1", 
          "vs_HLmelCap1", "vs_HLgulGul1", "vs_HLneoVis1", 
          "vs_HLpteBra1", "vs_HLlutLut1", "vs_enhLutKen1", 
          "vs_HLlycPic2", "vs_HLgloMel1", "vs_HLpepEle1", 
          "vs_HLturAdu1", "DolphinClade", "vs_orcOrc1", "vs_HLescRob1", "vs_HLlniGeo1", "amazonRiverDolphinFromYeast",
          "vs_HLbalEde1", "vs_HLmegNov1", "vs_HLcynGun1", 
          "vs_HLmerUng1", "vs_HLeulMon1", "vs_HLeulFul1", 
          "vs_eulMac1", "vs_HLeulFla1", "vs_ponAbe3", 
          "vs_panTro6", "vs_HLrhiRox2", "LangurClade", "vs_HLallNig1", 
          "vs_HLeryPat1", "vs_chlSab2", "geunonClade", "vs_HLtheGel1", 
          "vs_HLpapAnu5", "vs_HLmanSph1", "DrillMandrillClade", "vs_cerAty1", 
          "vs_HLtraJav1", "vs_mm10", "vs_HLmarFla1", "marmotClade"
         )',
         'p=c(
          "vs_HLellTal1", "vs_HLellLut1", "vs_HLarvAmp1", "voleClade",
          "vs_HLmusSpi1", "vs_HLmusCar1", "vs_HLmasCou1", "vs_HLmusPah1", "mouseClade",
          "vs_HLhysCri1", "vs_HLthrSwi1", "vs_HLpetTyp1", "vs_hetGla2", "vs_chiLan1", "vs_HLdinBra1", "vs_HLcteSoc1", "vs_octDeg1", "vs_HLcoePre1", "vs_HLdasPun1", "vs_HLdolPat1", "gundiGuineaPigClade",
          "vs_HLoryGaz1", "vs_HLbeaHun1", "vs_HLkobLecLec1", "vs_HLkobLecLec1", "vs_HLmadKir1", "vs_HLneoPyg1", "vs_HLphiMax1", "vs_HLoreOre1", "vs_HLneoMos1", "vs_HLaepMel1", "vs_HLtraImb1", "Bovidae",
          "vs_HLhydIne1", "vs_HLmunMun1", "Cervidae",
          "vs_HLtraKan1", "mouseDeerOtherIsKept",
          "vs_HLmurAurFea1", "outerVespert",
          "vs_HLmyoLuc1", "Nearctic",
          "vs_myoDav1", "Myotis",
          "vs_HLpipPip1", "vs_HLlasBor1", "vs_HLnycHum2", "Vespertilioninae",
          "vs_HLmacSob1", "FoxLongTounge",
          "vs_HLeidHel2", "outerPeropodidae",
          "vs_HLeonSpe1", "Roussetinae"
         )')
args = c('r=CategoricalDownsampledInsvertTree', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 
         'c=c(
              "C-Invertebrate-eater", "C-Endotherm-Carnivore", "C-Herpetivore", "C-Piscivore", "C-Nonspecific-Vertebrate-eater", "C-Scavenger", 
              "O-For Examination", "O-Scavenger", 
              "H-Frugivore", "H-Nectarivore", "H-Granivore", "H-Nonspecific-Herbivore", 
              "C-Terrestrial-vertebrates-eater", "C-All-vertebrate-eater", "C-All-Animals-Eater", 
              "H-High-sugar-plants-Eater", "H-Low-sugar-plants-Eater", "H-All-plants-Eater", 
              "O-Generalist", 
              "C-InsVertivore-Mixed", "C-InsVertivore-Piscivore", "C-InsVertivore-Insectivore","C-InsVertivore-Carnivore",
              "Insectivore", "Herpetivore", "Piscivore", "Vertivore", "InsVertivore", "Omnivore", "Frugivore", "Nectarivore", "Glucivore", "Herbivore", "Generalist"
            )', 
         'u=list(
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
          )', 
         'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01',
         'y=c(
          "vs_HLornAna3", "vs_HLtacAcu1", "vs_HLdidVir1", "vs_HLgymLea1", "vs_HLpseCup1", 
          "vs_HLmyrTri1", "vs_HLchoDid1", "vs_HLchoHof3", "vs_HLproCap3", 
          "vs_HLpanLeo1", "LionClade", "vs_HLpanOnc1", "vs_HLaciJub2", "CheetahClade", "vs_HLursThi1", 
          "vs_ursMar1", "vs_HLursArc1", "vs_lepWed1", "vs_HLailMel2", 
          "vs_HLmirAng2", "vs_HLphoVit1", "HarborSealGreySealClade", "vs_HLeriBar1", 
          "vs_HLodoRos1", "vs_HLcalUrs1", "vs_HLzalCal1", 
          "vs_HLmelCap1", "vs_HLgulGul1", "vs_HLneoVis1", 
          "vs_HLpteBra1", "vs_HLlutLut1", "vs_enhLutKen1", 
          "vs_HLlycPic2", "vs_HLgloMel1", "vs_HLpepEle1", 
          "vs_HLturAdu1", "DolphinClade", "vs_orcOrc1", "vs_HLescRob1", "vs_HLlniGeo1", "amazonRiverDolphinFromYeast",
          "vs_HLbalEde1", "vs_HLmegNov1", "vs_HLcynGun1", 
          "vs_HLmerUng1", "vs_HLeulMon1", "vs_HLeulFul1", 
          "vs_eulMac1", "vs_HLeulFla1", "vs_ponAbe3", 
          "vs_panTro6", "vs_HLrhiRox2", "LangurClade", "vs_HLallNig1", 
          "vs_HLeryPat1", "vs_chlSab2", "geunonClade", "vs_HLtheGel1", 
          "vs_HLpapAnu5", "vs_HLmanSph1", "DrillMandrillClade", "vs_cerAty1", 
          "vs_HLtraJav1", "vs_mm10", "vs_HLmarFla1", "marmotClade"
         )',
         'p=c(
          "vs_HLellTal1", "vs_HLellLut1", "vs_HLarvAmp1", "voleClade",
          "vs_HLmusSpi1", "vs_HLmusCar1", "vs_HLmasCou1", "vs_HLmusPah1", "mouseClade",
          "vs_HLhysCri1", "vs_HLthrSwi1", "vs_HLpetTyp1", "vs_hetGla2", "vs_chiLan1", "vs_HLdinBra1", "vs_HLcteSoc1", "vs_octDeg1", "vs_HLcoePre1", "vs_HLdasPun1", "vs_HLdolPat1", "gundiGuineaPigClade",
          "vs_HLoryGaz1", "vs_HLbeaHun1", "vs_HLkobLecLec1", "vs_HLkobLecLec1", "vs_HLmadKir1", "vs_HLneoPyg1", "vs_HLphiMax1", "vs_HLoreOre1", "vs_HLneoMos1", "vs_HLaepMel1", "vs_HLtraImb1", "Bovidae",
          "vs_HLhydIne1", "vs_HLmunMun1", "Cervidae",
          "vs_HLtraKan1", "mouseDeerOtherIsKept",
          "vs_HLmurAurFea1", "outerVespert",
          "vs_HLmyoLuc1", "Nearctic",
          "vs_myoDav1", "Myotis",
          "vs_HLpipPip1", "vs_HLlasBor1", "vs_HLnycHum2", "Vespertilioninae",
          "vs_HLmacSob1", "FoxLongTounge",
          "vs_HLeidHel2", "outerPeropodidae",
          "vs_HLeonSpe1", "Roussetinae"
         )')
}
{
args = c('r=Categorical4CategoryUnprunedTree', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 'v=T', 't=ER', 'n=ZoonomiaTip',
         'c=c(
              "C-Invertebrate-eater", "C-Endotherm-Carnivore", "C-Herpetivore", "C-Piscivore", "C-Nonspecific-Vertebrate-eater", "C-Scavenger", 
              "O-For Examination", "O-Scavenger", 
              "H-Frugivore", "H-Nectarivore", "H-Granivore", "H-Nonspecific-Herbivore", 
              "C-Terrestrial-vertebrates-eater", "C-All-vertebrate-eater", "C-All-Animals-Eater", 
              "H-High-sugar-plants-Eater", "H-Low-sugar-plants-Eater", "H-All-plants-Eater", 
              "O-Generalist", 
              "C-InsVertivore-Mixed", "C-InsVertivore-Piscivore", "C-InsVertivore-Insectivore","C-InsVertivore-Carnivore",
              "Insectivore", "Herpetivore", "Piscivore", "Vertivore", "InsVertivore", "Omnivore", "Frugivore", "Nectarivore", "Glucivore", "Herbivore", "Generalist"
            )', 
         'u=list(
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
          )')

args = c('r=CategoricalInsVertivoreTreeHeavyPrune', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.05',
         'c=c(
              "C-Invertebrate-eater", "C-Endotherm-Carnivore", "C-Herpetivore", "C-Piscivore", "C-Nonspecific-Vertebrate-eater", "C-Scavenger", 
              "O-For Examination", "O-Scavenger", 
              "H-Frugivore", "H-Nectarivore", "H-Granivore", "H-Nonspecific-Herbivore", 
              "C-Terrestrial-vertebrates-eater", "C-All-vertebrate-eater", "C-All-Animals-Eater", 
              "H-High-sugar-plants-Eater", "H-Low-sugar-plants-Eater", "H-All-plants-Eater", 
              "O-Generalist", 
              "C-InsVertivore-Mixed", "C-InsVertivore-Piscivore", "C-InsVertivore-Insectivore","C-InsVertivore-Carnivore",
              "Insectivore", "Herpetivore", "Piscivore", "Vertivore", "InsVertivore", "Omnivore", "Frugivore", "Nectarivore", "Glucivore", "Herbivore", "Generalist"
            )', 
         'u=list(
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
          )', 'y=c("vs_HLpanLeo1")')

args = c('r=CategoricalSlimMainInsVertivoreTree', 'm=data/CategoricalInsVertivoreMaintrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01',
         'c=c(
              "C-Invertebrate-eater", "C-Endotherm-Carnivore", "C-Herpetivore", "C-Piscivore", "C-Nonspecific-Vertebrate-eater", "C-Scavenger", 
              "O-For Examination", "O-Scavenger", 
              "H-Frugivore", "H-Nectarivore", "H-Granivore", "H-Nonspecific-Herbivore", 
              "C-Terrestrial-vertebrates-eater", "C-All-vertebrate-eater", "C-All-Animals-Eater", 
              "H-High-sugar-plants-Eater", "H-Low-sugar-plants-Eater", "H-All-plants-Eater", 
              "O-Generalist", 
              "C-InsVertivore-Mixed", "C-InsVertivore-Piscivore", "C-InsVertivore-Insectivore","C-InsVertivore-Carnivore",
              "Insectivore", "Herpetivore", "Piscivore", "Vertivore", "InsVertivore", "Omnivore", "Frugivore", "Nectarivore", "Glucivore", "Herbivore", "Generalist"
            )', 
         'u=list(
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
          )', 
         'y=c(
          "vs_HLornAna3", "vs_HLtacAcu1", "vs_HLdidVir1", "vs_HLgymLea1", "vs_HLpseCup1", 
          "vs_HLmyrTri1", "vs_HLchoDid1", "vs_HLchoHof3", "vs_HLproCap3", 
          "vs_HLpanLeo1", "LionClade", "vs_HLpanOnc1", "vs_HLaciJub2", "CheetahClade", "vs_HLursThi1", 
          "vs_ursMar1", "vs_HLursArc1", "vs_lepWed1", "vs_HLailMel2", 
          "vs_HLmirAng2", "vs_HLphoVit1", "HarborSealGreySealClade", "vs_HLeriBar1", 
          "vs_HLodoRos1", "vs_HLcalUrs1", "vs_HLzalCal1", 
          "vs_HLmelCap1", "vs_HLgulGul1", "vs_HLneoVis1", 
          "vs_HLpteBra1", "vs_HLlutLut1", "vs_enhLutKen1", 
          "vs_HLlycPic2", "vs_HLgloMel1", "vs_HLpepEle1", 
          "vs_HLturAdu1", "DolphinClade", "vs_orcOrc1", "vs_HLescRob1", "vs_HLlniGeo1", "amazonRiverDolphinFromYeast",
          "vs_HLbalEde1", "vs_HLmegNov1", "vs_HLcynGun1", 
          "vs_HLmerUng1", "vs_HLeulMon1", "vs_HLeulFul1", 
          "vs_eulMac1", "vs_HLeulFla1", "vs_ponAbe3", 
          "vs_panTro6", "vs_HLrhiRox2", "LangurClade", "vs_HLallNig1", 
          "vs_HLeryPat1", "vs_chlSab2", "geunonClade", "vs_HLtheGel1", 
          "vs_HLpapAnu5", "vs_HLmanSph1", "DrillMandrillClade", "vs_cerAty1", 
          "vs_HLtraJav1", "vs_mm10", "vs_HLmarFla1", "marmotClade"
         )')
}
args = c('r=CategoricalInsVertivoreTreePruningTest', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01',
         'c=c(
              "C-Invertebrate-eater", "C-Endotherm-Carnivore", "C-Herpetivore", "C-Piscivore", "C-Nonspecific-Vertebrate-eater", "C-Scavenger", 
              "O-For Examination", "O-Scavenger", 
              "H-Frugivore", "H-Nectarivore", "H-Granivore", "H-Nonspecific-Herbivore", 
              "C-Terrestrial-vertebrates-eater", "C-All-vertebrate-eater", "C-All-Animals-Eater", 
              "H-High-sugar-plants-Eater", "H-Low-sugar-plants-Eater", "H-All-plants-Eater", 
              "O-Generalist", 
              "C-InsVertivore-Mixed", "C-InsVertivore-Piscivore", "C-InsVertivore-Insectivore","C-InsVertivore-Carnivore",
              "Insectivore", "Herpetivore", "Piscivore", "Vertivore", "InsVertivore", "Omnivore", "Frugivore", "Nectarivore", "Glucivore", "Herbivore", "Generalist"
            )', 
         'u=list(
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
          )', 
         'y=c(
          "vs_HLornAna3", "vs_HLtacAcu1",  "PreserveMonotremeBranches"
          "vs_HLdidVir1", "vs_HLgymLea1", "vs_HLpseCup1", "MarsupialTransitionPreservation"
          "vs_HLmyrTri1", "vs_HLchoDid1", "vs_HLchoHof3", "vs_HLproCap3", "AfrotheriaPreserveTransitions"
          "vs_HLpanLeo1", "LionClade", "vs_HLpanOnc1", "vs_HLaciJub2", "CheetahClade", "PreserveBigCats" 
          "vs_HLursThi1", "vs_ursMar1", "vs_HLursArc1", "vs_HLailMel2", "UrsaPreserveTransition",
          "vs_lepWed1", "vs_HLmirAng2", "vs_HLphoVit1", "vs_HLeriBar1", "SealPreserveTransitions", 
          "vs_HLodoRos1", "vs_HLcalUrs1", "vs_HLzalCal1", "SealSeaLionPreserveTransitions"
          "vs_HLmelCap1", "vs_HLgulGul1", "vs_HLneoVis1", "MustelidPreserveTransitions"
          "vs_HLpteBra1", "vs_HLlutLut1", "vs_enhLutKen1", "OtterPreserveTransitions"
          "vs_HLlycPic2", "CanidPreserveTransision"
          "vs_HLgloMel1", "vs_HLpepEle1", "vs_HLturAdu1", "DolphinClade", "vs_orcOrc1", "vs_HLescRob1", "vs_HLlniGeo1", "amazonRiverDolphinFromYeast","vs_HLbalEde1", "vs_HLmegNov1", "vs_HLcynGun1", "CetaceaPreserveTransitions"
          "vs_HLmerUng1", "BankVoleTransition"
          "vs_HLeulMon1", "vs_HLeulFul1", "LemurTransition"
          "vs_eulMac1", "vs_HLeulFla1", "vs_ponAbe3", "PrimateTransitions"
          "vs_panTro6", "vs_HLrhiRox2", "LangurClade", "vs_HLallNig1", "PrimateTransitionsContinued"
          "vs_HLeryPat1", "vs_chlSab2", "geunonClade", "vs_HLtheGel1", "PrimateTransitionsContinuedAgain"
          "vs_HLpapAnu5", "vs_HLmanSph1", "DrillMandrillClade", "vs_cerAty1", "Drilltransitions"
          "vs_HLtraJav1", "MouseDeerUsedInOtherAnalysis
          "vs_mm10", "humans"
          "vs_HLmarFla1", "marmotClade" "DoormouseTransition"
         )',
         'p=c(
          "vs_HLellTal1", "vs_HLellLut1", "vs_HLarvAmp1", "voleClade",
          "vs_HLmusSpi1", "vs_HLmusCar1", "vs_HLmasCou1", "vs_HLmusPah1", "mouseClade",
          "vs_HLhysCri1", "vs_HLthrSwi1", "vs_HLpetTyp1", "vs_hetGla2", "vs_chiLan1", "vs_HLdinBra1", "vs_HLcteSoc1", "vs_octDeg1", "vs_HLcoePre1", "vs_HLdasPun1", "vs_HLdolPat1", "gundiGuineaPigClade",
          "vs_HLoryGaz1", "vs_HLbeaHun1", "vs_HLkobLecLec1", "vs_HLkobLecLec1", "vs_HLmadKir1", "vs_HLneoPyg1", "vs_HLphiMax1", "vs_HLoreOre1", "vs_HLneoMos1", "vs_HLaepMel1", "vs_HLtraImb1", "Bovidae",
          "vs_HLhydIne1", "vs_HLmunMun1", "Cervidae",
          "vs_HLtraKan1", "mouseDeerOtherIsKept",
          "vs_HLmurAurFea1", "outerVespert",
          "vs_HLmyoLuc1", "Nearctic",
          "vs_myoDav1", "Myotis",
          "vs_HLpipPip1", "vs_HLlasBor1", "vs_HLnycHum2", "Vespertilioninae",
          "vs_HLmacSob1", "FoxLongTounge",
          "vs_HLeidHel2", "outerPeropodidae",
          "vs_HLeonSpe1", "Roussetinae"
         )')

args = c('r=CategoricalInsVertivoreTreeLiamInference', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01',
         'c=c(
              "C-Invertebrate-eater", "C-Endotherm-Carnivore", "C-Herpetivore", "C-Piscivore", "C-Nonspecific-Vertebrate-eater", "C-Scavenger", 
              "O-For Examination", "O-Scavenger", 
              "H-Frugivore", "H-Nectarivore", "H-Granivore", "H-Nonspecific-Herbivore", 
              "C-Terrestrial-vertebrates-eater", "C-All-vertebrate-eater", "C-All-Animals-Eater", 
              "H-High-sugar-plants-Eater", "H-Low-sugar-plants-Eater", "H-All-plants-Eater", 
              "O-Generalist", 
              "C-InsVertivore-Mixed", "C-InsVertivore-Piscivore", "C-InsVertivore-Insectivore","C-InsVertivore-Carnivore",
              "Insectivore", "Herpetivore", "Piscivore", "Vertivore", "InsVertivore", "Omnivore", "Frugivore", "Nectarivore", "Glucivore", "Herbivore", "Generalist"
            )', 
         'u=list(
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
          )', 
         'y=c(
          "vs_HLornAna3", "vs_HLtacAcu1",  "PreserveMonotremeBranches",
          "vs_HLdidVir1", "vs_HLgymLea1", "vs_HLpseCup1", "MarsupialTransitionPreservation",
          "vs_HLmyrTri1", "vs_HLchoDid1", "vs_HLchoHof3", "vs_HLproCap3", "AfrotheriaPreserveTransitions",
          "vs_HLpanLeo1", "LionClade", "vs_HLpanOnc1", "vs_HLaciJub2", "CheetahClade", "PreserveBigCats", 
          "vs_HLursThi1", "vs_ursMar1", "vs_HLursArc1", "vs_HLailMel2", "UrsaPreserveTransition",
          "vs_lepWed1", "vs_HLmirAng2", "vs_HLphoVit1", "vs_HLeriBar1", "SealPreserveTransitions", 
          "vs_HLodoRos1", "vs_HLcalUrs1", "vs_HLzalCal1", "SealSeaLionPreserveTransitions",
          "vs_HLmelCap1", "vs_HLgulGul1", "vs_HLneoVis1", "MustelidPreserveTransitions",
          "vs_HLpteBra1", "vs_HLlutLut1", "vs_enhLutKen1", "OtterPreserveTransitions",
          "vs_HLlycPic2", "CanidPreserveTransision",
          "vs_HLgloMel1", "vs_HLpepEle1", "vs_HLturAdu1", "DolphinClade", "vs_orcOrc1", "vs_HLescRob1", "vs_HLlniGeo1", "amazonRiverDolphinFromYeast","vs_HLbalEde1", "vs_HLmegNov1", "vs_HLcynGun1", "CetaceaPreserveTransitions",
          "vs_HLmerUng1", "BankVoleTransition",
          "vs_HLeulMon1", "vs_HLeulFul1", "LemurTransition",
          "vs_eulMac1", "vs_HLeulFla1", "vs_ponAbe3", "PrimateTransitions",
          "vs_panTro6", "vs_HLrhiRox2", "LangurClade", "vs_HLallNig1", "PrimateTransitionsContinued",
          "vs_HLeryPat1", "vs_chlSab2", "geunonClade", "vs_HLtheGel1", "PrimateTransitionsContinuedAgain",
          "vs_HLpapAnu5", "vs_HLmanSph1", "DrillMandrillClade", "vs_cerAty1", "Drilltransitions",
          "vs_HLtraJav1", "MouseDeerUsedInOtherAnalysis",
          "vs_mm10", "humans",
          "vs_HLmarFla1", "marmotClade", "DoormouseTransition"
         )',
         'p=c(
          "vs_HLellTal1", "vs_HLellLut1", "vs_HLarvAmp1", "voleClade",
          "vs_HLmusSpi1", "vs_HLmusCar1", "vs_HLmasCou1", "vs_HLmusPah1", "mouseClade",
          "vs_HLhysCri1", "vs_HLthrSwi1", "vs_HLpetTyp1", "vs_hetGla2", "vs_chiLan1", "vs_HLdinBra1", "vs_HLcteSoc1", "vs_octDeg1", "vs_HLcoePre1", "vs_HLdasPun1", "vs_HLdolPat1", "gundiGuineaPigClade",
          "vs_HLoryGaz1", "vs_HLbeaHun1", "vs_HLkobLecLec1", "vs_HLkobLecLec1", "vs_HLmadKir1", "vs_HLneoPyg1", "vs_HLphiMax1", "vs_HLoreOre1", "vs_HLneoMos1", "vs_HLaepMel1", "vs_HLtraImb1", "Bovidae",
          "vs_HLhydIne1", "vs_HLmunMun1", "Cervidae",
          "vs_HLtraKan1", "mouseDeerOtherIsKept",
          "vs_HLmurAurFea1", "outerVespert",
          "vs_HLmyoLuc1", "Nearctic",
          "vs_myoDav1", "Myotis",
          "vs_HLpipPip1", "vs_HLlasBor1", "vs_HLnycHum2", "Vespertilioninae",
          "vs_HLmacSob1", "FoxLongTounge",
          "vs_HLeidHel2", "outerPeropodidae",
          "vs_HLeonSpe1", "Roussetinae"
         )')



# --- Standard start-up code ---
if(clusterRun){args = commandArgs(trailingOnly = TRUE)}
{  # Bracket used for collapsing purposes
  #File Prefix
  if(!is.na(cmdArgImport('r'))){
    filePrefix = cmdArgImport('r')
  }else{
    stop("THIS IS AN ISSUE MESSAGE; SPECIFY FILE PREFIX")
  }
  
  #  Output Directory 
  if(!dir.exists("Output")){                                      #Make output directory if it does not exist
    dir.create("Output")
  }
  outputFolderNameNoSlash = paste("Output/",filePrefix, sep = "") #Set the prefix sub directory
  if(!dir.exists(outputFolderNameNoSlash)){                       #create that directory if it does not exist
    dir.create(outputFolderNameNoSlash)
  }
  outputFolderName = paste("Output/",filePrefix,"/", sep = "")
  
  #  Force update argument
  forceUpdate = FALSE
  if(!is.na(cmdArgImport('v'))){                                 #Import if update being forced with argument 
    forceUpdate = cmdArgImport('v')
    forceUpdate = as.logical(forceUpdate)
  }else{
    message("Force update not specified, not forcing update")
  }
}

# --- Argument Imports ---
{ # Bracket used for collapsing purposes
# Defaults
mainTreesLocation = "/share/ceph/wym219group/shared/projects/MammalDiet/Zoonomia/RemadeTreesAllZoonomiaSpecies.rds"
spreadSheetLocation = "Data/manualAnnotationsSheet.csv"
annotColumn = NULL
categoryList = NULL
useScreen = F
screenColumn = NULL
modelType = "ER"
ancestralTrait = NULL
substitutions = NULL
nameColumn = "tipName"
usingPruning = F
usingAutoPruning = F
manualPruningProtections = NULL
pruningPrefrenceColumn = NA
pruningProtection = F
manualPruningSpecies = NULL

  #MainTrees Location
  if(!is.na(cmdArgImport('m'))){
    mainTreesLocation = cmdArgImport('m')
  }else{
    message("No maintrees arg, using default")
  }
  #read in the tree based on filetype extension
  if(file_ext(mainTreesLocation) == "rds"){
    if(!exists("mainTrees")){mainTrees = readRDS(mainTreesLocation)}
  }else{
    if(!exists("mainTrees")){mainTrees = readTrees(mainTreesLocation)} 
  }

  #spreadsheet File
  if(!is.na(cmdArgImport('d'))){
    spreadSheetLocation = cmdArgImport('d')
  }else{
    message("Using Data/manualAnnotationsSheet.csv spreadsheet")
  }
  
  #Annots Column
  if(!is.na(cmdArgImport('a'))){
    annotColumn = cmdArgImport('a')
  }else{
    stop("THIS IS AN ISSUE MESSAGE; SPECIFY ANNOTATION COLLUMN")
  }
  
  #Category list 
  if(!is.null(cmdArgImport('c'))){
    categoryList = cmdArgImport('c')
  }else{
    stop("THIS IS AN ISSUE MESSAGE; SPECIFY CATEGORIES")
  }
  
  #Screen Column
  if(!is.na(cmdArgImport('s'))){
    useScreen = T
    screenColumn = cmdArgImport('s')
  }else{
    message("No screen column used.")
  }
  
  #Model Type
  if(!is.na(cmdArgImport('t'))){
    modelType = cmdArgImport('t')
  }else{
    message("No model specified, using ER.")
  }
  
  #Ancestral Trait
  if(!is.na(cmdArgImport('g'))){
    ancestralTrait = cmdArgImport('g')
  }else{
    message("No ancestral trait specified, using NULL")
  }
  
  #Substitution list 
  if(!is.null(cmdArgImport('u'))){
    substitutions = cmdArgImport('u')
  }else{
    message("No substitutions provided")
  }
  #Merge list 
  if(!is.null(cmdArgImport('o'))){
    mergeOnlys = cmdArgImport('o')
  }else{
    message("No merges provided")
  }
  #Name Column
  if(!is.na(cmdArgImport('n'))){
    nameColumn = cmdArgImport('n')
  }else{
    message("Name Column not specified, using 'tipName'.")
  }

  #Pruning cutoff
  if(!is.na(cmdArgImport('z'))){
    usingPruning = T
    usingAutoPruning = T
    pruningCutoff = cmdArgImport('z')
  }else{
    message("Pruning Cutoff not specified, not pruning tree.")
  }
  
  #PruningPrefrenceColumn 
  if(!is.na(cmdArgImport('x'))){
    pruningPrefrenceColumn = cmdArgImport('x')
  }else{
    if(usingPruning){message("No pruning prefrence column specified")}
  }

  #ManualPruningProtections
  if(!all(is.na(cmdArgImport('y')))){
    manualPruningProtections = cmdArgImport('y')
  }else{
    if(usingPruning){message("No manually protected species specified")}
  }

  #ManualPruningSpecies
  if(!all(is.na(cmdArgImport('p')))){
    manualPruningSpecies = cmdArgImport('p')
    usingPruning = T
  }else{
    if(usingPruning){message("No manually pruned species specified")}
  }


}


#                   ------- Code Body --------  

manualAnnots = read.csv(spreadSheetLocation)                      #load the manual annotations file holding the phenotype data
manualAnnots[[annotColumn]] = trimws(manualAnnots[[annotColumn]])               #trim away whitespace to allow for better matching 

# - Merge hyrbid of either substituted phenotypes or merge-only phenotypes - 
if(!is.null(substitutions) & !all(is.na(substitutions))){                                                    #Consider species with multiple combined categories as the merged category
  for( i in 1:length(substitutions)){                                           #Eg if [X] is replaced with [Y], [X/Y] becomes [Y]
    substitutePhenotypes = substitutions[[i]]
    message(paste("Combining", substitutePhenotypes[1], "/", substitutePhenotypes[2]))
    entriesWithPhen1 = grep(substitutePhenotypes[1], manualAnnots[[annotColumn]])
    entriesWithPhen2 = grep(substitutePhenotypes[2], manualAnnots[[annotColumn]])
    combineEntries = which(entriesWithPhen1 %in% entriesWithPhen2)
    combineIndexes = entriesWithPhen1[combineEntries]
    manualAnnots[[annotColumn]][combineIndexes] = substitutePhenotypes[2]
  }
}

if(!is.null(mergeOnlys) & ! all(is.na(mergeOnlys))){                                                    #Consider species with multiple combined categories as the merged category
  for( i in 1:length(mergeOnlys)){                                           #Eg if [X] is replaced with [Y], [X/Y] becomes [Y]
    substitutePhenotypes = mergeOnlys[[i]]
    message(paste("Merging Hybrids of", substitutePhenotypes[1], "/", substitutePhenotypes[2], "to", substitutePhenotypes[2]))
    entriesWithPhen1 = grep(substitutePhenotypes[1], manualAnnots[[annotColumn]])
    entriesWithPhen2 = grep(substitutePhenotypes[2], manualAnnots[[annotColumn]])
    combineEntries = which(entriesWithPhen1 %in% entriesWithPhen2)
    combineIndexes = entriesWithPhen1[combineEntries]
    manualAnnots[[annotColumn]][combineIndexes] = substitutePhenotypes[2]
  }
}


# - Species Filter - 
speciesFilterFilename = paste(outputFolderName, filePrefix, "SpeciesFilter.rds",sep="") #set a filename for the species filter based on the prefix 

if(!file.exists(speciesFilterFilename) | forceUpdate){                          #if no filter exists or update is forced, make a filter 
  # --- subset the manual annots to only those with data in the categories used, and optionally by the screen column
  relevantSpecies = manualAnnots[manualAnnots[[annotColumn]] %in% categoryList,]#remove all species which are not part of the specified categories
  if(useScreen){                                                                #if using a screening collumn 
    relevantSpecies = relevantSpecies[ relevantSpecies[[screenColumn]] %in% 1, ]  #remove all species not positive for that collumn 
  }
  relevantSpecies = relevantSpecies[!relevantSpecies[[nameColumn]] %in% "", ]          #remove any species without an FA name (not on the master tree)
  speciesFilter = relevantSpecies[[nameColumn]]                                       #make a list of the master tree tip labels of the included species
  if(usingPruning){
    if(usingAutoPruning){
      source("Src/Reu/autoPruner.R")
      pruningProtectionSpecies = NULL
      if(!is.na(pruningPrefrenceColumn)){
        if(all(is.logical(manualAnnots[[pruningPrefrenceColumn]]))){
          pruningProtection = T
        }else{ 
          pruningProtection = F
        }
        
        pruningProtectionRows = manualAnnots[which(as.logical(manualAnnots[[pruningPrefrenceColumn]])),]
        pruningProtectionSpecies = pruningProtectionRows[[nameColumn]]
      }
      allProtectedSpecies = append(pruningProtectionSpecies, manualPruningProtections)
      if(all(is.null(allProtectedSpecies))){allProtectedSpecies = ""}
      
      workingTree = mainTrees$masterTree
      workingTree = drop.tip(workingTree, which(!workingTree$tip.label %in% speciesFilter))
      
      fewGeneSpecies = dropFewGeneSpecies(mainTrees, workingTree, nameConversionColumn = nameColumn, nameConversionData = spreadSheetLocation)
      if(length(which(fewGeneSpecies %in% allProtectedSpecies)) > 0){
        fewGeneSpecies = fewGeneSpecies[-which(fewGeneSpecies %in% allProtectedSpecies)]
      }
      workingTree = drop.tip(workingTree, fewGeneSpecies)
      names(fewGeneSpecies)[1:length(fewGeneSpecies)] = "fewGenes"
      
      pruningFilename = paste(outputFolderName, filePrefix, "PruningTree.pdf", sep="")
      pdf(pruningFilename, width = 16, height = length(workingTree$tip.label)/8)
      prunedTree = autopruner(workingTree, dropValue = pruningCutoff, tipsToKeep = allProtectedSpecies, nameConversionColumn = nameColumn, nameConversionData = spreadSheetLocation, preDroppedTips = fewGeneSpecies)
      if(!pruningProtection){
        prunedTree = autopruner(prunedTree, dropValue = pruningCutoff, tipsToKeep = manualPruningProtections, nameConversionColumn = nameColumn, nameConversionData = spreadSheetLocation, preDroppedTips = droppedTips, originalTree = workingTree)
      }
      dev.off()
    }
    if(!all(is.null(manualPruningSpecies))){
      prunedTree = drop.tip(prunedTree, manualPruningSpecies)
      names(manualPruningSpecies)[1:length(manualPruningSpecies)] = "manualDrop"
      droppedTips = append(droppedTips, manualPruningSpecies)
    }
    
    
    prunedSpecies = speciesFilter[!speciesFilter %in% prunedTree$tip.label]
    speciesFilter = speciesFilter[-which(speciesFilter %in% prunedSpecies)]
    
    prunedSpeciesFilename = paste(outputFolderName, filePrefix, "prunedSpecies.rds",sep="")
    saveRDS(droppedTips, prunedSpeciesFilename)
    prunedSpeciesTextFilename = file(paste(outputFolderName, filePrefix, "prunedSpecies.txt",sep=""))
    writeLines(print(droppedTips),prunedSpeciesTextFilename)
    close(prunedSpeciesTextFilename)
  }
  
  
  
  saveRDS(speciesFilter, file = speciesFilterFilename)                          #save that as the species filter
  
  irrelevantSpecies = manualAnnots[! manualAnnots[[nameColumn]] %in% speciesFilter,]
}else{ #if not, use the existing one 
  relevantSpecieslist = readRDS(speciesFilterFilename)                          #if not, use the existing list 
  speciesFilter = relevantSpecieslist                                           #make the speciesFilter object for later 
  relevantSpecies = manualAnnots[ manualAnnots[[nameColumn]] %in% relevantSpecieslist,] #and select the manual annotations entries in that list (useful if the list is more restrictive than it would be by default) 
  irrelevantSpecies = manualAnnots[! manualAnnots[[nameColumn]] %in% relevantSpecieslist,]
}

# - Phenotype Vector - 
relevantSpecies = relevantSpecies[relevantSpecies[[nameColumn]] %in% speciesFilter, ] 
speciesNames = relevantSpecies[[nameColumn]]                                         #Exract the tip name of each species
speciesCategories = relevantSpecies[[annotColumn]]                              #extract the category of each species (in same order)

phenotypeVector = speciesCategories                                             #combine those into⌄
names(phenotypeVector) = speciesNames                                           #the format the functions expect
if(!is.null(substitutions) & !all(is.na(substitutions))){
  for( i in 1:length(substitutions)){
    substitutePhenotypes = substitutions[[i]]
    message(paste("replacing", substitutePhenotypes[1], "with", substitutePhenotypes[2]))
    phenotypeVector = gsub(substitutePhenotypes[1], substitutePhenotypes[2], phenotypeVector)
  }
}

phenotypeVectorFilename = paste(outputFolderName, filePrefix, "CategoricalPhenotypeVector.rds",sep="") #make a filename based on the prefix
saveRDS(phenotypeVector, file = phenotypeVectorFilename)                        #save the phenotype vector

# - Make common name versions of objects (used in visualization) - 
commonMainTrees = mainTrees
commonMainTrees$masterTree = ZoonomTreeNameToCommon(commonMainTrees$masterTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)
commonPhenotypeVector = phenotypeVector
names(commonPhenotypeVector) = ZonomNameConvertVectorCommon(names(commonPhenotypeVector), annotationLocation = spreadSheetLocation, tipColumn = nameColumn)
commonSpeciesFilter = ZonomNameConvertVectorCommon(speciesFilter, annotationLocation = spreadSheetLocation, tipColumn = nameColumn)

# - Categorical Tree - 
treeImageFilename = paste(outputFolderName, filePrefix, "CategoricalTree.pdf", sep="") #make a filename based on the prefix
pdf(treeImageFilename, height = length(phenotypeVector)/18, width = 10)                     #make a pdf to store the plot, sized based on tree size
  commonCategoricalTree = char2TreeCategorical(commonPhenotypeVector, commonMainTrees, commonSpeciesFilter, model = modelType, anctrait = ancestralTrait, plot = T)

  categoricalTree = char2TreeCategorical(phenotypeVector, mainTrees, speciesFilter, model = modelType, anctrait = ancestralTrait, plot = T) #use the phenotype vector to make a tree
dev.off()                                                                       #save the plot to the pdf
 
scientificCategoricalTree = ZoonomTreeNameToCommon(commonCategoricalTree, manualAnnotLocation = spreadSheetLocation, tipCol = "CommonName", scientific = T, scientificCol = "Scientific_Binomial", plot = F)

categoricalTreeFilename = paste(outputFolderName, filePrefix, "CategoricalTree.rds", sep="") #make a filename based on the prefix
saveRDS(categoricalTree, categoricalTreeFilename)                               #save the tree
categoricalCommonTreeFilename = paste(outputFolderName, filePrefix, "CategoricalCommonTree.rds", sep="") #make a filename based on the prefix
saveRDS(commonCategoricalTree, categoricalCommonTreeFilename)
scientificCategoricalTreeFilename = paste(outputFolderName, filePrefix, "CategoricalScientificTree.rds", sep="") #make a filename based on the prefix
saveRDS(scientificCategoricalTree, scientificCategoricalTreeFilename)



# - Paths - 
#pathsFilename = paste(outputFolderName, filePrefix, "CategoricalPathsFile.rds", sep= "") #make a filename based on the prefix
#paths = char2PathsCategorical(phenotypeVector, mainTrees, speciesFilter, model = modelType, anctrait = ancestralTrait) #make a path based on the phenotype vector
#saveRDS(paths, file = pathsFilename)                                            #save the path 






# -- Convert Tree to Binary (Manual only) --
convertToBinary = T
convertToBinary = F
convertToBinaryStandalone = T
convertToBinaryStandalone = F
foreground = "Insectivore"

if(convertToBinaryStandalone){
  categoricalTreeFilename = paste(outputFolderName, filePrefix, "CategoricalTree.rds", sep="") #make a filename based on the prefix
  categoricalTree = readRDS(categoricalTreeFilename)
  phenotypeVectorFilename = paste(outputFolderName, filePrefix, "CategoricalPhenotypeVector.rds",sep="") #make a filename based on the prefix
  speciesFilterFilename = paste(outputFolderName, filePrefix, "SpeciesFilter.rds",sep="") #set a filename for the species filter based on the prefix 
  relevantSpecieslist = readRDS(speciesFilterFilename)                          #if not, use the existing list 
  speciesFilter = relevantSpecieslist                                           #make the speciesFilter object for later 
}

if(convertToBinary){
  
  binaryTree = categoricalTree
  phenotypeVector = readRDS(phenotypeVectorFilename)                            #load in the phenotype vector 
  categories = map_to_state_space(phenotypeVector) 
  categoryNames = categories$name2index                                         #store the length-phenotype connection
  foregroundInt = categoryNames[which(names(categoryNames) == foreground)]
  binaryTree$edge.length[-(which(binaryTree$edge.length == foregroundInt))] = 0
  binaryTree$edge.length[(which(binaryTree$edge.length == foregroundInt))] = 1
  
  binaryTreeImageFilename = paste(outputFolderName, filePrefix, "BinaryTree.pdf", sep="") #make a filename based on the prefix
  pdf(binaryTreeImageFilename, height = length(phenotypeVector)/7)                     #make a pdf to store the plot, sized based on tree size
  source("Src/Reu/plotBinaryTree.R")
  plotBinaryTree(mainTrees, binaryTree, mainTitle = paste(filePrefix, "Binary", "Foreground", "Tree"), tipColumn = nameColumn)
  plotBinaryTree(mainTrees, binaryTree, convertNames = F, mainTitle = paste(filePrefix, "Binary", "Foreground", "Tree"), tipColumn = nameColumn)
  plotTree(binaryTree)
  dev.off()                                                                       #save the plot to the pdf
  
  binaryTreeFilename = paste(outputFolderName, filePrefix, "BinaryTree.rds", sep="") #make a filename based on the prefix
  saveRDS(binaryTree, binaryTreeFilename)                               #save the tree
  
  # - Paths - 
  binaryPathsFilename = paste(outputFolderName, filePrefix, "PathsFile.rds", sep= "") #make a filename based on the prefix
  binaryPaths = tree2Paths(binaryTree, mainTrees, binarize = T, speciesFilter) #make a path based on the phenotype vector
  saveRDS(binaryPaths, file = binaryPathsFilename) 
}
