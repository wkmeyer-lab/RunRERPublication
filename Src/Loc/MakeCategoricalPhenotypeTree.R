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

#When running on the command line, note that the arguments should be formatted as:
#r=Prefix v=forceUpdateValue etc
#Rather than being in quotes and separated by commas. 
#Terminal uses spaces to split indexes in a vector, where r uses indepent quotation and commas.


# -- Command arguments list
# r = filePrefix                                         This is a prefix used to organize and separate files by analysis run. Always required. 
# v = <T or F>                                           This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# m = mainTreeFilename.txt or .rds                       This sets the location of the maintrees file
# d = spreadSheetFilename.csv                            This sets the spreadsheet to read the data from 
# a = "annotCollumn"                                     This is the column in the manual annotations spreadsheet to use
# n = "nameColumn"                                       This sets the column with the tip names as they appear in the maintrees file. 
# c = <c("nameOfCategory1,"nameOfCategory2")>            This is the list of category names 
# u = list(c("replace1", "with1"),c("replace2, with2"))
# o = list(c("phenotype1", "intophen1"), c("2", "i2"))   This causes combination phenotypes to be merged into the second phenotype, but does not replace standalone phenotypes
# s = "screenCollumn"                                    This is a collumn which must have a value of 1 for the species to be included. 
# t = <ER or SYM or ARD>                                 This sets the model type used to estimate ancestral branches 
# g = "ancestralTrait"                                   This can be used to set all non-terminal branches to this category. Use be one of the categories in the list. 
# z = <minimum branch length>                            This sets the minimum branch length for terminal branches in the master tree. Branches shorter than this will be removed. 
# x = "pruningPrefrenceColumn"                           This sets a column, where if the value is 1, the tip will be preferentially kept. If the value is TRUE, the tip will never be pruned.
# y = "c('unprunedtip1', 'unprunedtip2')"                This allows you to add a list of specific tips to not be dropped during pruning. Must use the tip name, not common name. 
# p = "c('prunedtip1', 'prunedtip2')"                    This allows you to manually specify additional branches to be pruned
# l = <T or F>                                           DO NOT USE UNLESS SPECIFIED -- This determines if the liam infrence nodes should be added 
# e = <T or F>                                           This sets if alternate species sets should be created. 
#----------------

{
  args = c('r=CladeBinaryBovidae', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'v=T', 't=ER', 'n=ZoonomiaTip', 'l=F', 'e=T',
           'c=c(0, 1)',
           'a=isBovidae','s=isInAnalysisWithFullFamilies') 
  args = c('r=CladeBinaryCervidae', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'v=T', 't=ER', 'n=ZoonomiaTip', 'l=F', 'e=T',
           'c=c(0, 1)',
           'a=isCervidae','s=isInAnalysisWithFullFamilies') 
  args = c('r=CladeBinaryCricetidae', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'v=T', 't=ER', 'n=ZoonomiaTip', 'l=F', 'e=T',
           'c=c(0, 1)',
           'a=isCricetidae','s=isInAnalysisWithFullFamilies') 
  args = c('r=CladeBinaryHystricognathi', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'v=T', 't=ER', 'n=ZoonomiaTip', 'l=F', 'e=T',
           'c=c(0, 1)',
           'a=isHystricognathi','s=isInAnalysisWithFullFamilies') 
  args = c('r=CladeBinaryVespertilionidae', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'v=T', 't=ER', 'n=ZoonomiaTip', 'l=F', 'e=T',
           'c=c(0, 1)',
           'a=isVespertilionidae','s=isInAnalysisWithFullFamilies') 
  args = c('r=CladeBinaryPeropdidae', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'v=T', 't=ER', 'n=ZoonomiaTip', 'l=F', 'e=T',
           'c=c(0, 1)',
           'a=isPeropdidae','s=isInAnalysisWithFullFamilies') 
  args = c('r=CladeBinaryPeropdidae', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'v=T', 't=ER', 'n=ZoonomiaTip', 'l=F', 'e=T',
           'c=c(0, 1)',
           'a=isCervidae','s=isInAnalysisWithFullFamilies') 
}
{ 
  args = c('r=ComplexDietCentralAnalysisSimplifyAllSpecies', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=SimplifiedDietConvert', 'v=T', 't=ER', 'n=ZoonomiaTip', 'l=T',
           'c=c(
              "Vertivore", "Omnivore", "Herbivore", "Invertivore"
            )'
           )
  
  args = c('r=ComplexDietCentralAnalysisSimplify2', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=SimplifiedDietConvert', 'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01', 'l=T', 'e=T',
           'c=c(
              "Vertivore", "Omnivore", "Herbivore", "Invertivore"
            )', 
           'y=c(
          "vs_HLornAna3", "vs_HLtacAcu1",  "PreserveMonotremeBranches",
          "vs_HLdidVir1", "vs_HLgymLea1", "vs_HLpseCup1", "MarsupialTransitionPreservation",
          "vs_HLmyrTri1", "vs_HLchoDid1", "vs_HLchoHof3", "vs_HLproCap3", "AfrotheriaPreserveTransitions",
          "vs_HLursThi1", "vs_ursMar1", "vs_HLursArc1", "vs_HLailMel2", "UrsaPreserveTransition",
          "vs_lepWed1", "vs_HLmirAng2", "vs_HLphoVit1", "vs_HLeriBar1", "SealPreserveTransitions",
          "vs_HLodoRos1", "vs_HLcalUrs1", "vs_HLzalCal1", "SealSeaLionPreserveTransitions",
          "vs_HLmelCap1", "vs_HLgulGul1", "vs_HLneoVis1", "MustelidPreserveTransitions",
          "vs_HLlycPic2", "CanidPreserveTransision",
          "vs_HLgloMel1", "vs_HLpepEle1", "vs_HLturAdu1", "DolphinClade", "vs_orcOrc1", "vs_HLescRob1", "vs_HLlniGeo1", "amazonRiverDolphinFromYeast","vs_HLbalEde1", "vs_HLmegNov1", "vs_HLcynGun1", "CetaceaPreserveTransitions",
          "vs_HLmerUng1", "BankVoleTransition",
          "vs_HLeulMon1", "vs_HLeulFul1", "LemurTransition",
          "vs_panTro6", "vs_HLrhiRox2", "LangurClade", "vs_HLallNig1", "PrimateTransitionsContinued",
          "vs_HLeryPat1", "vs_chlSab2", "geunonClade", "vs_HLtheGel1", "PrimateTransitionsContinuedAgain",
          "vs_HLpapAnu5", "vs_HLmanSph1", "DrillMandrillClade", "vs_cerAty1", "Drilltransitions",
          "vs_HLmarFla1", "marmotClade", "DoormouseTransition",
          "vs_eulMac1", "vs_ponAbe3", "PrimateTransitions"
          
         )',
         'p=c(
          "vs_HLellTal1", "vs_HLellLut1", "vs_HLarvAmp1", "voleClade",
          "vs_HLhysCri1", "vs_HLthrSwi1", "vs_HLpetTyp1", "vs_hetGla2", "vs_chiLan1", "vs_HLdinBra1", "vs_HLcteSoc1", "vs_octDeg1", "vs_HLcoePre1", "vs_HLdasPun1", "vs_HLdolPat1", "gundiGuineaPigClade",
          "vs_HLoryGaz1", "vs_HLbeaHun1", "vs_HLkobLecLec1", "vs_HLkobLecLec1", "vs_HLmadKir1", "vs_HLneoPyg1", "vs_HLphiMax1", "vs_HLoreOre1", "vs_HLneoMos1", "vs_HLaepMel1", "vs_HLtraImb1", "Bovidae",
          "vs_HLhydIne1", "vs_HLmunMun1", "Cervidae",
          "vs_HLmurAurFea1", "outerVespert",
          "vs_HLmyoLuc1", "Nearctic",
          "vs_myoDav1", "Myotis",
          "vs_HLpipPip1", "vs_HLlasBor1", "vs_HLnycHum2", "Vespertilioninae",
          "vs_HLmacSob1", "FoxLongTounge",
          "vs_HLeidHel2", "outerPeropodidae",
          "vs_HLeonSpe1", "Roussetinae"
         )')
  
  args = c('r=ComplexDietCentralAnalysisSimplifyExclusivePredAllSpecies', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=SimplifiedDietConvertStrictPred', 'v=T', 't=ER', 'n=ZoonomiaTip', 'l=T', 'e=T',
           'c=c(
              "Vertivore", "Omnivore", "Herbivore", "Invertivore"
            )')
  args = c('r=ComplexDietCentralAnalysisSimplifyExclusivePred', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=SimplifiedDietConvertStrictPred', 'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01', 'l=T', 'e=T',
           'c=c(
              "Vertivore", "Omnivore", "Herbivore", "Invertivore"
            )', 
           'y=c(
          "vs_HLornAna3", "vs_HLtacAcu1",  "PreserveMonotremeBranches",
          "vs_HLdidVir1", "vs_HLgymLea1", "vs_HLpseCup1", "MarsupialTransitionPreservation",
          "vs_HLmyrTri1", "vs_HLchoDid1", "vs_HLchoHof3", "vs_HLproCap3", "AfrotheriaPreserveTransitions",
          "vs_HLursThi1", "vs_ursMar1", "vs_HLursArc1", "vs_HLailMel2", "UrsaPreserveTransition",
          "vs_lepWed1", "vs_HLmirAng2", "vs_HLphoVit1", "vs_HLeriBar1", "SealPreserveTransitions",
          "vs_HLodoRos1", "vs_HLcalUrs1", "vs_HLzalCal1", "SealSeaLionPreserveTransitions",
          "vs_HLmelCap1", "vs_HLgulGul1", "vs_HLneoVis1", "MustelidPreserveTransitions",
          "vs_HLlycPic2", "CanidPreserveTransision",
          "vs_HLgloMel1", "vs_HLpepEle1", "vs_HLturAdu1", "DolphinClade", "vs_orcOrc1", "vs_HLescRob1", "vs_HLlniGeo1", "amazonRiverDolphinFromYeast","vs_HLbalEde1", "vs_HLmegNov1", "vs_HLcynGun1", "CetaceaPreserveTransitions",
          "vs_HLmerUng1", "BankVoleTransition",
          "vs_HLeulMon1", "vs_HLeulFul1", "LemurTransition",
          "vs_panTro6", "vs_HLrhiRox2", "LangurClade", "vs_HLallNig1", "PrimateTransitionsContinued",
          "vs_HLeryPat1", "vs_chlSab2", "geunonClade", "vs_HLtheGel1", "PrimateTransitionsContinuedAgain",
          "vs_HLpapAnu5", "vs_HLmanSph1", "DrillMandrillClade", "vs_cerAty1", "Drilltransitions",
          "vs_HLmarFla1", "marmotClade", "DoormouseTransition",
          "vs_eulMac1", "vs_ponAbe3", "PrimateTransitions"
          
         )',
         'p=c(
          "vs_HLellTal1", "vs_HLellLut1", "vs_HLarvAmp1", "voleClade",
          "vs_HLhysCri1", "vs_HLthrSwi1", "vs_HLpetTyp1", "vs_hetGla2", "vs_chiLan1", "vs_HLdinBra1", "vs_HLcteSoc1", "vs_octDeg1", "vs_HLcoePre1", "vs_HLdasPun1", "vs_HLdolPat1", "gundiGuineaPigClade",
          "vs_HLoryGaz1", "vs_HLbeaHun1", "vs_HLkobLecLec1", "vs_HLkobLecLec1", "vs_HLmadKir1", "vs_HLneoPyg1", "vs_HLphiMax1", "vs_HLoreOre1", "vs_HLneoMos1", "vs_HLaepMel1", "vs_HLtraImb1", "Bovidae",
          "vs_HLhydIne1", "vs_HLmunMun1", "Cervidae",
          "vs_HLmurAurFea1", "outerVespert",
          "vs_HLmyoLuc1", "Nearctic",
          "vs_myoDav1", "Myotis",
          "vs_HLpipPip1", "vs_HLlasBor1", "vs_HLnycHum2", "Vespertilioninae",
          "vs_HLmacSob1", "FoxLongTounge",
          "vs_HLeidHel2", "outerPeropodidae",
          "vs_HLeonSpe1", "Roussetinae"
         )')
  
  args = c('r=ComplexDietCentralAnalysisSimplifyAltReconstruction', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=SimplifiedDietConvert', 'v=F', 't=ER', 'n=ZoonomiaTip', 'z=0.01', 'l=F', 'e=T',
           'c=c(
              "Vertivore", "Omnivore", "Herbivore", "Invertivore"
            )', 
           'y=c(
          "vs_HLornAna3", "vs_HLtacAcu1",  "PreserveMonotremeBranches",
          "vs_HLdidVir1", "vs_HLgymLea1", "vs_HLpseCup1", "MarsupialTransitionPreservation",
          "vs_HLmyrTri1", "vs_HLchoDid1", "vs_HLchoHof3", "vs_HLproCap3", "AfrotheriaPreserveTransitions",
          "vs_HLursThi1", "vs_ursMar1", "vs_HLursArc1", "vs_HLailMel2", "UrsaPreserveTransition",
          "vs_lepWed1", "vs_HLmirAng2", "vs_HLphoVit1", "vs_HLeriBar1", "SealPreserveTransitions",
          "vs_HLodoRos1", "vs_HLcalUrs1", "vs_HLzalCal1", "SealSeaLionPreserveTransitions",
          "vs_HLmelCap1", "vs_HLgulGul1", "vs_HLneoVis1", "MustelidPreserveTransitions",
          "vs_HLlycPic2", "CanidPreserveTransision",
          "vs_HLgloMel1", "vs_HLpepEle1", "vs_HLturAdu1", "DolphinClade", "vs_orcOrc1", "vs_HLescRob1", "vs_HLlniGeo1", "amazonRiverDolphinFromYeast","vs_HLbalEde1", "vs_HLmegNov1", "vs_HLcynGun1", "CetaceaPreserveTransitions",
          "vs_HLmerUng1", "BankVoleTransition",
          "vs_HLeulMon1", "vs_HLeulFul1", "LemurTransition",
          "vs_panTro6", "vs_HLrhiRox2", "LangurClade", "vs_HLallNig1", "PrimateTransitionsContinued",
          "vs_HLeryPat1", "vs_chlSab2", "geunonClade", "vs_HLtheGel1", "PrimateTransitionsContinuedAgain",
          "vs_HLpapAnu5", "vs_HLmanSph1", "DrillMandrillClade", "vs_cerAty1", "Drilltransitions",
          "vs_HLmarFla1", "marmotClade", "DoormouseTransition",
          "vs_eulMac1", "vs_ponAbe3", "PrimateTransitions"
          
         )',
         'p=c(
          "vs_HLellTal1", "vs_HLellLut1", "vs_HLarvAmp1", "voleClade",
          "vs_HLhysCri1", "vs_HLthrSwi1", "vs_HLpetTyp1", "vs_hetGla2", "vs_chiLan1", "vs_HLdinBra1", "vs_HLcteSoc1", "vs_octDeg1", "vs_HLcoePre1", "vs_HLdasPun1", "vs_HLdolPat1", "gundiGuineaPigClade",
          "vs_HLoryGaz1", "vs_HLbeaHun1", "vs_HLkobLecLec1", "vs_HLkobLecLec1", "vs_HLmadKir1", "vs_HLneoPyg1", "vs_HLphiMax1", "vs_HLoreOre1", "vs_HLneoMos1", "vs_HLaepMel1", "vs_HLtraImb1", "Bovidae",
          "vs_HLhydIne1", "vs_HLmunMun1", "Cervidae",
          "vs_HLmurAurFea1", "outerVespert",
          "vs_HLmyoLuc1", "Nearctic",
          "vs_myoDav1", "Myotis",
          "vs_HLpipPip1", "vs_HLlasBor1", "vs_HLnycHum2", "Vespertilioninae",
          "vs_HLmacSob1", "FoxLongTounge",
          "vs_HLeidHel2", "outerPeropodidae",
          "vs_HLeonSpe1", "Roussetinae"
         )')
}





#args = c('r=Demo', 'v=T', 'm=C:/Users/mit221/AppData/Local/R/win-library/4.2/RERconverge/extdata/SubsetMammalGeneTrees.txt', 'd=Results/DemoMergedData.csv', 'a=DemoDietPhenotype', 'c=c("Carnivore", "Herbivore", "Omnivore")', 'n=demoTreeTipName')


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
useLiam = F
generateAlternates = F

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
    pruningCutoff = as.numeric(pruningCutoff)
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

  #use Liam Infrence
  if(!is.na(cmdArgImport('l'))){
    useLiam = as.logical(cmdArgImport('l'))
    message("Using Liam infrence -- DO NOT DO THIS UNLESS USING A SPECIFIC LIAM PHENOTYPE")
  }else{
    message("Liam infrence not used.")
  }

  #generate alternates
  if(!is.na(cmdArgImport('e'))){
    generateAlternates = as.logical(cmdArgImport('e'))
    message("Generating Alternates")
  }else{
    message("Not Generating Alternates.")
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
  relevantSpecies = manualAnnots[manualAnnots[[annotColumn]] %in% categoryList,]#remove all species which are not part of the specified categories
  if(useScreen){                                                                #if using a screening collumn 
    relevantSpecies = relevantSpecies[ relevantSpecies[[screenColumn]] %in% 1, ]  #remove all species not positive for that collumn 
  }
  relevantSpecies = relevantSpecies[!relevantSpecies[[nameColumn]] %in% "", ]   
}

# - Phenotype Vector - 
filteredSpecies = relevantSpecies[relevantSpecies[[nameColumn]] %in% speciesFilter, ] 
speciesNames = filteredSpecies[[nameColumn]]                                         #Exract the tip name of each species
speciesCategories = filteredSpecies[[annotColumn]]                              #extract the category of each species (in same order)

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


if(!useLiam){
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
}else{
  
  masterTree = mainTrees$masterTree
  
  nodesToAdd = c(455, 457, 471, 650, 492)
  names(nodesToAdd) = c("Mammalia", "Marsupalia", "Placentalia", "Chiroptera", "Primates")
  #phenToAdd = c("Insectivore", "Insectivore", "Insectivore", "Insectivore", "Omnivore")
  phenToAdd = c("Invertivore", "Invertivore", "Invertivore", "Invertivore", "Omnivore")
  #nodesToAdd = c(785)
  #names(nodesToAdd) = c("Placentalia")
  #phenToAdd = c("Carnivore")
  
  masterTreeAdded = masterTree
  for(i in 1:length(nodesToAdd)){
    M<-matchNodes(masterTree,masterTreeAdded,method="distances")
    masterTreeAdded<-bind.tip(masterTreeAdded,names(nodesToAdd)[i],edge.length=0,
                              where=M[which(M[,1]==as.numeric(nodesToAdd[i])),2])
  }
  mainTrees$masterTree = masterTreeAdded
  
  
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
  palette(c( "darkgreen", "darkblue","black", "red", "yellow"))
  
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
  
  #plotTreeCategorical(commonCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = commonMainTrees$masterTree)
  #plotTreeCategorical(categoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = mainTrees$masterTree)
  
  #plotTreeCategorical(commonCategoricalTreeExtraTip, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = commonMasterAdded)
  #plotTreeCategorical(categoricalTreeExtraTip, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = masterTreeAdded)
  
  plotTreeCategorical(commonCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore", "zMixedPredator"), master = commonMainTrees$masterTree)
  plotTreeCategorical(categoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore", "zMixedPredator"), master = mainTrees$masterTree)
  
  plotTreeCategorical(commonCategoricalTreeExtraTip, c("Herbivore", "Insectivore", "Omnivore", "Vertivore", "zMixedPredator"), master = commonMasterAdded)
  plotTreeCategorical(categoricalTreeExtraTip, c("Herbivore", "Insectivore", "Omnivore", "Vertivore", "zMixedPredator"), master = masterTreeAdded)
  dev.off()  
  
}


 
scientificCategoricalTree = ZoonomTreeNameToCommon(commonCategoricalTree, manualAnnotLocation = spreadSheetLocation, tipCol = "CommonName", scientific = T, scientificCol = "Scientific_Binomial", plot = F)

categoricalTreeFilename = paste(outputFolderName, filePrefix, "CategoricalTree.rds", sep="") #make a filename based on the prefix
saveRDS(categoricalTree, categoricalTreeFilename)                               #save the tree
categoricalCommonTreeFilename = paste(outputFolderName, filePrefix, "CategoricalCommonTree.rds", sep="") #make a filename based on the prefix
saveRDS(commonCategoricalTree, categoricalCommonTreeFilename)
scientificCategoricalTreeFilename = paste(outputFolderName, filePrefix, "CategoricalScientificTree.rds", sep="") #make a filename based on the prefix
saveRDS(scientificCategoricalTree, scientificCategoricalTreeFilename)

phenotypeVectorSaving = data.frame(names(phenotypeVector), names(commonPhenotypeVector), phenotypeVector)
directReadablePhenotypeVectorFilename = paste(outputFolderName, filePrefix, "CategoricalPhenotypeVector.csv",sep="")
write.csv(phenotypeVectorSaving, directReadablePhenotypeVectorFilename)



if(generateAlternates){
  useCladeLimit = F
  useStartDrop = T
  
    #Load in necessary script
    source("Src/Reu/CategoricalDropTip.R")
    mainTrees = readRDS(mainTreesLocation) #refreshes the main tree from any changes made during pruning step 
    
    #Load in the full tree with all species
    fullTreePrefix = paste0(filePrefix, "AllSpecies")
    fullTreeOutputFolder = paste0(outputFolderNameNoSlash, "AllSpecies", "/")
    
    #Used to convert v2 back to v1 to reuse all species
    #fullTrePrefixStub = substr(filePrefix, 1, (nchar(filePrefix)-1))
    #fullTreeOutputFolderStub = substr(outputFolderNameNoSlash, 1, (nchar(outputFolderNameNoSlash)-1))
    #fullTreePrefix = paste0(fullTrePrefixStub, "AllSpecies")
    #fullTreeOutputFolder = paste0(fullTreeOutputFolderStub, "AllSpecies", "/")
    
    
    tryCatch({
      allSpeciesTreeFilename = paste(fullTreeOutputFolder, fullTreePrefix, "CategoricalTree.rds", sep="") #make a filename based on the prefix
      allSpeciesTree = readRDS(allSpeciesTreeFilename)
      haveAllSpeciesTree = T
    }, 
    error=function(i){
      haveAllSpeciesTree = F
      message("No All Species tree found. Generating alternates, but not creating trees for them. ")
    }
    )
    
    #Figure out the number of short terminal branches in the main analysis 
    masterTree = mainTrees$masterTree
    mainAnalysisTree = mainTrees$masterTree
    tipsDroppedForMain = mainAnalysisTree$tip.label[which(!mainAnalysisTree$tip.label %in% categoricalTree$tip.label)]
    mainAnalysisTree = drop.tip(mainAnalysisTree, tipsDroppedForMain)
    terminalEdges = which(mainAnalysisTree$edge[,2] <= length(mainAnalysisTree$tip.label))
    mainNumShortBranches = length(which(mainAnalysisTree$edge.length[terminalEdges] <= pruningCutoff))
    
    
    masterTerminalEdges = which(masterTree$edge[,2] <= length(masterTree$tip.label))
    masterShortBranches = which(masterTree$edge.length <= pruningCutoff)
    
    masterShortTerminals = which((masterTree$edge[,2] <= length(masterTree$tip.label)) & (masterTree$edge.length <= pruningCutoff))
    masterShortTips = masterTree$edge[masterShortTerminals,2]
    
    
    
    phenotypeSizes = table(phenotypeVector)
    
    fullSpeciesNames = relevantSpecies[[nameColumn]]                                         #Exract the tip name of each species
    fullSpeciesCategories = relevantSpecies[[annotColumn]]                              #extract the category of each species (in same order)
    fullDataPhenotype = fullSpeciesCategories
    names(fullDataPhenotype) = fullSpeciesNames 
    if(!is.null(substitutions) & !all(is.na(substitutions))){
      for( i in 1:length(substitutions)){
        substitutePhenotypes = substitutions[[i]]
        message(paste("replacing", substitutePhenotypes[1], "with", substitutePhenotypes[2]))
        fullDataPhenotype = gsub(substitutePhenotypes[1], substitutePhenotypes[2], fullDataPhenotype)
      }
    }
    length(fullDataPhenotype)
    
    fullPhenotypeSizes = table(fullDataPhenotype)
    
    mastertreeDataPhenotype = fullDataPhenotype[names(fullDataPhenotype) %in% mainTrees$masterTree$tip.label]
    length(mastertreeDataPhenotype)
    masterTreePhenotypeSizes = table(mastertreeDataPhenotype)
    

    
    
    
    speciesSets = list()
    for(j in 1:length(phenotypeSizes)){
      speciesSets[[j]] = names(mastertreeDataPhenotype)[which(mastertreeDataPhenotype == names(phenotypeSizes)[j])]
    }
    numTips = sum(phenotypeSizes) # This is the number of tips, calculated here so it is done once instead of in the loop
    
    
    stableSpeciesSets = speciesSets
    stablePhenotypeSizes = phenotypeSizes
    
    

    
    alternateSets = list()
    i=1
    minShort = length(mainAnalysisTree$edge.length) 
    while(length(alternateSets) < 100){
      
      randomizedSpeciesSet = character()
      if(useCladeLimit){
        #reset the selection data to undo previous cuts
        phenotypeSizes = stablePhenotypeSizes
        speciesSets = stableSpeciesSets
        
        {
          # add preselected
          numPerClade = 3
          
          Bovidae = manualAnnots[manualAnnots$MSWC_Family %in% "Bovidae",][[nameColumn]]
          Pteropodidae = manualAnnots[manualAnnots$MSWC_Family %in% "Pteropodidae",][[nameColumn]]
          Cervidae = manualAnnots[manualAnnots$MSWC_Family %in% "Cervidae",][[nameColumn]]
          Cricetidae = manualAnnots[manualAnnots$MSWC_Family %in% "Cricetidae",][[nameColumn]]
          HystricognathiFamilies =  c("Caviidae", "Chinchillidae", "Ctenomyidae", "Dasyproctidae", "Dinomyidae", "Caviidae", "Bathyergidae", "Caviidae", "Hystricidae", "Myocastoridae", "Octodontidae", "Petromuridae", "Thryonomyidae", "Erethizontidae")
          Hystricognathi = manualAnnots[manualAnnots$MSWC_Family %in% HystricognathiFamilies,][[nameColumn]]
          
          Vespertilionidae = manualAnnots[manualAnnots$MSWC_Family %in% "Vespertilionidae",][[nameColumn]]
          speciesInLargeFamilies = c(Bovidae, Pteropodidae, Cervidae, Cricetidae, Hystricognathi, Vespertilionidae)

        selectedBovids = sample(Bovidae, numPerClade)
        selectedPteropodidae = sample(Pteropodidae, numPerClade)
        selectedCervidae = sample(Cervidae, numPerClade)
        selectedCricetidae = sample(Cricetidae, numPerClade)
        selectedHystricognathi = sample(Hystricognathi, numPerClade)
        
        selectedVespertilionidae = sample(Vespertilionidae, numPerClade)
        
        selectedHerbs = c(selectedBovids, selectedPteropodidae, selectedCervidae, selectedCricetidae, selectedHystricognathi)
        selectedInvert = c(selectedVespertilionidae)
        
        
        
        
        
        preselectedSpecies = c(selectedHerbs, selectedInvert)
        speciesToNotAutoselect = speciesInLargeFamilies
        
        
        }
        
        #Handle the preselcted species 
        {
          #Add to sample
          randomizedSpeciesSet = append(randomizedSpeciesSet, preselectedSpecies)
          
          #reduce category sizes 
          for(j in 1:length(speciesSets)){
            numPreselectedInCategory = length(which(preselectedSpecies %in% speciesSets[[j]]))
            phenotypeSizes[j] = phenotypeSizes[j] - numPreselectedInCategory
          }
          #remove from selectable list 
          for(j in 1:length(speciesSets)){
          speciesSets[[j]] = speciesSets[[j]][-!(speciesSets[[j]] %in% speciesToNotAutoselect)]
          }
        }
    }
       
      if(useStartDrop){
        #message(i)
        keptFraction = 2
        validRemoval = F
        regectedRemovals = 0
        
        while(validRemoval == F){
          speciesSets = stableSpeciesSets
          startRemoval = sample(masterShortTips, length(masterShortTips)/keptFraction)
          startRemoveTips = masterTree$tip.label[startRemoval]
          validSize = logical()
          
          #remove from selectable list 
          for(j in 1:length(speciesSets)){
            speciesSets[[j]] = speciesSets[[j]][!(speciesSets[[j]] %in% startRemoveTips)]
            #length(speciesSets[[j]])
            if(length(speciesSets[[j]]) > phenotypeSizes[j]){
              validSize[j] = T
            }else{
              validSize[j] = F
              #message(names(phenotypeSizes[j]))
              #message(paste(length(speciesSets[[j]]), "<", phenotypeSizes[j]))
            }
          }
          removePhens = fullDataPhenotype[names(fullDataPhenotype) %in% startRemoveTips] 
          #print(masterTreePhenotypeSizes)
          #print(phenotypeSizes)
          #print(table(removePhens))
          #length(speciesSets[[j]])
          if(all(validSize)){
            validRemoval = T
            #message(paste("Removal Accepted after rejections:", regectedRemovals))
            regectedRemovals = 0
          }else{
            regectedRemovals = regectedRemovals +1
            #message("Removal Rejected")
          }
          
        }

        
      }
      
      
      for(j in 1:length(phenotypeSizes)){
        #numberOfSpecies = phenotypeSizes[j]
        #speciesSet = mastertreeDataPhenotype[which(mastertreeDataPhenotype == names(phenotypeSizes)[j])]
        chosenSpecies = sample(speciesSets[[j]], phenotypeSizes[j])
        randomizedSpeciesSet = append(randomizedSpeciesSet, chosenSpecies)
      } 
      
      
      

      
      
      tipsToDrop = mainTrees$masterTree$tip.label[!mainTrees$masterTree$tip.label %in% randomizedSpeciesSet]
      testTree = drop.tip(mainTrees$masterTree, tipsToDrop)
      i=i+1
      if(i %% 10000 == 0){message(i)}
      terminalEdges = which(testTree$edge[,2] <= numTips)
      
      numberOfShortBranches = length(which(testTree$edge.length[terminalEdges] <= pruningCutoff))
      if(numberOfShortBranches < minShort){
        minShort = numberOfShortBranches
        message(i)
        message(numberOfShortBranches)
      }
      
      
      if(numberOfShortBranches <= mainNumShortBranches*1.5){
        message("Found Valid Alternate")
        message(i)
        alternateTips = testTree$tip.label
        alternateTips = list(alternateTips)
        alternateSets = append(alternateSets, alternateTips)
      }
    }
  
  saveRDS(alternateSets, paste0(outputFolderName, filePrefix, "AlternatePruningSpecies.rds"))
    
  
  if(!dir.exists(paste0("Output/", filePrefix, "/Alternates"))){                                      #Make output directory if it does not exist
    dir.create(paste0("Output/", filePrefix, "/Alternates"))
  }
  
  
  if(!haveAllSpeciesTree){
    for(i in 1:length(alternateSets)){
    currentSet = alternateSets[[i]]
    alternateFilePrefix = paste0("/Alternates/Alternate", i)
    
    
    currentPhenotypeVector = fullDataPhenotype[which(names(fullDataPhenotype) %in% currentSet)]
    currentCommonPhenotypeVector = currentPhenotypeVector
    names(currentCommonPhenotypeVector) = ZonomNameConvertVectorCommon(names(currentCommonPhenotypeVector), annotationLocation = spreadSheetLocation, tipColumn = nameColumn)
    
    
    currentPhenotypeVectorSaving = data.frame(names(currentPhenotypeVector), names(currentCommonPhenotypeVector), currentPhenotypeVector)
    currentDirectReadablePhenotypeVectorFilename = paste(outputFolderName, alternateFilePrefix, filePrefix, "CategoricalPhenotypeVector.csv",sep="")
    #write.csv(currentPhenotypeVectorSaving, currentDirectReadablePhenotypeVectorFilename) 
    currentPhenotypeVectorFilename = paste(outputFolderName, alternateFilePrefix, filePrefix, "CategoricalPhenotypeVector.rds",sep="") #make a filename based on the prefix
    saveRDS(currentPhenotypeVector, file = currentPhenotypeVectorFilename)                        #save the phenotype vector
    
    
    currentSpeciesFilterFilename =  paste(outputFolderName, alternateFilePrefix, filePrefix, "SpeciesFilter.rds",sep="") #set a filename for the species filter based on the prefix 
    currentSpeciesFilter = currentSet
    saveRDS(currentSpeciesFilter, currentSpeciesFilterFilename)
  }
  }
  
  if(haveAllSpeciesTree){
    
    
    commonMainTrees = mainTrees
    commonMainTrees$masterTree = ZoonomTreeNameToCommon(commonMainTrees$masterTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)
    
    for(i in 1:length(alternateSets)){
      message(i)
      currentSet = alternateSets[[i]]
      alternateFilePrefix = paste0("/Alternates/Alternate", i)
      
      
      currentPhenotypeVector = fullDataPhenotype[which(names(fullDataPhenotype) %in% currentSet)]
      currentCommonPhenotypeVector = currentPhenotypeVector
      names(currentCommonPhenotypeVector) = ZonomNameConvertVectorCommon(names(currentCommonPhenotypeVector), annotationLocation = spreadSheetLocation, tipColumn = nameColumn)
      
      
      currentPhenotypeVectorSaving = data.frame(names(currentPhenotypeVector), names(currentCommonPhenotypeVector), currentPhenotypeVector)
      currentDirectReadablePhenotypeVectorFilename = paste(outputFolderName, alternateFilePrefix, filePrefix, "CategoricalPhenotypeVector.csv",sep="")
      #write.csv(currentPhenotypeVectorSaving, currentDirectReadablePhenotypeVectorFilename)  
      
      currentPhenotypeVectorFilename = paste(outputFolderName, alternateFilePrefix, filePrefix, "CategoricalPhenotypeVector.rds",sep="") #make a filename based on the prefix
      saveRDS(currentPhenotypeVector, file = currentPhenotypeVectorFilename)                        #save the phenotype vector
      
      currentSpeciesFilterFilename =  paste(outputFolderName, alternateFilePrefix, filePrefix, "SpeciesFilter.rds",sep="") #set a filename for the species filter based on the prefix 
      currentSpeciesFilter = currentSet
      saveRDS(currentSpeciesFilter, currentSpeciesFilterFilename)
      
      
      tipsToRemove = allSpeciesTree$tip.label[!allSpeciesTree$tip.label %in% currentSet]
      currentTree = categoricalDropTip(allSpeciesTree, tipsToRemove)
      currentCommonSpeciesFilter = ZonomNameConvertVectorCommon(currentSet, annotationLocation = spreadSheetLocation, tipColumn = nameColumn)
      
      # - Categorical Tree - 
      currentCategoricalTree = currentTree
      currentCommonCategoricalTree = currentTree
      currentCommonCategoricalTree$tip.label = ZonomNameConvertVectorCommon(currentCommonCategoricalTree$tip.label, annotationLocation = spreadSheetLocation, tipColumn = nameColumn)
      
      if(useCladeLimit){
      
        treeImageFilename = paste(outputFolderName, alternateFilePrefix, filePrefix,"CategoricalTree.pdf", sep="") #make a filename based on the prefix
        pdf(treeImageFilename, height = length(currentPhenotypeVector)/18, width = 10)                     #make a pdf to store the plot, sized based on tree size
        palette(c( "darkgreen", "darkblue","black", "red"))

        plotTreeCategorical(currentCommonCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = commonMainTrees$masterTree)
        plotTreeCategorical(currentCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = mainTrees$masterTree)
        
        dev.off()  
      }
      
      
      currentCategoricalTreeFilename = paste(outputFolderName, alternateFilePrefix, filePrefix, "CategoricalTree.rds", sep="") #make a filename based on the prefix
      saveRDS(currentCategoricalTree, currentCategoricalTreeFilename)                               #save the tree
      currentCategoricalCommonTreeFilename = paste(outputFolderName, alternateFilePrefix, filePrefix, "CategoricalCommonTree.rds", sep="") #make a filename based on the prefix
      saveRDS(currentCommonCategoricalTree, currentCategoricalCommonTreeFilename)
      
    }
  }

  
  #Try to test for broken paths 
  alternateSets = readRDS(paste0(outputFolderName, filePrefix, "AlternatePruningSpecies.rds"))
  brokenPaths = integer()
  if(!dir.exists(paste0(outputFolderName, "/Alternates/Broken"))){                                      #Make output directory if it does not exist
    dir.create(paste0(outputFolderName, "/Alternates/Broken"))
  }
  for(i in 1:length(alternateSets)){
    message(i)
    currentSet = alternateSets[[i]]
    alternateFilePrefix = paste0("/Alternates/Alternate", i)
    currentCategoricalTreeFilename = paste(outputFolderName, alternateFilePrefix, filePrefix, "CategoricalTree.rds", sep="") #make a filename based on the prefix
    currentCommonCategoricalTreeFilename = paste(outputFolderName, alternateFilePrefix, filePrefix, "CategoricalCommonTree.rds", sep="") #make a filename based on the prefix
    currentPhenotypeVectorFilename = paste(outputFolderName, alternateFilePrefix, filePrefix, "CategoricalPhenotypeVector.rds",sep="") #make a filename based on the prefix
    currentSpeciesFilterFilename =  paste(outputFolderName, alternateFilePrefix, filePrefix, "SpeciesFilter.rds",sep="") #set a filename for the species filter based on the prefix 
    
    currentPhenotypeTree = readRDS(currentCategoricalTreeFilename)
    currentSpeciesFilter = readRDS(currentSpeciesFilterFilename)
    
    currentPathsObject = tree2Paths(currentPhenotypeTree, mainTrees, useSpecies = currentSpeciesFilter, categorical = TRUE)
    currentPathsFilename =  paste(outputFolderName, alternateFilePrefix, filePrefix, "CategoricalPathsFile.rds",sep="") #set a filename for the species filter based on the prefix 
    saveRDS(currentPathsObject, currentPathsFilename)
    
    print(table(currentPathsObject))
    
    if(length(table(currentPathsObject)) != length(phenotypeSizes)){
      message("Alternate is broken")
      file.rename(currentCategoricalTreeFilename, paste(outputFolderName, "Alternates/Broken/Alternate", i, filePrefix, "CategoricalTree.rds", sep=""))
      file.rename(currentCommonCategoricalTreeFilename, paste(outputFolderName, "Alternates/Broken/Alternate", i, filePrefix, "CategoricalCommonTree.rds", sep=""))
      file.rename(currentPhenotypeVectorFilename, paste(outputFolderName, "Alternates/Broken/Alternate", i, filePrefix, "CategoricalPhenotypeVector.rds", sep=""))
      file.rename(currentSpeciesFilterFilename, paste(outputFolderName, "Alternates/Broken/Alternate", i, filePrefix, "SpeciesFilter.rds", sep=""))
      
      brokenPaths = append(brokenPaths, i)
      
      
    }else{
      message("Paths fine, maintaining")
    }
  
  }
  brokenAlternateSet = alternateSets[brokenPaths]
  names(brokenAlternateSet) = brokenPaths
  brokenAlternatesFilename = paste(outputFolderName, "Alternates/Broken/", filePrefix, "BrokenAlternateSets.rds", sep="")
  saveRDS(brokenAlternateSet, brokenAlternatesFilename)
  
  backupAlternatesFilename = paste(outputFolderName, "Alternates/Broken/", filePrefix, "OriginalAlternateSets.rds", sep="")
  saveRDS(alternateSets, backupAlternatesFilename)
  
  
  #Exit the script if being run automatically before making any changes 
  if(!length(brokenPaths) ==0){
  stop("SOME OF THE ALTERNATES ARE PRODUCING BROKEN PATHS! Check the Broken folder to investigate further.")
  alternateSets = alternateSets[-brokenPaths]
  }
}

