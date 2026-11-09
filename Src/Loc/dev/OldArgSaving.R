#This is a file to save old copies of arguments, to maintain records while keeping the main scripts they are a part of cleaner. 



#---------------------------------------------------------------------
# --- New Merged Argument IDs  --- 
# --------------------------------------------------------------------
#Core arguments, must be supplied
# r = filePrefix                                                                REQUIRED: This is a prefix used to organize and separate files by analysis run. Always required. 
# v = <T or F>                                                                  This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# m = mainTreeFilename.txt or .rds                                              This sets the location of the maintrees file
# p = <T or F or C or A>                                                        This sets if the code should use permulated or unpermulated values. If set to "C" or "A", will use permulations for categorical values, which are stored in the main file. If set to "A", it will add permulations data, but use the unpermulationed data for significance. 
# g = gmtFileLocation.gmt                                                       This is the location of the main gmt file
# a = < ["b" or "binary"] or ["c" or "continuous"] or ["g" or "categorical"]>   This prefix is used to set the type of phenotype being supplied; it will detect automatically if not specified.



#Optional arguments, can typically self-determine if not specified
# s = "subdirectoryName"                                                        This is used to specify subdirectory(s) for the analysis to be run in. Primarily used for the components of categorical results. 
# n = < "diff" or "mean" or "last" >                                            This is used for continuous traits, to determine if the metic should be the difference between the nodes (diff), the mean(mean of the two nodes), or last(the downstream value). Note that Mean and Last are not phylogenetically independent, and do not have downstream processing. 
# l = <min.sp value>                                                            This sets the min.sp value to be used in the correlation. 
# w = saveData                                                                  This sets if the data should be saved, potentially overwriting old combinations
# c = significance cutoff                                                       This sets the cutoff for inclusion in the significant columns
# e = c("phenotypes", "to", "merge)                                             This sets which phenotype should be combined into a single phenotype in a merged analysis. 
# t = c("nameOfMergedPhenotype")                                                This sets the name of the now merged phenotype



#Overrides, which break from the normal flow and should only be used during manual runs or debugging. 
# f = speciesFilterText                                                         OPTIONAL OVERRIDE: This can be used to manually specify a species filer; leave blank for automatic
# u = phenotypeTreeFilename.txt or .rds                                         OPTIONAL OVERRIDE: This can be used to manually override the phenotype tree being used. For continuous analyses, this is the location of the trait vector.
# y = "correlationFileOverride.rds"                                             This can be used to manually set the correlation file. Used primarily to target pairwise comparisons of categorical phenotypes.
# z = "permulationPvalueFileLocation.rds"                                       This is a manual override to specify the script use a specific Permulation p-value file.









#---------------------------------------------------------------------
# --- Original Argument IDs   --- 
# --------------------------------------------------------------------
RunRERAndCorrelation.R
# -- Command arguments list
# -- Command arguments list
# r = filePrefix                                                               REQUIRED: This is a prefix used to organize and separate files by analysis run. Always required. 
# m = mainTreeFilename.txt or .rds                                             REQUIRED: This sets the location of the maintrees file
# v = <T or F>                                                                 This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# s = < ["b" or "binary"] or ["c" or "continuous"] or ["g" or "categorical"]>  This prefix is used to set the type of phenotype being supplied; it will detect automatically if not specified.
# c = < "diff" or "mean" or "last" >                                           This is used for continuous traits, to determine if the metic should be the difference between the nodes (diff), the mean(mean of the two nodes), or last(the downstream value). Note that Mean and Last are not phylogenetically independent, and do not have downstream processing. 
# l = <min.sp value>                                                           This sets the min.sp value to be used in the correlation. 
# f = speciesFilterText                                                        OPTIONAL OVERRIDE: This can be used to manually specify a species filer; leave blank for automatic
# p = phenotypeTreeFilename.txt or .rds                                        OPTIONAL OVERRIDE: This can be used to manually override the phenotype tree being used. For continuous analyses, this is the location of the trait vector.



EnrichmentAnalysisNew.R
# -- Command arguments list
# -- Command arguments list
# r = filePrefix                            This is a prefix used to organize and separate files by analysis run. Always required. 
# v = <T or F>                              This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# m = gmtFileLocation.gmt                   This is the location of the main gmt file
# p = <T or F or C>                         This sets if the code should use permulated or unpermulated values. If set to "C", will use permulations for categorical values, which are stored in the main file.
# s = "subdirectoryName"                    This is used to specify a subdirectory for the analysis to be run in. Primarily used for the components of categorical results. 
# c = "correlationFileOverride.rds"         This can be used to manually set the correlation file. Used primarily to target pairwise comparisons of categorical phenotypes.
# f = "permulationPvalueFileLocation.rds"   This is a manual override to specify the script use a specific Permulation p-value file.
#If using any file other than "CombinedPrunedFastAll" with no run instance number, it must be specified manually.


CategoricalCombinePairwiseResults.R
# -- Command arguments list
# r = filePrefix                            This is a prefix used to organize and separate files by analysis run. Always required. 
# v = <T or F>                              This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# g = gene set                              This sets the gene set (name, no file) that the script should be run on, or "gene" for use on genes. 
# s = saveData                                                                  This sets if the data should be saved, potentially overwriting old combinations
# c = significance cutoff                                                       This sets the cutoff for inclusion in the significant columns
# l = <T or F or C or A>                                                        This sets if the code should use permulated or unpermulated values. If set to "C" or "A", will use permulations for categorical values, which are stored in the main file. If set to "A", it will add permulations data, but use the unpermulationed data for significance. 








MakeCategoricalPhenotype.R
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




#---------------------------------------------------------------------
# --- MakeCategoricalPhenotype.R (new) --- 
# --------------------------------------------------------------------


{
  args = c('r=ComplexDietCentralAnalysisSimplifyStrictPred2', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=SimplifiedDietConvertStrictPred', 'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01', 'l=T', 'e=T',
           'c=c(
              "Vertivore", "Omnivore", "Herbivore", "Invertivore", "zMixedPredator"
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
  
  args = c('r=ComplexDietViewHerbAllSpecies', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=ComplexHerbDiet', 'v=T', 't=ER', 'n=ZoonomiaTip', 'l=T', 'e=T',
           'c=c(
              "Vertivore", "Omnivore", "Invertivore", "Frugivore", "Generalist Herbivore", "Granivore", "Mixed Herbivore", "Nectivore", "Strict Frugivore", "Strict Nectivore"
            )',
           'u=list(
            c("Mixed Herbivore", "Generalist Herbivore"),
            c("Strict Frugivore", "Frugivore"),
            c("Strict Nectivore", "Nectivore")
           )'
  )
}
{
  
  args = c('r=ComplexDietHighSugarAllSpecies', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=ComplexHerbDiet', 'v=T', 't=ER', 'n=ZoonomiaTip', 'l=T', 'e=F',
           'c=c(
              "Vertivore", "Omnivore", "Invertivore", "Frugivore", "Generalist Herbivore", "Granivore", "Mixed Herbivore", "Nectivore", "Strict Frugivore", "Strict Nectivore"
            )',
           'u=list(
            c("Mixed Herbivore", "HerbivoreGeneralist"), c("Generalist Herbivore", "HerbivoreGeneralist"), c("Granivore", "HerbivoreGeneralist"), 
            c("Strict Frugivore", "SugarHerbviore"), c("Strict Nectivore", "SugarHerbviore"),c("Frugivore", "SugarHerbviore"), c("Nectivore", "SugarHerbviore")
           )',
           'b=c("darkgreen", "darkblue","black", "yellow", "red")'
  )
  
  args = c('r=ComplexDietHighSugar', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=ComplexHerbDiet', 'v=T', 't=ER', 'n=ZoonomiaTip', 'l=T', 'e=T','z=0.01',
           'c=c(
              "Vertivore", "Omnivore", "Invertivore", "Frugivore", "Generalist Herbivore", "Granivore", "Mixed Herbivore", "Nectivore", "Strict Frugivore", "Strict Nectivore"
            )',
           'u=list(
            c("Mixed Herbivore", "HerbivoreGeneralist"), c("Generalist Herbivore", "HerbivoreGeneralist"), c("Granivore", "HerbivoreGeneralist"), 
            c("Strict Frugivore", "SugarHerbviore"), c("Strict Nectivore", "SugarHerbviore"),c("Frugivore", "SugarHerbviore"), c("Nectivore", "SugarHerbviore")
           )',
           'b=c("darkgreen", "darkblue","black", "yellow", "red")',
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
         )'
  )
}


#---------------------------------------------------------------------
# --- MakeCategoricalPhenotype.R  --- 
# --------------------------------------------------------------------
{
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
{
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
}
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

{args = c('r=CategoricalInsVertivoreTreePruningTest', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01',
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
}
{
  args = c('r=CategoricalInsVertivoreTreeLiamInference', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01', 'l=T',
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
}

{
  
  args = c('r=CategoricalInsVertivoreTreeDuplicateLiamInference', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01', 'l=T',
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
          "vs_HLgloMel1", "vs_HLpepEle1", "vs_HLturAdu1", "DolphinTransistions", "vs_orcOrc1", "vs_HLescRob1", "vs_HLbalEde1", "vs_HLmegNov1", "vs_HLcynGun1", "CetaceaPreserveTransitions",
          "vs_HLmerUng1", "BankVoleTransition",
          "vs_HLeulMon1", "vs_HLeulFul1", "LemurTransition",
          "vs_eulMac1", "vs_HLeulFla1", "vs_ponAbe3", "PrimateTransitions",
          "vs_panTro6", "vs_HLrhiRox2", "LangurClade", "vs_HLallNig1", "PrimateTransitionsContinued",
          "vs_HLeryPat1", "vs_chlSab2", "geunonClade", "vs_HLtheGel1", "PrimateTransitionsContinuedAgain",
          "vs_HLpapAnu5", "vs_HLmanSph1", "DrillMandrillClade", "vs_cerAty1", "Drilltransitions",
          "vs_HLmarFla1", "marmotClade", "DoormouseTransition",
          "vs_HLtraJav1", "MouseDeerUsedInOtherAnalysis",
          "vs_mm10", "humans",
          "vs_HLlniGeo1", "amazonRiverDolphinFromYeast",
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
  
  args = c('r=CategoricalInsVertivoreTreeFullFamilyLiamInference', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01', 'l=T',
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
          "vs_HLgloMel1", "vs_HLpepEle1", "vs_HLturAdu1", "DolphinTransistions", "vs_orcOrc1", "vs_HLescRob1", "vs_HLbalEde1", "vs_HLmegNov1", "vs_HLcynGun1", "CetaceaPreserveTransitions",
          "vs_HLmerUng1", "BankVoleTransition",
          "vs_HLeulMon1", "vs_HLeulFul1", "LemurTransition",
          "vs_eulMac1", "vs_HLeulFla1", "vs_ponAbe3", "PrimateTransitions",
          "vs_panTro6", "vs_HLrhiRox2", "LangurClade", "vs_HLallNig1", "PrimateTransitionsContinued",
          "vs_HLeryPat1", "vs_chlSab2", "geunonClade", "vs_HLtheGel1", "PrimateTransitionsContinuedAgain",
          "vs_HLpapAnu5", "vs_HLmanSph1", "DrillMandrillClade", "vs_cerAty1", "Drilltransitions",
          "vs_HLmarFla1", "marmotClade", "DoormouseTransition"
         )',
         'p=c(
          "vs_HLellTal1", "vs_HLellLut1", "vs_HLarvAmp1","vs_HLmicAgr2", "vs_HLmyoGla2", "vs_HLondZib1", "voleClade", "6Herb",
          "vs_HLmusSpi1", "vs_HLmusCar1", "vs_HLmasCou1", "vs_HLmusPah1", "vs_HLratNor7", "vs_HLarvNil1", "vs_mm10", "mouseClade", "7Omni,1Herb",
          "vs_HLhysCri1", "vs_HLthrSwi1", "vs_HLpetTyp1", "vs_hetGla2", "vs_chiLan1", "vs_HLdinBra1", "vs_HLcteSoc1", "vs_octDeg1", "vs_HLcoePre1", "vs_HLdasPun1", "vs_HLdolPat1", "vs_HLmyoCoy1", "vs_HLhydHyd1", "vs_HLcavTsc1", "gundiGuineaPigClade", "14Herb",
          "vs_HLoryGaz1", "vs_HLbeaHun1", "vs_HLkobLecLec1", "vs_HLkobLecLec1", "vs_HLmadKir1", "vs_HLneoPyg1", "vs_HLphiMax1", "vs_HLoreOre1", "vs_HLneoMos1", "vs_HLaepMel1", "vs_HLtraImb1",  "vs_bisBis1", "vs_HLoviNivLyd1", "vs_HLproPrz1", "Bovidae", "13Herb",
          "vs_HLhydIne1", "vs_HLmunMun1", "vs_HLodoHem1", "vs_HLantAme1", "vs_HLgirCam1" "Cervidae", "5Herb",
          "vs_HLmurAurFea1", "Murina", "vs_HLmyoLuc1", "Nearctic", "vs_myoDav1", "vs_HLmyoMyo6", "vs_HLmyoSep1", "vs_HLmyoLuc1", "Myotis", "vs_HLpipPip1", "vs_HLlasBor1", "vs_HLnycHum2", "vs_eptFus1", "Vespertilioninae","Vespertilionidae", "10Inse",
          "vs_HLmacSob1", "vs_HLpteGig1", FoxLongTounge", "vs_HLeidHel2", "vs_HLcynBra1", outerPeropodidae", "vs_HLeonSpe1", "vs_HLrouLes1", "Roussetinae", "Peropodidae", "6Herb",
         )')
  
  args = c('r=CategoricalInsVertivoreTreeFamilyAgnostictLiamInference', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01', 'l=T',
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
          "vs_HLgloMel1", "vs_HLpepEle1", "vs_HLturAdu1", "DolphinTransistions", "vs_orcOrc1", "vs_HLescRob1", "vs_HLbalEde1", "vs_HLmegNov1", "vs_HLcynGun1", "CetaceaPreserveTransitions",
          "vs_HLmerUng1", "BankVoleTransition",
          "vs_HLeulMon1", "vs_HLeulFul1", "LemurTransition",
          "vs_eulMac1", "vs_HLeulFla1", "vs_ponAbe3", "PrimateTransitions",
          "vs_panTro6", "vs_HLrhiRox2", "LangurClade", "vs_HLallNig1", "PrimateTransitionsContinued",
          "vs_HLeryPat1", "vs_chlSab2", "geunonClade", "vs_HLtheGel1", "PrimateTransitionsContinuedAgain",
          "vs_HLpapAnu5", "vs_HLmanSph1", "DrillMandrillClade", "vs_cerAty1", "Drilltransitions",
          "vs_HLmarFla1", "marmotClade", "DoormouseTransition"
         )')  
  
  
  
  args = c('r=CategoricalInsVertivoreTreeNoYeastLiamInference', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01', 'l=T',
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
          "vs_HLmusSpi1", "vs_HLmusCar1", "vs_HLmasCou1", "vs_HLmusPah1", "mouseClade",
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
  
  args = c('r=CategoricalInsVertivoreTreeNoYeastCarnivoreLiamInference', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01', 'l=T',
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
          "vs_HLmusSpi1", "vs_HLmusCar1", "vs_HLmasCou1", "vs_HLmusPah1", "mouseClade",
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
  
  args = c('r=CategoricalInsVertivoreTreeNoManualLiamInference', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01', 'l=T',
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
  
  args = c('r=CategoricalInsVertivoreTreeNoManualCarnivoreLiamInference', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01', 'l=T',
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
          )')
  
  args = c('r=PredatorFullTree', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 'v=T', 't=ER', 'n=ZoonomiaTip', 'l=T',
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
  
  args = c('r=DuplicatePredatorFullTree', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 'v=T', 't=ER', 'n=ZoonomiaTip', 'l=T',
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
  
  args = c('r=TestFullTreeFiltered', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 'v=T', 't=ER', 'n=ZoonomiaTip', 'l=T',
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
  
  args = c('r=SplitBothFullTree', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 'v=T', 't=ER', 'n=ZoonomiaTip', 'l=T',
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
            c("H-Frugivore", "H-High-sugar-plants-Eater"), 
            c("H-Nectarivore", "H-High-sugar-plants-Eater"), 
            c("H-High-sugar-plants-Eater", "SugarHigh"),
            c("H-Granivore", "H-Low-sugar-plants-Eater"), c("H-Nonspecific-Herbivore", "H-Low-sugar-plants-Eater"), 
            c("H-Low-sugar-plants-Eater", "SugarLow"), 
            c("H-All-plants-Eater", "Omnivore"),
            c("O-Generalist", "Omnivore")
          )')
  args = c('r=SplitBothPruned', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 'v=T', 't=ER', 'n=ZoonomiaTip', 'l=T','z=0.01',
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
            c("H-Frugivore", "H-High-sugar-plants-Eater"), 
            c("H-Nectarivore", "H-High-sugar-plants-Eater"), 
            c("H-High-sugar-plants-Eater", "SugarHigh"),
            c("H-Granivore", "H-Low-sugar-plants-Eater"), c("H-Nonspecific-Herbivore", "H-Low-sugar-plants-Eater"), 
            c("H-Low-sugar-plants-Eater", "SugarLow"), 
            c("H-All-plants-Eater", "SugarHigh"),
            c("O-Generalist", "Omnivore")
          )')
}

{  args = c('r=ComplexDietCentralAnalysisAllSpecies', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 'v=T', 't=ER', 'n=ZoonomiaTip', 'l=T',
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
  
  args = c('r=ComplexDietCentralAnalysis', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01', 'l=T', 'e=T',
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
  
  
  args = c('r=ComplexDietCentralAnalysisNoRodent', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01', 'l=T', 'e=T',
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
          "vs_HLellTal1", "vs_HLellLut1", "vs_HLarvAmp1", "voleClade", "Cricetidae",
          "vs_HLoryGaz1", "vs_HLbeaHun1", "vs_HLkobLecLec1", "vs_HLkobLecLec1", "vs_HLmadKir1", "vs_HLneoPyg1", "vs_HLphiMax1", "vs_HLoreOre1", "vs_HLneoMos1", "vs_HLaepMel1", "vs_HLtraImb1", "Bovidae",
          "vs_HLhydIne1", "vs_HLmunMun1", "Cervidae",
          "vs_HLmurAurFea1", "outerVespert", "vs_HLmyoLuc1", "Nearctic", "vs_myoDav1", "Myotis", "vs_HLpipPip1", "vs_HLlasBor1", "vs_HLnycHum2", "Vespertilioninae",
          "vs_HLmacSob1", "FoxLongTounge", "vs_HLeidHel2", "outerPeropodidae", "vs_HLeonSpe1", "Roussetinae", "Peropodidae"
         )')
  
  args = c('r=ComplexDietCentralAnalysisNoFamilyPrune', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01', 'l=T', 'e=T',
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
          
         )')
}
{
  args = c('r=CladeBinaryBovidae', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'v=T', 't=ER', 'n=ZoonomiaTip', 'l=F', 'e=T',
           'c=c(0, 1)',
           'a=isBovidae','s=isInAnalysisWithFullFamilies') 
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
  args = c('r=CladeBinaryCervidae', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'v=T', 't=ER', 'n=ZoonomiaTip', 'l=F', 'e=T',
           'c=c(0, 1)',
           'a=isCervidae','s=isInAnalysisWithFullFamilies') 
  
  args = c('r=ComplexDietCentralAnalysisEqualLengthInference', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01', 'l=F', 'e=T',
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
{args = c('r=HarshalCategoricalRERNew', 'm=Data/HarshalFakeMainTrees.rds', 
          'd=Data/VGP_mammals_Diet.csv', 
          'n=Accession',
          'a=trophic_level',
          'c=c("Herbivore", "Carnivore", "Omnivore")',
          'v=T', 't=ER', 'l=T')
  
  args = c('r=PosterTreeNew', 'm=Data/LeahFakeMainTrees.rds', 
           'd=Data/LeahPosterData.csv', 
           'n=ScientificName',
           'a=Four_Diet',
           'c=c("Herbivore", "Vertivore", "Invertivore", "Omnivore")',
           'v=T', 't=ER', 'l=F')
}
{
  args = c('r=ComplexDietCentralAnalysis', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01', 'l=T', 'e=T',
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
{ 
  args = c('r=ComplexDietCentralAnalysisSimplifyAllSpecies', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=SimplifiedDietConvert', 'v=T', 't=ER', 'n=ZoonomiaTip', 'l=T',
           'c=c(
              "Vertivore", "Omnivore", "Herbivore", "Invertivore"
            )'
  )
  
  args = c('r=ComplexDietCentralAnalysisSimplify2', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=SimplifiedDietConvert', 'v=F', 't=ER', 'n=ZoonomiaTip', 'z=0.01', 'l=T', 'e=T',
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
  
  args = c('r=ComplexDietCentralAnalysisSimplifyNoScav', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=SimplifiedDietConvertNoScav', 'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01', 'l=T', 'e=T',
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
  
  args = c('r=ComplexDietCentralAnalysisPiscFix', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=TranslatedDerek5050Fix', 'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01', 'l=T', 'e=T',
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
  
  args = c('r=ComplexDietCentralAnalysisNewAlternates', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=DerekDietClassification90InsVertivoreSorting', 'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01', 'l=T', 'e=T',
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
  
  
  args = c('r=ComplexDietCentralAnalysisNoFishbat', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=TranslatedDerekNoFishbat', 'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01', 'l=T', 'e=T',
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
  
  
  args = c('r=ComplexDietCentralAnalysisSimplifyEqualDropAllSpecies', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=SimplifiedDietConvertEqualDrop', 'v=T', 't=ER', 'n=ZoonomiaTip', 'l=T',
           'c=c(
              "Vertivore", "Omnivore", "Herbivore", "Invertivore", "zMixedPredator"
            )'
  )
  
  args = c('r=ComplexDietCentralAnalysisSimplifyEqualDrop', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=SimplifiedDietConvertEqualDrop', 'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01', 'l=T', 'e=T',
           'c=c(
              "Vertivore", "Omnivore", "Herbivore", "Invertivore", "zMixedPredator"
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
  
  args = c('r=ComplexDietCentralAnalysisSimplifyStrictPredAllSpecies', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=SimplifiedDietConvertStrictPred', 'v=T', 't=ER', 'n=ZoonomiaTip', 'l=T',
           'c=c(
              "Vertivore", "Omnivore", "Herbivore", "Invertivore", "zMixedPredator"
            )'
  )
  
  args = c('r=ComplexDietCentralAnalysisSimplifyStrictPred', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=SimplifiedDietConvertStrictPred', 'v=T', 't=ER', 'n=ZoonomiaTip', 'z=0.01', 'l=T', 'e=T',
           'c=c(
              "Vertivore", "Omnivore", "Herbivore", "Invertivore", "zMixedPredator"
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
  
  args = c('r=ComplexDietCentralAnalysisSimplifySubT70AllSpecies', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=SimplifiedDietConvert', 'v=F', 't=ER', 'n=ZoonomiaTip', 'l=T', 'e=T',
           'c=c(
              "Vertivore", "Omnivore", "Herbivore", "Invertivore", "zMixedPredator"
            )')
  
  args = c('r=ComplexDietCentralAnalysisSimplifySubT70', 'm=data/zoonomiaAllMammalsTrees.rds', 'd=Data/mergedData.csv', 'a=SimplifiedDietConvert', 'v=F', 't=ER', 'n=ZoonomiaTip', 'z=0.01', 'l=T', 'e=T',
           'c=c(
              "Vertivore", "Omnivore", "Herbivore", "Invertivore", "zMixedPredator"
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
}




#---------------------------------------------------------------------
# --- CombineCategoricalPairwiseResults.R  --- 
# --------------------------------------------------------------------

args = c("r=CategoricalInsvertivoreTree", "p=NULL", "g=gene")
args = c("r=CategoricalInsvertivoreTree", "p=NULL", "g=KeggReactome")

args = c("r=CategoricalInsvertivoreTreeLiamInference", "p=NULL", "g=gene")
args = c("r=CategoricalInsvertivoreTreeLiamInference", "p=NULL", "g=KeggReactome")

args = c("r=CategoricalInsvertivoreTreeCarnivoreLiamInference", "p=NULL", "g=gene", "c=0.05", "s=F")
args = c("r=CategoricalInsvertivoreTreeCarnivoreLiamInference", "p=NULL", "g=KeggReactome", "c=0.1", "s=T")

args = c("r=CategoricalInsvertivoreTreeLiamInference", "p=NULL", "g=MGI_Mammalian_Phenotype_Level_4", "c=0.05", "s=T")
args = c("r=CategoricalInsvertivoreTreeLiamInference", "p=NULL", "g=GO_Biological_Process_2023", "c=0.05", "s=T")
args = c("r=CategoricalInsvertivoreTreeLiamInference", "p=NULL", "g=DisGeNET", "c=0.05", "s=T")
args = c("r=CategoricalInsvertivoreTreeLiamInference", "p=NULL", "g=tissue_specific", "c=0.05", "s=T")
args = c("r=CategoricalInsvertivoreTreeLiamInference", "p=NULL", "g=EnrichmentHsSymbolsFile2", "c=0.05", "s=T")
args = c("r=CategoricalInsvertivoreTreeLiamInference", "p=NULL", "g=gene", "s=T", "l=T")
args = c("r=CategoricalInsvertivoreTreeLiamInference", "p=NULL", "g=KeggReactome", "s=T", "l=F")
args = c("r=CategoricalInsvertivoreTreeLiamInference", "p=NULL", "g=KeggReactome", "s=T", "l=T")
args = c("r=CategoricalInsvertivoreTreeLiamInference", "p=NULL", "g=KeggReactome", "s=T", "l=T")


args = c("r=CategoricalInsvertivoreTreeNoYeastLiamInference", "p=NULL", "g=gene", "s=T", "l=F")
args = c("r=CategoricalInsvertivoreTreeNoYeastLiamInference", "p=NULL", "g=KeggReactome", "s=T", "l=F")

args = c("r=CategoricalInsvertivoreTreeNoManualLiamInference", "p=NULL", "g=gene", "s=T", "l=F")
args = c("r=CategoricalInsvertivoreTreeNoManualLiamInference", "p=NULL", "g=KeggReactome", "s=T", "l=F")

args = c("r=PredatorFullTree", "p=NULL", "g=gene", "s=T", "l=F")
args = c("r=PredatorFullTree", "p=NULL", "g=KeggReactome", "s=T", "l=F")


args = c("r=CategoricalInsvertivoreTreeFamilyAgnostictLiamInference", "p=NULL", "g=gene", "s=T", "l=F")
args = c("r=CategoricalInsvertivoreTreeFamilyAgnostictLiamInference", "p=NULL", "g=KeggReactome", "s=T", "l=F")


args = c("r=ComplexDietCentralAnalysis", "p=NULL", "g=gene", "s=T", "l=A")
args = c("r=ComplexDietCentralAnalysis", "p=NULL", "g=KeggReactome", "s=T", "l=A")


args = c("r=ComplexDietCentralAnalysisSimplify", "p=NULL", "g=gene", "s=T", "l=F")
args = c("r=ComplexDietCentralAnalysisSimplify", "p=NULL", "g=KeggReactome", "s=T", "l=F")

args = c("r=ComplexDietCentralAnalysisSimplifyNoScav", "p=NULL", "g=gene", "s=T", "l=F")
args = c("r=ComplexDietCentralAnalysisSimplifyNoScav", "p=NULL", "g=KeggReactome", "s=T", "l=F")

args = c("r=ComplexDietCentralAnalysisPiscFix", "p=NULL", "g=gene", "s=T", "l=F")
args = c("r=ComplexDietCentralAnalysisPiscFix", "p=NULL", "g=KeggReactome", "s=T", "l=F")

args = c("r=ComplexDietCentralAnalysisNoFishbat", "p=NULL", "g=gene", "s=T", "l=F")
args = c("r=ComplexDietCentralAnalysisNoFishbat", "p=NULL", "g=KeggReactome", "s=T", "l=F")


args = c("r=ComplexDietCentralAnalysisSimplifyEqualDrop", "g=KeggReactome", "s=T", "l=F", 'p=c("Invertivore-Vertivore","Herbivore-Vertivore","Invertivore-Omnivore","Herbivore-Omnivore","Herbivore-Invertivore","Omnivore-Vertivore","Carnivore-Herbivore","Carnivore-Omnivore","Overall")')

args = c("r=ComplexDietCentralAnalysisSimplifyEqualDrop", "g=gene", "s=T", "l=F")
args = c("r=ComplexDietCentralAnalysisSimplifyEqualDrop", "g=KeggReactome", "s=T", "l=F")


args = c("r=ComplexDietCentralAnalysisSimplify2", "g=MGI_Mammalian_Phenotype_Level_4", "c=0.05", "s=T", 'p=c("Invertivore-Vertivore","Herbivore-Vertivore","Invertivore-Omnivore","Herbivore-Omnivore","Herbivore-Invertivore", "Carnivore-Herbivore")')
args = c("r=ComplexDietCentralAnalysisSimplify2", "g=GO_Biological_Process_2023", "c=0.05", "s=T", 'p=c("Invertivore-Vertivore","Herbivore-Vertivore","Invertivore-Omnivore","Herbivore-Omnivore","Herbivore-Invertivore", "Carnivore-Herbivore")')
args = c("r=ComplexDietCentralAnalysisSimplify2", "g=DisGeNET", "c=0.05", "s=T", 'p=c("Invertivore-Vertivore","Herbivore-Vertivore","Invertivore-Omnivore","Herbivore-Omnivore","Herbivore-Invertivore", "Carnivore-Herbivore")')
args = c("r=ComplexDietCentralAnalysisSimplify2", "g=EnrichmentHsSymbolsFile2", "c=0.05", "s=T", 'p=c("Invertivore-Vertivore","Herbivore-Vertivore","Invertivore-Omnivore","Herbivore-Omnivore","Herbivore-Invertivore", "Carnivore-Herbivore")')

args = c("r=ComplexDietCentralAnalysisSimplify2", "g=tissue_specific", "c=0.05", "s=T", 'p=c("Invertivore-Vertivore","Herbivore-Vertivore","Invertivore-Omnivore","Herbivore-Omnivore","Herbivore-Invertivore", "Carnivore-Herbivore")')




#---------------------------------------------------------------------
# --- EnrichmentAnalysisNew.R  --- 
# --------------------------------------------------------------------



#'m=c("Data/tissue_specific.gmt", "Data/GSEA-c5-HsSymbols.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")'
#'m=c("Data/DisGeNET.gmt", "Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt")'
args = c('r=CategoricalDiet4Phen', "s=_Omnivore-Carnivore", 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'v=T', 'p=T') #This is a debug argument set. It is used to set arguments locally, when not running the code through a bash script.
args = c('r=LiverExpression',  'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'v=T', 'p=F') #This is a debug argument set. It is used to set arguments locally, when not running the code through a bash script.
args = c('r=EcholocationUpdate',  'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'v=T', 'p=T') #This is a debug argument set. It is used to set arguments locally, when not running the code through a bash script.
args = c('r=CategoricalDiet4Phen', 's=c("_Omnivore-Herbivore", "Carnivore-Herbivore", "_Omnivore-Insectivore", "Carnivore-Insectivore", "Herbivore-Insectivore", "_Omnivore-carnivore")', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'v=T', 'p=F') #This is a debug argument set. It is used to set arguments locally, when not running the code through a bash script.
args = c('r=CategoricalDiet5Phen', 's=c("_Omnivore-Herbivore", "Carnivore-Herbivore", "_Omnivore-Insectivore", "Carnivore-Insectivore", "Herbivore-Insectivore", "_Omnivore-Piscivore", "Carnivore-Piscivore", "Herbivore-Piscivore", "Insectivore-Piscivore", "_Omnivore-carnivore")', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'v=T', 'p=F') #This is a debug argument set. It is used to set arguments locally, when not running the code through a bash script.
args = c('r=CategoricalDiet3Phen', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 's=c("_Omnivore-Herbivore", "Carnivore-Herbivore", "_Omnivore-Carnivore", "Overall")', 'p=T')
args = c('r=LiverExpression3', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=T', "c=LiverExpression3CorrelationDataPermulatedNamesConverted.rds")
args = c('r=EcholocationUpdate2',  'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'v=T', 'p=T') #This is a debug argument set. It is used to set arguments locally, when not running the code through a bash script.
args = c('r=CVHRemake',  'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'v=T', 'p=T') #This is a debug argument set. It is used to set arguments locally, when not running the code through a bash script.
args = c('r=CategoricalDiet3Phen', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 's=c("_Omnivore-Herbivore", "Carnivore-Herbivore", "_Omnivore-Carnivore", "Overall")', 'p=C')
args = c('r=CVHRemake',  'm=c("Data/customGeneSet.gmt")', 'v=T', 'p=T') #This is a debug argument set. It is used to set arguments locally, when not running the code through a bash script.
args = c('r=CategoricalDiet3Phen', 'm=c("Data/customGeneSet.gmt")', 's=c("_Omnivore-Herbivore", "Carnivore-Herbivore", "_Omnivore-Carnivore", "Overall")', 'p=C')
args = c('r=LiverExpression3', 'm=c("Data/customGeneSet.gmt")', 'p=T', "c=LiverExpression3CorrelationDataPermulatedNamesConverted.rds")


args = c('r=NewHiller4Phen', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("_Omnivore-Carnivore", "_Omnivore-Herbivore", "Carnivore-Herbivore", "_Omnivore-Insectivore", "Carnivore-Insectivore", "Herbivore-Insectivore", "Overall")')
args = c('r=NewHillerTestSupraPrimates', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("0-1", "Overall")')
args = c('r=NewHillerTestSupraPrimates', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("0-1", "Overall")')
args = c('r=NewHillerTestSupraPrimates', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("0-1", "Overall")')

args = c('r=MaturityLifespanPercent', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F')
args = c('r=MaturityLifespanPercent', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=T', 'f=MaturityLifespanPercentPermulationPValue.rds')
args = c('r=MaturityLogRaw', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F')
args = c('r=PankajBodysize', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F')


args = c('r=MaturityLifespanPercent', 'm=c("Data/YifanGenesets.gmt")', 'p=T', 'f=MaturityLifespanPercentPermulationPValue.rds')
args = c('r=MaturityLogRaw', 'm=c("Data/YifanGenesets.gmt")', 'p=F')



args = c('r=CategoricalMobivoreTree', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("Omnivore-Vertivore", "Insectivore-Vertivore", "Herbivore-Vertivore", "Insectivore-Omnivore", "Herbivore-Omnivore", "Herbivore-Insectivore", "Overall")' )
args = c('r=CategoricalCarnivoreTree', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("Carnivore-Herbivore", "Carnivore-Omnivore", "Herbivore-Omnivore", "Overall")' )

args = c('r=CategoricalInsVertivoreTree', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("Omnivore-Vertivore", "Insectivore-Vertivore", "Herbivore-Vertivore", "Insectivore-Omnivore", "Herbivore-Omnivore", "Herbivore-Insectivore", "Overall")' )
args = c('r=CategoricalPrunedCarnivoreTree', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("Carnivore-Herbivore", "Carnivore-Omnivore", "Herbivore-Omnivore", "Overall")' )
args = c('r=CategoricalBinaryInsectivoreTree', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("Background-Insectivore", "Overall")' )
args = c('r=CategoricalBinaryVertivoreTree', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("Background-Vertivore", "Overall")' )
args = c('r=CategoricalBinaryHerbivoreTree', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("Background-Herbivore", "Overall")' )
args = c('r=CategoricalBinaryOmnivoreTree', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("Background-Omnivore", "Overall")' )


args = c('r=CategoricalInsVertivoreTree', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("Omnivore-Vertivore", "Insectivore-Vertivore", "Herbivore-Vertivore", "Insectivore-Omnivore", "Herbivore-Omnivore", "Herbivore-Insectivore", "Overall")' )

args = c('r=CategoricalInsVertivoreTree', 'm=c("Data/KeggReactome.gmt")', 'p=F', 's=c("Omnivore-Vertivore", "Insectivore-Vertivore", "Herbivore-Vertivore", "Insectivore-Omnivore", "Herbivore-Omnivore", "Herbivore-Insectivore", "Overall")' )
args = c('r=CategoricalPrunedCarnivoreTree', 'm=c("Data/KeggReactome.gmt")', 'p=F', 's=c("Carnivore-Herbivore", "Carnivore-Omnivore", "Herbivore-Omnivore", "Overall")' )

args = c('r=meanTemp', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt", "Data/KeggReactome.gmt")', 'p=F')


args = c('r=CategoricalInsVertivoreTreeLiamInference', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("Omnivore-Vertivore", "Insectivore-Vertivore", "Herbivore-Vertivore", "Insectivore-Omnivore", "Herbivore-Omnivore", "Herbivore-Insectivore", "Overall")' )
args = c('r=CategoricalInsVertivoreTreeLiamInference', 'm=c("Data/KeggReactome.gmt")', 'p=F', 's=c("Omnivore-Vertivore", "Insectivore-Vertivore", "Herbivore-Vertivore", "Insectivore-Omnivore", "Herbivore-Omnivore", "Herbivore-Insectivore", "Overall")' )

args = c('r=CategoricalInsVertivoreTreeCarnivoreLiamInference', 'm=c("Data/KeggReactome.gmt")', 'p=F', 's=c("Carnivore-Omnivore", "Carnivore-Herbivore", "Herbivore-Omnivore", "Overall")' )
args = c('r=CategoricalInsVertivoreTreeCarnivoreLiamInference', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("Carnivore-Omnivore", "Carnivore-Herbivore", "Herbivore-Omnivore", "Overall")' )

args = c('r=Demo', 'v=F', 'm=c("Data/KeggReactome.gmt", "Data/DisGeNET.gmt")', 's=c("Carnivore-Herbivore", "Carnivore-Omnivore", "Herbivore-Omnivore", "Overall")')


args = c('r=CategoricalInsVertivoreTreeLiamInference', 'm=c("Data/KeggReactome.gmt", "Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=C', 's=c("Insectivore-Vertivore", "Herbivore-Vertivore", "Insectivore-Omnivore", "Herbivore-Omnivore", "Herbivore-Insectivore", "Overall")' )

args = c('r=CategoricalInsVertivoreTreeDuplicateLiamInference', 'm=c("Data/KeggReactome.gmt")', 'p=F', 's=c("Insectivore-Vertivore", "Herbivore-Vertivore", "Insectivore-Omnivore", "Herbivore-Omnivore", "Herbivore-Insectivore", "Overall")' )

args = c('r=ComplexDietCentralAnalysis', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("Insectivore-Vertivore", "Herbivore-Vertivore", "Insectivore-Omnivore", "Herbivore-Omnivore", "Herbivore-Insectivore", "Carnivore-Herbivore", "Overall")' )


args = c('r=ComplexDietCentralAnalysis', 'm=c("Data/KeggReactome.gmt")', 'p=C', 's=c("Insectivore-Vertivore", "Herbivore-Vertivore", "Insectivore-Omnivore", "Herbivore-Omnivore", "Herbivore-Insectivore", "Overall")' )
args = c('r=ComplexDietCentralAnalysisSimplify', 'm=c("Data/KeggReactome.gmt")', 'p=F', 's=c("Invertivore-Vertivore", "Herbivore-Vertivore", "Invertivore-Omnivore", "Herbivore-Omnivore", "Herbivore-Invertivore", "Overall", "Carnivore-Herbivore", "Carnivore-Omnivore")' )

args = c('r=ComplexDietCentralAnalysisSimplifyNoScav', 'm=c("Data/KeggReactome.gmt")', 'p=F', 's=c("Invertivore-Vertivore", "Herbivore-Vertivore", "Invertivore-Omnivore", "Herbivore-Omnivore", "Herbivore-Invertivore", "Overall", "Carnivore-Herbivore", "Carnivore-Omnivore")' )
args = c('r=ComplexDietCentralAnalysisPiscFix', 'm=c("Data/KeggReactome.gmt")', 'p=F', 's=c("Invertivore-Vertivore", "Herbivore-Vertivore", "Invertivore-Omnivore", "Herbivore-Omnivore", "Herbivore-Invertivore", "Overall", "Carnivore-Herbivore", "Carnivore-Omnivore")' )
args = c('r=ComplexDietCentralAnalysisNoFishbat', 'm=c("Data/KeggReactome.gmt")', 'p=F', 's=c("Invertivore-Vertivore", "Herbivore-Vertivore", "Invertivore-Omnivore", "Herbivore-Omnivore", "Herbivore-Invertivore", "Overall", "Carnivore-Herbivore", "Carnivore-Omnivore")' )
args = c('r=ComplexDietCentralAnalysisSimplifyEqualDrop', 'm=Data/KeggReactome.gmt', 'p=F') 

args = c('r=ComplexDietCentralAnalysisSimplify2', 'm=c("Data/MGI_Mammalian_Phenotype_Level_4.gmt", "Data/GO_Biological_Process_2023.gmt", "Data/DisGeNET.gmt", "Data/tissue_specific.gmt", "Data/EnrichmentHsSymbolsFile2.gmt")', 'p=F', 's=c("Invertivore-Vertivore", "Herbivore-Vertivore", "Invertivore-Omnivore", "Herbivore-Omnivore", "Herbivore-Invertivore", "Carnivore-Herbivore", "Overall")' )


args = c('r=ComplexDietCentralAnalysisSimplify2', 'm=c("Data/KeggReactome.gmt")', 'p=C', 's=c("Invertivore-Vertivore", "Herbivore-Vertivore", "Herbivore-Invertivore")' )











#---------------------------------------------------------------------
# --- RunRERAndCorrelation.R --- 
# --------------------------------------------------------------------




args = c('r=CVO', 'm=data/RemadeTreesAllZoonomiaSpecies.rds', 'v=F', 's=b') #This is a debug argument set. It is used to set arguments locally, when not running the code through a bash script.
args = c("r=EcholocationUpdate2", "m=data/RemadeTreesAllZoonomiaSpecies.rds", "s=b")
args = c("m=Data/NoSignExpressionTreesRound3.rds", "r=LiverExpression2", "v=F", "s=b")
args = c("m=data/RemadeTreesAllZoonomiaSpecies.rds", "r=CategoricalDiet4Phen", "v=T", "s=c")
args = c("m=data/RemadeTreesAllZoonomiaSpecies.rds", "r=CVHRemake", "s=b")
args = c("m=data/RemadeTreesAllZoonomiaSpecies.rds", "r=BinaryCVHApplesToApples", "s=b")
args = c('r=RubyRegenARD',   'm=data/mam120aa_trees.rds', 'v=F', 't=ARD')
args = c('r=RubyRegenER',   'm=data/mam120aa_trees.rds', 'v=F', 't=ER')
args = c('r=ThreePhenLikeihoodTest', 'm=data/mam120aa_trees.rds', 'v=F', 's=g')
args = c('r=HMGRelaxTest', 'm=data/mam120aa_trees.rds', 'v=F', 's=g')
args = c('r=IPCRelaxTest', 'm=data/mam120aa_trees.rds', 'v=F', 's=g')

args = c('r=NewHiller4Phen', 'm=data/newHillerMainTrees.rds', 'v=F', 's=g')
args = c('r=NewHiller2Phen', 'm=data/newHillerMainTrees.rds', 'v=T', 's=g')
args = c('r=NewHillerTestSupraPrimates', 'm=data/newHillerMainTrees.rds', 'v=T', 's=g')

args = c('r=MaturityLifespanPercent', 'm=data/newHillerMainTrees.rds', 's=c','v=F')
args = c('r=MaturityLogRaw', 'm=data/newHillerMainTrees.rds', 's=c','v=F')
args = c('r=PankajBodysize', 'm=data/newHillerMainTrees.rds', 's=c','v=F')
args = c('r=CVHNew', 's=b', 'v=T', 'm=data/RemadeTreesAllZoonomiaSpecies.rds')
args = c('r=CIvAllZoonomia', 's=b', 'v=T', 'm=data/RemadeTreesAllZoonomiaSpecies.rds')

args = c("r=CVHNew", 'm=data/RemadeTreesAllZoonomiaSpecies.rds', "s=b")
args = c("r=CategoricalMobivoreTree", 'm=data/zoonomiaAllMammalsTrees.rds', "s=g", "v=F")
args = c("r=RyanPermualtionsDemo", 'm=data/zoonomiaAllMammalsTrees.rds', "s=g", "v=F")

args = c("r=CategoricalCarnivoreTree", 'm=data/zoonomiaAllMammalsTrees.rds', "s=g", "v=F")
args = c("r=CategoricalInsVertivoreTree", 'm=data/zoonomiaAllMammalsTrees.rds', "s=g", "v=F")
args = c("r=CategoricalInsVertivoreTree", 'm=data/zoonomiaAllMammalsTrees.rds', "s=g", "v=F", "l=170")
args = c("r=CategoricalPrunedCarnivoreTree", 'm=data/zoonomiaAllMammalsTrees.rds', "s=g", "v=F", "l=170")
args = c("r=CategoricalBinaryInsectivoreTree", 'm=data/zoonomiaAllMammalsTrees.rds', "s=g", "v=F", "l=50")
args = c("r=CategoricalBinaryVertivoreTree", 'm=data/zoonomiaAllMammalsTrees.rds', "s=g", "v=F", "l=170")
args = c("r=CategoricalBinaryHerbivoreTree", 'm=data/zoonomiaAllMammalsTrees.rds', "s=g", "v=F", "l=50")
args = c("r=CategoricalBinaryOmnivoreTree", 'm=data/zoonomiaAllMammalsTrees.rds', "s=g", "v=F", "l=170")
args = c("r=CategoricalBinaryCarnivoreTree", 'm=data/zoonomiaAllMammalsTrees.rds', "s=g", "v=F", "l=170")


args = c("r=CategoricalNoMegabranchInsvertTree", 'm=data/zoonomiaAllMammalsTrees.rds', "s=g", "v=F", "l=170")
args = c("r=CategoricalDownsampledInsvertTree", 'm=data/zoonomiaAllMammalsTrees.rds', "s=g", "v=F", "l=170")

args = c('r=meanTemp', 'm=data/zoonomiaAllMammalsTrees.rds', 's=c','v=T', 'n=ZoonomiaTip')

args = c("r=CategoricalSlimMainInsVertivoreTree", 'm=data/CategoricalInsVertivoreMaintrees.rds', "s=g", "v=F", "l=170")
args = c("r=CategoricalInsVertivoreTreeLiamInference", 'm=data/zoonomiaAllMammalsTrees.rds', "s=g", "v=F", "l=170")
args = c("r=CategoricalInsVertivoreTreeCarnivoreLiamInference", 'm=', "s=g", "v=F", "l=170")

args = c('r=Demo', 'v=F', 'm=demoTreeFile', 's=g')

args = c('r=Demo', 'v=T', 'm=demoTreeFile', 's=g')

args = c("r=PredatorFullTree", 'm=', "s=g", "v=F", "l=170")

args = c("r=ComplexDietCentralAnalysisEqualLengthInference", 'm=data/zoonomiaAllMammalsTrees.rds', "s=g", "v=F", "l=170")

args = c("r=ComplexDietBinaryInvertivore", 'm=data/zoonomiaAllMammalsTrees.rds', "s=g", "v=F", "l=170")
args = c("r=ComplexDietBinaryVertivore", 'm=data/zoonomiaAllMammalsTrees.rds', "s=g", "v=F", "l=170")
args = c("r=ComplexDietBinaryHerbivore", 'm=data/zoonomiaAllMammalsTrees.rds', "s=g", "v=F", "l=170")
args = c("r=ComplexDietBinaryOmnivore", 'm=data/zoonomiaAllMammalsTrees.rds', "s=g", "v=F", "l=170")
args = c("r=ComplexDietBinaryCarnivore", 'm=data/zoonomiaAllMammalsTrees.rds', "s=g", "v=F", "l=170")

args = c("r=ComplexDietCentralAnalysisSimplifyAltReconstruction", 'm=data/zoonomiaAllMammalsTrees.rds', "s=g", "v=F", "l=170")




#---------------------------------------------------------------------
# --- RumMergedAnalysis.R --- 
# --------------------------------------------------------------------



args = c('r=CategoricalInsVertivoreTreeDuplicateLiamInference', 'm=data/zoonomiaAllMammalsTrees.rds')
args = c('r=CategoricalInsVertivoreTreeNoYeastLiamInference', 'm=data/zoonomiaAllMammalsTrees.rds')

args = c('r=ComplexDietCentralAnalysis', 'm=data/zoonomiaAllMammalsTrees.rds')
args = c('r=ComplexDietCentralAnalysisSimplify2', 'm=data/zoonomiaAllMammalsTrees.rds')
args = c('r=ComplexDietCentralAnalysisSimplifyAllSpecies', 'm=data/zoonomiaAllMammalsTrees.rds')

args = c('r=ComplexDietCentralAnalysisSimplifyAltReconstruction', 'm=data/zoonomiaAllMammalsTrees.rds')
