clusterRun = F
clusterRun = T

if(clusterRun){.libPaths("/share/ceph/wym219group/shared/libraries/R4")} #add path to custom libraries to searched locations
library(GSEABase)
source("Src/Reu/cmdArgImport.R")

args = c("r=CategoricalInsVertivoreTree", "s=Carnivore-Herbivore", 
         "a=CategoricalBinaryCarnivoreTree", "b=Carnivore", "c=CategoricalBinaryHerbivoreTree", "d=Herbivore",
         "o=T", 'l=c("HV", "HI", "CH")', "k=H", "i=T")
args = c("r=CategoricalInsVertivoreTree", "s=Herbivore-Vertivore", 
         "a=CategoricalBinaryHerbivoreTree", "b=Herbivore", "c=CategoricalBinaryVertivoreTree", "d=Vertivore",
         "o=T", 'l=c("HV", "HI", "CH")', "k=H", "i=T")
args = c("r=CategoricalInsVertivoreTree", "s=Herbivore-Insectivore", 
         "a=CategoricalBinaryHerbivoreTree", "b=Herbivore", "c=CategoricalBinaryInsectivoreTree", "d=Insectivore",
         "o=T", 'l=c("HV", "HI", "CH")', "k=H", "i=T")

args = c("r=CategoricalInsVertivoreTree", "s=Insectivore-Vertivore", 
         "a=CategoricalBinaryInsectivoreTree", "b=Insectivore", "c=CategoricalBinaryVertivoreTree", "d=Vertivore",
         "o=T", 'l=c("HV", "HI", "CH")', "k=H", "i=T")



args = c("r=CategoricalInsVertivoreTreeLiamInference", "s=Carnivore-Herbivore", 
         "a=CategoricalBinaryCarnivoreTree", "b=Carnivore", "c=CategoricalBinaryHerbivoreTree", "d=Herbivore",
         "o=T", 'l=c("HV", "HI", "CH")', "k=H", "i=T")
args = c("r=CategoricalInsVertivoreTreeLiamInference", "s=Herbivore-Vertivore", 
         "a=CategoricalBinaryHerbivoreTree", "b=Herbivore", "c=CategoricalBinaryVertivoreTree", "d=Vertivore",
         "o=T", 'l=c("HV", "HI", "CH")', "k=H", "i=T")
args = c("r=CategoricalInsVertivoreTreeLiamInference", "s=Herbivore-Insectivore", 
         "a=CategoricalBinaryHerbivoreTree", "b=Herbivore", "c=CategoricalBinaryInsectivoreTree", "d=Insectivore",
         "o=T", 'l=c("HV", "HI", "CH")', "k=H", "i=T")

args = c("r=CategoricalInsVertivoreTreeLiamInference", "s=Insectivore-Vertivore", 
         "a=CategoricalBinaryInsectivoreTree", "b=Insectivore", "c=CategoricalBinaryVertivoreTree", "d=Vertivore",
         "o=T", 'l=c("HV", "HI", "CH")', "k=H", "i=T")

{
  # --- Standard start-up code ---
  if(clusterRun)args = commandArgs(trailingOnly = TRUE)
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
  
  # --- Import arguments --- 
  subdirectory = "Carnivore-Herbivore"
  BinaryTreeOne = "CategoricalBinaryHerbivoreTree"
  BinaryPhenotype = "Herbivore"
  BinaryTreeTwo = "CategoricalBinaryCarnivoreTree"
  BinaryPhenotypeTwo = "Carnivore"
  GOGroup = "KeggReactome"
  
  useOverlap = T
  overlapValues = c("HV", "HI", "CH")
  overlapBackground = ("H")
  pvalueCutoff = 0.05
  useDriver = T
  
  {
  #subdirectory
  if(!is.na(cmdArgImport('s'))){
    subdirectory = cmdArgImport('s')
  }else{
    stop("Missing Arugment")
  }
  
  #Binary Location one 
  if(!is.na(cmdArgImport('a'))){
    BinaryTreeOne = cmdArgImport('a')
  }else{
    stop("Missing Arugment")
  }
  
  #Binary Prefix one
  if(!is.na(cmdArgImport('b'))){
    BinaryPhenotype = cmdArgImport('b')
  }else{
    stop("Missing Arugment")
  }
  
  #Binary Location two 
  if(!is.na(cmdArgImport('c'))){
    BinaryTreeTwo = cmdArgImport('c')
  }else{
    stop("Missing Arugment")
  }
  
  #Binary Prefix two 
  if(!is.na(cmdArgImport('d'))){
    BinaryPhenotypeTwo = cmdArgImport('d')
  }else{
    stop("Missing Arugment")
  }
  
  #GOGroup 
  if(!is.na(cmdArgImport('g'))){
    GOGroup = cmdArgImport('g')
  }else{
    message("No GO arg, using KeggReactome")
  }
  
  #Use Overlap 
  if(!is.na(cmdArgImport('o'))){
    useOverlap = cmdArgImport('o')
  }else{
    message("No overlap Arg, using overlap")
  }
  
  #Overlap list
  if(!all(is.na(cmdArgImport('l')))){
    overlapValues = cmdArgImport('l')
  }else{
    if(useOverlap){
      stop("Missing Arugment")
    }
  }
  
  #Overlap Background
  if(!is.na(cmdArgImport('k'))){
    overlapBackground = cmdArgImport('k')
  }else{
    if(useOverlap){
      stop("Missing Arugment")
    }
  }
  
  #p value cutoff 
  if(!is.na(cmdArgImport('p'))){
    pvalueCutoff = cmdArgImport('p')
  }else{
    if(useOverlap){
      message("No pvalue cutoff, using 0.05")
    }
  }
  
  #use driver
    #Use Overlap 
    if(!is.na(cmdArgImport('i'))){
      useDriver = cmdArgImport('i')
    }else{
      message("No driver Arg, using driver")
    }
  }
}



# -- setup output files ---
cytoscapeDirectory = paste0("Output/", filePrefix, "/", subdirectory, "/", "Cytoscape")
if(!dir.exists(cytoscapeDirectory)){                       #create that directory if it does not exist
  dir.create(cytoscapeDirectory)
}

gmtFilename = paste0("Data/", GOGroup, ".gmt")
file.copy(gmtFilename, paste0(cytoscapeDirectory, "/", GOGroup, ".gmt")) #make a copy in the cytoscape direcotry for easy cytoscape work 




# ------ Driver analysis ---- 
if(useDriver){
  #------ Run driver assessment code if not already run  --- 
  driverTableFilename = paste0("Output/", filePrefix, "/", subdirectory, "/", filePrefix, subdirectory, "DriverTable")
  if(!file.exists(paste0(driverTableFilename, ".rds"))){
    source("Src/Reu/AssessRERDriver.R")
    driverTable = AssessRERDriver(filePrefix,subdirectory , BinaryTreeOne, BinaryPhenotype,  BinaryTreeTwo, BinaryPhenotypeTwo)
    write.csv(driverTable, paste0(driverTableFilename, ".csv"))
    saveRDS(driverTable, paste0(driverTableFilename, ".rds"))
  }else{
    driverTable = readRDS(paste0(driverTableFilename, ".rds"))
  }
  
  goDriverTableFilename = paste0("Output/", filePrefix, "/", subdirectory, "/", filePrefix, subdirectory, "GoDriverTable")
  if(!file.exists(paste0(goDriverTableFilename, ".csv"))){
    source("Src/Reu/AssessGOCategoryDriver.R")
    GODriver = AssessGoCategoryDriver(filePrefix, subdirectory, GOGroup, 0.1, 1, F)
    View(GODriver)
    write.csv(GODriver, paste0(goDriverTableFilename, ".csv"))
  }
  
  #--- split GO Driver into positive and negative
  
  goDriverTableFilename = paste0("Output/", filePrefix, "/", subdirectory, "/", filePrefix, subdirectory, "GoDriverTable")
  GODriver = read.csv(paste0(goDriverTableFilename, ".csv"))
  
  
  GoDriverPositive = GODriver[which(GODriver$stat > 0),]
  GoDriverNegative = GODriver[which(GODriver$stat < 0),]
  
  GODriverPositiveColored = GODriver
  GODriverPositiveColored$p.adj[GODriverPositiveColored$stat < 0] = 0.11
  GODriverPositiveColored$pval[GODriverPositiveColored$stat < 0] = 0.11
  
  

  
  goPositiveGoDriverTableFilename = paste0("Output/", filePrefix, "/", subdirectory, "/", "Cytoscape/", filePrefix, subdirectory, "GoDriverPositiveTable")
  goNegativeGoDriverTableFilename = paste0("Output/", filePrefix, "/", subdirectory, "/", "Cytoscape/", filePrefix, subdirectory, "GoDriverNegativeTable")
  goColoredGoDriverTableFilename = paste0("Output/", filePrefix, "/", subdirectory, "/", "Cytoscape/", filePrefix, subdirectory, "GoDriverPositiveColoredTable")
  
  write.table(GODriverPositiveColored, paste0(goColoredGoDriverTableFilename, ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.csv(GoDriverPositive, paste0(goPositiveGoDriverTableFilename, ".csv"), row.names = F)
  write.csv(GoDriverNegative, paste0(goNegativeGoDriverTableFilename, ".csv"), row.names = F)
  write.csv(GODriverPositiveColored, paste0(goColoredGoDriverTableFilename, ".csv"), row.names = F)
  
  GODriver
  
  
  GOCytoscape = GODriver[,c(1,7,3,4,2,6)]
  colnames(GOCytoscape) = c("GO.ID", "Description", "pVal", "p.adj", "Phenotype", "Gene.vals")
  GOCytoscape$Phenotype = sign(GOCytoscape$Phenotype)
  #GOCytoscape$Phenotype[GOCytoscape$Phenotype == 1] = "1"
  write.table(GOCytoscape, paste0(cytoscapeDirectory, "/CytoscapeInput.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  GOCytoscapePositive = GOCytoscape[which(GOCytoscape$Phenotype >0),]
  GOCytoscapePositiveFilename = paste0(cytoscapeDirectory, "/Cytoscape", BinaryPhenotype,"FasterInput.txt")
  write.table(GOCytoscapePositive, GOCytoscapePositiveFilename, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  GOCytoscapeNegative = GOCytoscape[which(GOCytoscape$Phenotype <0),]
  GOCytoscapeNegativeFilename = paste0(cytoscapeDirectory, "/Cytoscape", BinaryPhenotypeTwo,"FasterInput.txt")
  write.table(GOCytoscapeNegative, GOCytoscapeNegativeFilename, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  #Use that cytoscape input file as the data in Cytoscape, on the generic setting
  
  # -- save as gmt file --- 
  gmtData = getGmt(gmtFilename)
  
  gmtNames = names(gmtData)
  gmtDriver = GODriver$Driver[match(gmtNames, GODriver$X)]
  gmtUpdate = data.frame(gmtNames, gmtDriver)
  write.csv(gmtUpdate, paste0(cytoscapeDirectory, "/GmtDirectionColumn.csv"), row.names = F)
  #This creates a column with the driving phenotype information. 
  #Open the gmt file in excel, and replace the description column with the produced driver column. 
  #This is set up this way because R is bad at handling variable row lengths, so it is easier to edit the gmt file in excel.
  
}


# -- Make overlap gmt column instead ---
if(useOverlap){
  combinedGOFilename = paste0("Output/", filePrefix, "/", filePrefix, "combinedGOResults-", GOGroup, ".rds")
  combinedGOData = readRDS(combinedGOFilename)
  
  combinedGOData$overlapName = NA
  
  testColumns = paste0(overlapValues, "-p.adj")
  for(i in 1:nrow(combinedGOData)){
    currentOverlapName = ""
    for(j in testColumns){
      currentTest = which(colnames(combinedGOData) == j)
      currentValue = combinedGOData[i, currentTest]
      if(currentValue < pvalueCutoff){
        overlapNameAddition = j
        overlapNameAddition = gsub("-p.adj", "", overlapNameAddition)
        overlapNameAddition = gsub(overlapBackground, "", overlapNameAddition)
        currentOverlapName = paste0(currentOverlapName, overlapNameAddition)
      }
    }
    combinedGOData$overlapName[i] = currentOverlapName
  }
 
  
  gmtData = getGmt(gmtFilename)
  gmtNames = names(gmtData)
  gmtOverlap = combinedGOData$overlapName[match(gmtNames, rownames(combinedGOData))]
  gmtUpdate = data.frame(gmtNames, gmtOverlap)
  write.csv(gmtUpdate, paste0(cytoscapeDirectory, "/GmtOverlapColumn.csv"), row.names = F)
   
}


# --- Guide on how to convert these outputs into a cytoscape figure: ----

#Run the NewCytoscapeMaker script
#Copy the Driver column over to the local gmt file 
#Make sure to remove the first row holding the column names before copying
#Make a new Enrichment map (At either 0.1 cutoff or 1 cutoff, see below)
#Input the gmt file as a shared file
#Input cytoscape input at dataset 1
#Set the style on the enrichment map page to chartdata = None
#Make style changes
  #Change border paint mapping
    #Set mapping to fdr_qvalue
    #Set mapping type to continuous
    #Set min to 0 and max to 0.1
    #Set 0.1 side to white
      #From the black-to-white column, not the orange column (the orange column top box is not pure white)
    #Set 0.05 to mid-orange (Two down from top)
    #Set 0 to dark orange
  #Change border width to 7.0
  #Change Fill color settings 
    #Set to be based on Discrete mapping
    #Set to be based on description
    #Use colors based on phenotypes 
  #Change label to be based on Name
  #Set shape to be based on driver 
    #Set to continuous mapping
    #Set to be based on colouring (stand in for phenotype) 
    #Set min to -1
    #Set max to 1
    #Add node (this defaults to zero, which is correct) 
    #Set the shape on one side to different from the other side
  #Use auto-annotate to generate clusters
    #Using a p value cutoff of 1 (all sets) in the original enrichmentMap run is useful for positioning pathways reasonably. 
    #However, the annotations produced are much less helpful.
  #Manually rename autogenerated clusters based on cluster contents 
          