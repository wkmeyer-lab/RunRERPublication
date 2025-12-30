library(ggtree)
library(ggimage)
library(rphylopic)
library(RERconverge)
library(utils)
library(tidytree)
library(ape)
source("Src/Reu/ZoonomTreeNameToCommon.R")
source("Src/Reu/cmdArgImport.R")
nodeid.tbl_tree <- utils::getFromNamespace("nodeid.tbl_tree", "tidytree")
rootnode.tbl_tree <- utils::getFromNamespace("rootnode.tbl_tree", "tidytree")
offspring.tbl_tree <- utils::getFromNamespace("offspring.tbl_tree", "tidytree")
offspring.tbl_tree_item <- utils::getFromNamespace(".offspring.tbl_tree_item", "tidytree")
child.tbl_tree <- utils::getFromNamespace("child.tbl_tree", "tidytree")
parent.tbl_tree <- utils::getFromNamespace("parent.tbl_tree", "tidytree")


args =c("r=CategoricalInsVertivoreTree", 'p=c("darkgreen", "darkblue", "black", "red")', 'c=c("Herbivore", "Invertivore", "Omnivore", "Vertivore")', 'n=ZoonomiaTip', "l=Diet", "i=F" )
args =c("r=makeLalithaTree", 'p=c("darkgreen", "darkblue", "black")', 'c=c("1", "2", "3")', 'n=ZoonomiaTip', "l=Phen", "i= T")

args =c("r=CategoricalInsVertivoreTreeLiamInference", 'p=c("darkgreen", "darkblue", "black", "red")', 'c=c("Herbivore", "Invertivore", "Omnivore", "Vertivore")', 'n=ZoonomiaTip', "l=Diet", "i=F" )



# -- Standard Startup code -- 
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

# -- Argument imports --- 
{
  mainTreesLocation = "Data/zoonomiaAllMammalsTrees.rds"
  palette(c( "darkgreen", "darkblue", "black", "red"))
  spreadSheetLocation = "Data/mergedData.csv"
  nameColumn = "tipName"
  CategoryReplacements = NULL
  legendText = NA
  imageAllTips = T
  
  #MainTrees Location
  if(!is.na(cmdArgImport('m'))){
    mainTreesLocation = cmdArgImport('m')
  }else{
    message("No maintrees arg, using Data/zoonomiaAllMammalsTrees.rds")
  }
  
  #Pallette
  if(!all(is.na(cmdArgImport('p')))){
    paletteValues = cmdArgImport('p')
    palette(paletteValues)
  }else{
    message("No Palette Provided, using: darkgreen, darkblue, black, red")
  }
  
  #spreadsheet File
  if(!is.na(cmdArgImport('d'))){
    spreadSheetLocation = cmdArgImport('d')
  }else{
    message("Using Data/mergedData.csv spreadsheet")
  }
  
  #Name Column
  if(!is.na(cmdArgImport('n'))){
    nameColumn = cmdArgImport('n')
  }else{
    message("Name Column not specified, using 'tipName'.")
  }
  
  #Category replacements 
  if(!all(is.na(cmdArgImport('c')))){
    CategoryReplacements = cmdArgImport('c')
  }else{
    message("No category replacements provided, using category branch lengths as labels")
  }
  
  #legend Label
  if(!is.na(cmdArgImport('l'))){
    legendText = cmdArgImport('l')
  }else{
    message("No legendText provided, legend will be unlabeled")
  }
  
  #imageAllTips
  if(!is.na(cmdArgImport('i'))){
    imageAllTips = as.logical(cmdArgImport('i'))
  }else{
    message("Use of clade labels not specified, labeling all tips")
  }
}





if(!exists("mainTrees")){mainTrees = readRDS(mainTreesLocation)}
#commonMainTrees = mainTrees
#commonMainTrees$masterTree = ZoonomTreeNameToCommon(commonMainTrees$masterTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn, plot = F)
#commonMasterTree = commonMainTrees$masterTree
commonMasterTree = ZoonomTreeNameToCommon(mainTrees$masterTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn, plot = F)

categoricalCommonTreeFilename = paste(outputFolderName, filePrefix, "CategoricalCommonTree.rds", sep="") #make a filename based on the prefix
commonCategoricalTree = readRDS(categoricalCommonTreeFilename)
#scientificCategoricalTreeFilename = paste(outputFolderName, filePrefix, "CategoricalScientificTree.rds", sep="") #make a filename based on the prefix
#scientificCategoricalTree = readRDS(scientificCategoricalTreeFilename)
scientificCategoricalTree = ZoonomTreeNameToCommon(commonCategoricalTree, manualAnnotLocation = spreadSheetLocation, tipCol = "CommonName", scientific = T, scientificCol = "Scientific_Binomial", plot = F)

commonMasterTrimmed = drop.tip(commonMasterTree, commonMasterTree$tip.label[!commonMasterTree$tip.label %in% commonCategoricalTree$tip.label])
if(!is.rooted(commonCategoricalTree)){commonMasterTrimmed = unroot(commonMasterTrimmed)}
scientificMasterTrimmed = ZoonomTreeNameToCommon(commonMasterTrimmed, manualAnnotLocation = spreadSheetLocation, tipCol = "CommonName", scientific = T, scientificCol = "Scientific_Binomial", plot = F)
if(!is.rooted(scientificCategoricalTree)){scientificMasterTrimmed = unroot(scientificMasterTrimmed)}

commonCategoricalTreeEdgeLengths = commonCategoricalTree$edge.length
commonCategoricalTreeEdgeLengths = as.character(commonCategoricalTreeEdgeLengths)
edge=data.frame(commonCategoricalTree$edge, edge_num=1:nrow(commonCategoricalTree$edge))
colnames(edge)=c("parent", "node", "edge_num")
edge$Categorylength = commonCategoricalTree$edge.length
edge$CategorylengthChar = as.character(edge$Categorylength)
if(!is.null(CategoryReplacements)){
  for(i in 1:length(unique(edge$CategorylengthChar))){
    edge$CategorylengthChar[edge$CategorylengthChar == i] = CategoryReplacements[i]
  }
}

commonCategoricalTree$edge.length = commonMasterTrimmed$edge.length
scientificCategoricalTree$edge.length = scientificMasterTrimmed$edge.length




phylopicNames = NULL
uuidList = NULL
missingPictures = NA

uuidListFilename =  paste(outputFolderName, filePrefix, "UuidList.rds", sep="") #make a filename based on the prefix
if(!file.exists(uuidListFilename) | forceUpdate){                             #if it does not exist, or update is forced 
  inTips = scientificCategoricalTree$tip.label
  for(i in 1:length(inTips)){
    print(i)
    phylopicNames[i] = tryCatch({autocomplete_name(inTips[i])[1,2]}, error = function(msg) {return(inTips[i])})
    
    uuidList[i] = tryCatch(
      {get_uuid(phylopicNames[i])}, 
      error = function(msg){
        if(length(grep(" ", phylopicNames[i])) > 0){genusName = strsplit(phylopicNames[i], " ")[[1]][1]}else{
          if(length(grep("_", phylopicNames[i])) > 0){genusName = strsplit(phylopicNames[i], "_")[[1]][1]}
        }
        tryCatch(
          {get_uuid(genusName)},
          error = function(msg){
            message(paste("No pic found for ", phylopicNames[i]))
            missingReport = phylopicNames[i]
            names(missingReport) = i 
            append(missingPictures, missingReport)
            return("NULL")
          }
        )
      }
    )
  }
  saveRDS(uuidList, uuidListFilename)
  }else{                                                                          #Otherwise
    uuidList = readRDS(uuidListFilename)                                              #Use the existing ones
}


# Replace any missing UUIDs with tardigrades as debug images 

#uuidList[21] = get_uuid("tardigrades")

#uuidList[uuidList == "NULL"] = NULL

tip_data = data.frame(
  scientificlabel = scientificCategoricalTree$tip.label,
  uuid = uuidList,
  node = 1:length(scientificCategoricalTree$tip.label),
  stringsAsFactors = FALSE
)
#tip_data$uuid[tip_data$uuid == "NULL"] = NULL


ggTreeOut = ggtree(commonCategoricalTree, layout = "circular") +scale_color_manual(values=palette()) 
ggTreeOut = ggTreeOut %<+% edge + aes(color=CategorylengthChar)
ggTreeOut = ggTreeOut %<+% tip_data 
#ggTreeOut$data$label = paste(ggTreeOut$data$label, "-", ggTreeOut$data$node, sep="")
#ggTreeOut = ggTreeOut + geom_tiplab()
#ggTreeOut = ggTreeOut + geom_tiplab(geom = "phylopic", aes(image = uuid))
#ggTreeOut + geom_phylopic(aes(uuid = uuid), color = "black", alpha = 1, size = 0.08)
ggTreeOut

#make the ggtree object temporarily for refrencing in later code
ggTreeClades = ggtree(commonCategoricalTree, layout = "circular") +scale_color_manual(values=palette()) 


if(imageAllTips){
  
  
  missingPhylopics = which(uuidList == "NULL")
  for(i in 1:length(missingPhylopics)){
    
    message( "------- WARNING -------")
    cat( "------- WARNING ------- \n")
    cat( "Missing phylopic for species:", commonCategoricalTree$tip.label[missingPhylopics][i], "\n")
    cat( "This image has been replaced with a tardigrade to stop the code from breaking. \n Either find an appropriate phylopic and manually replace it, \n or remove the tardigrade manually from the final output. \n")
    cat( "Your total number of tardigrades is:", length(missingPhylopics))
  }
  
  
  clades = data.frame(
    node = 1:length(commonCategoricalTree$tip.label),
    cladelabel = commonCategoricalTree$tip.label,
    uuid = tip_data$uuid,
    phylopic = tip_data$uuid
    
  )
  
  clades$x_new = NA
  clades$y_new = NA
  baseOffset = 0.005
  for(i in 1:nrow(clades)){
    cladeNode = clades$node[i]
    
    #message(clades$cladelabel[i])
    #message(cladeNode)
    
    coords = ggTreeClades$data[ggTreeClades$data$node == cladeNode, ]
    
    #Determine the longest branch length from the MRCA node, based on the longest tip species 
    tipDescendants = Descendants(commonCategoricalTree, cladeNode, type = "tip")[[1]]
    totalDistances = NULL
    
    for(j in 1:length(tipDescendants)){
      totalLength = 0 
      currentAncestor = tipDescendants[j]
      #message("Child path:", j)
      while(currentAncestor != cladeNode){
        #message(paste("processing Node:", currentAncestor))
        currentLength = commonCategoricalTree$edge.length[which(commonCategoricalTree$edge[, 2] == currentAncestor)]
        #message(paste("Length:", currentLength))
        if(currentLength > totalLength/100){totalLength = totalLength + currentLength}; message("adding")
        #message(paste("Total Length:", totalLength))
        currentAncestor = Ancestors(commonCategoricalTree, currentAncestor, "parent")
        if(length(currentAncestor) == 0){break}
      }
      #message(paste("Total Length:", totalLength))
      #message("-")
      totalDistances = append(totalDistances, totalLength)
      
    }
    maxLength = max(totalDistances)
    
    #message(paste("Max Length:", maxLength))
    x_new = coords$x + maxLength + baseOffset
    y_new = coords$y
    clades$x_new[i] = x_new
    clades$y_new[i] = y_new
    #message(x_new)
    #message(y_new)
    #message("------")
    
  }
}else{
  {
    collapsedClades = data.frame()
    collapsedClades[1,] = NA
    
    collapsedClades$Platypus = MRCA(commonCategoricalTree, 195)
    collapsedClades$Opossums = MRCA(commonCategoricalTree, c(194,192,193))
    collapsedClades$Koala = MRCA(commonCategoricalTree, c(190,191))
    collapsedClades$Kangaroos = MRCA(commonCategoricalTree, c(186,187,189,188))
    collapsedClades$Anteaters = MRCA(commonCategoricalTree, c(183,182))
    collapsedClades$Sloths = MRCA(commonCategoricalTree, c(181,180))
    collapsedClades$Elephant= MRCA(commonCategoricalTree, c(178,177,176))
    collapsedClades$Aardvark = MRCA(commonCategoricalTree, c(175,174,171,172,173))
    collapsedClades$Strepsirrhini = MRCA(commonCategoricalTree, c(28,29,27,26,24,25,23,21,22,20,19))
    collapsedClades$Atelidae = MRCA(commonCategoricalTree, c(6,5,4,3,2,1))
    collapsedClades$Chimpanze = MRCA(commonCategoricalTree, c(18,17,16,15,14,12,13,10,9,8,7,11))
    collapsedClades$Hares = MRCA(commonCategoricalTree, c(71,70))
    collapsedClades$Squirrels = MRCA(commonCategoricalTree, c(69,67,68,61,66,65,63,62,64))
    collapsedClades$Capybara = MRCA(commonCategoricalTree, c(33,32,31,30))
    collapsedClades$Beaver = MRCA(commonCategoricalTree, c(36,35,34))
    collapsedClades$Jerboa = MRCA(commonCategoricalTree, c(59,58,57))
    collapsedClades$Deomyinae = MRCA(commonCategoricalTree, c(41,39,40,38,37,44,43,42))
    collapsedClades$Vole = MRCA(commonCategoricalTree, c(48,47,46,45))
    collapsedClades$Neotominae = MRCA(commonCategoricalTree, c(54,49,50,52,51,53))
    collapsedClades$`African Hedgehogs` = MRCA(commonCategoricalTree, c(168,169))
    collapsedClades$`Talpa europaea` = MRCA(commonCategoricalTree, c(167,165,164,166))
    collapsedClades$`Flying Fox` = MRCA(commonCategoricalTree, c(163,161,162))
    collapsedClades$Rhinolophidae = MRCA(commonCategoricalTree, c(155,156,158,157,159,160))
    collapsedClades$`Big Brown Bat` = MRCA(commonCategoricalTree, c(143,142,139,140,141))
    collapsedClades$Phyllostomidae = MRCA(commonCategoricalTree, c(144,145,147,146))
    collapsedClades$Noctilio = MRCA(commonCategoricalTree, c(154))
    collapsedClades$Horse = MRCA(commonCategoricalTree, c(100,99,98))
    collapsedClades$Pig = MRCA(commonCategoricalTree, c(95,96))
    collapsedClades$`bos bison` = MRCA(commonCategoricalTree, c(87,89,88))
    collapsedClades$`Humpback Whale` = MRCA(commonCategoricalTree, c(75,74,72,73))
    collapsedClades$`Dolphins` = MRCA(commonCategoricalTree, c(82,81,80,77,76))
    collapsedClades$`Pangolin` = MRCA(commonCategoricalTree, c(138,137))
    collapsedClades$`Lion` = MRCA(commonCategoricalTree, c(131,130,129))
    collapsedClades$`Meerkat` = MRCA(commonCategoricalTree, c(135,134))
    collapsedClades$`Dog` = MRCA(commonCategoricalTree, c(101,102))
    collapsedClades$`Brown Bear` = MRCA(commonCategoricalTree, c(125,128,127))
    collapsedClades$`Odobenus rosmarus` = MRCA(commonCategoricalTree, c(120,119,118,117))
    collapsedClades$`Phocidae` = MRCA(commonCategoricalTree, c(124,123,121,122))
    collapsedClades$`Procyon lotor` = MRCA(commonCategoricalTree, c(114,113,112))
    collapsedClades$`Lontra provocax` = MRCA(commonCategoricalTree, c(108,107,106,105,104))
    collapsedClades$`Tasmanian Devil` = MRCA(commonCategoricalTree, c(184,185))
  }
  
  
  cladeUuids = NA
  for(i in 1:length(names(collapsedClades))){
    message(i)
    autoName = autocomplete_name(names(collapsedClades)[i])[1,2]
    message(autoName)
    cladeUuids[i] = get_uuid(autoName)
  }
  
  
  clades = data.frame(
    node = unlist(as.vector(collapsedClades[1,])),
    cladelabel = names(collapsedClades),
    uuid = cladeUuids,
    phylopic = cladeUuids
    
  )
  
  clades$x_new = NA
  clades$y_new = NA
  baseOffset = 0.04
  for(i in 1:nrow(clades)){
    cladeNode = clades$node[i]
    
    message(clades$cladelabel[i])
    message(cladeNode)
  
    coords = ggTreeClades$data[ggTreeClades$data$node == cladeNode, ]
    
    #Determine the longest branch length from the MRCA node, based on the longest tip species 
    tipDescendants = Descendants(commonCategoricalTree, cladeNode, type = "tip")[[1]]
    totalDistances = NULL
    
    for(j in 1:length(tipDescendants)){
      totalLength = 0 
      currentAncestor = tipDescendants[j]
      message("Child path:", j)
      while(currentAncestor != cladeNode){
        message(paste("processing Node:", currentAncestor))
        currentLength = commonCategoricalTree$edge.length[which(commonCategoricalTree$edge[, 2] == currentAncestor)]
        message(paste("Length:", currentLength))
        if(currentLength > totalLength/100){totalLength = totalLength + currentLength}; message("adding")
        message(paste("Total Length:", totalLength))
        currentAncestor = Ancestors(commonCategoricalTree, currentAncestor, "parent")
        if(length(currentAncestor) == 0){break}
      }
      message(paste("Total Length:", totalLength))
      message("-")
      totalDistances = append(totalDistances, totalLength)
      
    }
    maxLength = max(totalDistances)
  
    message(paste("Max Length:", maxLength))
    x_new = coords$x + maxLength + baseOffset
    y_new = coords$y
    clades$x_new[i] = x_new
    clades$y_new[i] = y_new
    message(x_new)
    message(y_new)
    message("------")
    
  }
}

ggTreeClades = ggtree(commonCategoricalTree, layout = "circular") +scale_color_manual(values=palette()) 
ggTreeClades = ggTreeClades %<+% edge + aes(color=CategorylengthChar) +labs(color = legendText)
ggTreeClades = ggTreeClades %<+% tip_data 
ggTreeClades$data$label = paste(ggTreeClades$data$label, "-", ggTreeClades$data$node, sep="")

ggTreeClades= ggTreeClades %<+% clades

if(!imageAllTips){
  for(i in 1:ncol(collapsedClades)){
    #ggTreeClades = ggTreeClades + geom_cladelabel(collapsedClades[1,i], names(collapsedClades)[i])
    #ggTreeClades = ggTreeClades + geom_cladelabel(collapsedClades[1,i], label = phylopiclist[i], color = "darkgray", barsize = 1, geom = "label", parse = T)
    ggTreeClades = ggTreeClades + geom_cladelabel(collapsedClades[1,i], label = NA, color = "darkgray", barsize = 1, geom = "label", parse = T)
  }
}

#output the tree
treeOutputLocation = paste0(outputFolderName, filePrefix, "RadialDisplayTree.pdf")
pdf(treeOutputLocation)
#ggTreeClades + rphylopic::geom_phylopic(data = ggTreeClades$data, aes(uuid = phylopic, x = x_new, y= y_new),size = 0.02)
ggTreeClades + ggimage:: geom_phylopic(data = ggTreeClades$data, aes(image = phylopic, x = x_new, y= y_new),size = 0.02)
dev.off()

treeOutputLocation = paste0(outputFolderName, filePrefix, "RadialDisplayTree.png")
png(treeOutputLocation, 3000, 3000)
#ggTreeClades + rphylopic::geom_phylopic(data = ggTreeClades$data, aes(uuid = phylopic, x = x_new, y= y_new),size = 0.02)
ggTreeClades + ggimage:: geom_phylopic(data = ggTreeClades$data, aes(image = phylopic, x = x_new, y= y_new),size = 0.02)
dev.off()





# Make code for using the clade style tip labels on individual tips 



