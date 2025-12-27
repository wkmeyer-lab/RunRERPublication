?foreground2Tree

mainTrees = readTrees("Data/DemoUNICORNs.txt") 

annotationsTable = read.csv("Data/MyDentalAnnotations.csv")

treeOutputFile = "Output/thisSpecificOUtputTree.rds"

# Lets pick one collumn to use first 

phenotypeColumn = "toothShape1"
phenotypeColumn = "toothShape2"



#Let's only use the species in that column for now 

relevantSpecies = annotationsTable[annotationsTable[[phenotypeColumn]] %in% c(0,1)  ]

relevantSpecies # this is basically our whole data 


#Let's make a speceis fliter, so that we don't waste time claculating RERs on species we aren't using 
speciesFilter = relevantSpecies$nameAsItAppearsOnTheTree


#Now, we need to pick which of the species is in the foreground 
foregroundSpecies = relevantSpecies[ relevantSpecies[[phenotypeColumn]] == 1]

foregroundNames = foregroundSpecies$nameAsItAppearsOnTheTree



binaryPhenotypeTree = foreground2Tree(
  foreground = foregroundNames, 
  treesObj = mainTrees, 
  clade = "all", 
  weighted = F, 
  transition = "bidirectional", 
  useSpecies = speciesFilter
  )


#save the tree
saveRDS(binaryPhenotypeTree, treeOutputFile)

#plot the tree to check it 
pdf("output/TestPDF.pdf")
foreground2Tree(
  foreground = foregroundNames, 
  treesObj = mainTrees, 
  clade = "all", 
  weighted = F, 
  transition = "bidirectional", 
  useSpecies = speciesFilter
)
dev.off()


#------ 

mainTrees$mastertrees$tip.label 
#make sure you have a column in your phenotype speadsheet with these for each species 