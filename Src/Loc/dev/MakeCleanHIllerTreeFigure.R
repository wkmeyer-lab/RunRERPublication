library(ggtree)

hiller4MasterCommon = commonMainTrees$masterTree

commonCategoricalTree = char2TreeCategorical(commonPhenotypeVector, commonMainTrees, commonSpeciesFilter, model = modelType, anctrait = ancestralTrait, plot = T)
nodelabels(col = "red", adj = c(0, -0), frame = "none")
tiplabels(col = "blue", frame = "none")
edgelabels(col = "darkgreen", frame = "none")



hiller4MasterCommonTrimmed = drop.tip(hiller4MasterCommon, hiller4MasterCommon$tip.label[!hiller4MasterCommon$tip.label %in% commonCategoricalTree$tip.label])

commonCategoricalTreeEdgeLengths = commonCategoricalTree$edge.length
commonCategoricalTreeEdgeLengths = append(commonCategoricalTreeEdgeLengths, 1)

commonCategoricalTreeEdgeLengths = as.character(commonCategoricalTreeEdgeLengths)

length(commonCategoricalTree$edge.length)


edge=data.frame(hiller4MasterCommonTrimmed$edge, edge_num=1:nrow(hiller4MasterCommonTrimmed$edge))
colnames(edge)=c("parent", "node", "edge_num")
edge$Categorylength = commonCategoricalTree$edge.length
edge$CategorylengthChar = as.character(edge$Categorylength)

catEdge=data.frame(commonCategoricalTree$edge, edge_num=1:nrow(commonCategoricalTree$edge))
colnames(catEdge)=c("parent", "node", "edge_num")

all.equal(edge, catEdge)

edge$length = hiller4MasterCommonTrimmed$edge.length



edge$length = hiller4MasterCommonTrimmed$edge.length

rm(ggTreeOut)

ggTreeOut = ggtree(hiller4MasterCommonTrimmed) + geom_tiplab() + scale_color_manual(values=c("black", "firebrick", "darkgreen", "steelblue"))

ggTreeOut = ggTreeOut %<+% edge + aes(color=CategorylengthChar)

ggTreeOut


# Collapse some clades
{
  collapsedClades = data.frame()
  collapsedClades[1,] = NA
  
  collapsedClades$Bovidae = MRCA(hiller4MasterCommonTrimmed, c("Zebu", "Asian Water Buffalo"))
  
  collapsedClades$Caprids = MRCA(hiller4MasterCommonTrimmed, c("Bighorn sheep", "Tibetan antelope"))
  
  collapsedClades$Deer = MRCA(hiller4MasterCommonTrimmed, c("Central European Red Deer", "White-tailed Deer"))
  
  collapsedClades$Cetacea = MRCA(hiller4MasterCommonTrimmed, c("Bottle-nose Dolphin", "Sperm whale"))
  
  collapsedClades$Camelidae = MRCA(hiller4MasterCommonTrimmed, c("Ferus Camel", "Alpaca"))
  
  collapsedClades$Pinnipeds = MRCA(hiller4MasterCommonTrimmed, c("Hawaiian Monk Seal", "Pacific walrus"))
  
  collapsedClades$Mustelidae = MRCA(hiller4MasterCommonTrimmed, c("Ferret", "Red Panda "))
  
  collapsedClades$Bears = MRCA(hiller4MasterCommonTrimmed, c("Polar Bear", "Panda"))
  
  collapsedClades$Dogs = MRCA(hiller4MasterCommonTrimmed, c("African Hunting Dog ", "Dog"))
  
  collapsedClades$Cats = MRCA(hiller4MasterCommonTrimmed, c("Tiger", "Cat"))
  
  collapsedClades$Pangolins = MRCA(hiller4MasterCommonTrimmed, c("Chinese pangolin", "Sunda pangolin"))
  
  collapsedClades$Horses = MRCA(hiller4MasterCommonTrimmed, c("Horse", "Wild Donkey"))
  
  collapsedClades$Megabats = MRCA(hiller4MasterCommonTrimmed, c("Black flying-fox", "Rousette Fruit Bats"))
  
  collapsedClades$leaftnosedBats = MRCA(hiller4MasterCommonTrimmed, c("Horseshoe Bats", "Old World Leaf-nosed Bats"))
  
  collapsedClades$bigBrownBat = MRCA(hiller4MasterCommonTrimmed, c("Microbat", "Long-winged bats"))
  
  collapsedClades$Eulipotyphla = MRCA(hiller4MasterCommonTrimmed, c("Hedgehog", "Shrew"))
  
  collapsedClades$Maccaca = MRCA(hiller4MasterCommonTrimmed, c("Rhesus Macaca", "Green monkey"))
  
  collapsedClades$Langur = MRCA(hiller4MasterCommonTrimmed, c("Proboscis Monkey", "Black-and-white Colobus Monkey"))
  
  collapsedClades$Apes = MRCA(hiller4MasterCommonTrimmed, c("Chimp", "Gibbon"))
  
  collapsedClades$Monkeys = MRCA(hiller4MasterCommonTrimmed, c("Squirrel monkey", "Night Monkey"))
  
  collapsedClades$Lemur = MRCA(hiller4MasterCommonTrimmed, c("Sifakas", "Bushbaby"))
  
  collapsedClades$Rat = MRCA(hiller4MasterCommonTrimmed, c("Ryukyu mouse", "Rat"))
  
  collapsedClades$Hamster = MRCA(hiller4MasterCommonTrimmed, c("Chinese hamster", "Deer mouse"))
  
  collapsedClades$Beaver = MRCA(hiller4MasterCommonTrimmed, c("Ord Kangaroo Rat", "Beaver"))
  
  collapsedClades$Chinchilla = MRCA(hiller4MasterCommonTrimmed, c("Brush-tailed rat", "Domestic guinea pig"))
  
  collapsedClades$Molerats = MRCA(hiller4MasterCommonTrimmed, c("Naked mole-rat", "Damara mole rat"))
  
  collapsedClades$Squirrel = MRCA(hiller4MasterCommonTrimmed, c("Ground Squirrel ", "Marmot"))
  
  collapsedClades$Rabbit = MRCA(hiller4MasterCommonTrimmed, c("Rabbit", "Pika"))
  
  collapsedClades$Elephant = MRCA(hiller4MasterCommonTrimmed, c("Manatee", "Hydrax"))
  
  collapsedClades$Aardvark = MRCA(hiller4MasterCommonTrimmed, c("Cape elephant shrew", "Aardvark"))
  
  collapsedClades$Marsupials = MRCA(hiller4MasterCommonTrimmed, c("Tasmanian devil", "Opossum"))
}
#

ggTreeClades = ggtree(hiller4MasterCommonTrimmed, layout = "circular", linewidth = 1) + scale_color_manual(values=c("black", "firebrick", "lightgreen", "steelblue"))

ggTreeClades = ggTreeClades %<+% edge + aes(color=CategorylengthChar)

for(i in 1:ncol(collapsedClades)){
#  ggTreeClades = ggTreeClades + geom_cladelabel(collapsedClades[1,i], names(collapsedClades)[i])
  ggTreeClades = ggTreeClades + geom_cladelabel(collapsedClades[1,i], NA)
}
ggTreeClades

png("Output/NewHiller4Phen/NewHiller4PhenCicleTree.png", width = 2000, height = 2000)
ggTreeClades
dev.off()

?png
?geom_cladelabel

hiller4MasterCommonTrimmed$tip.label

collapsedClades = marsupials

ggTreeOut2 = ggtree(hiller4MasterCommonTrimmed) + geom_tree()

ggTreeOut 
