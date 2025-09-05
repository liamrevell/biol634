## install phytools from CRAN
install.packages("phytools")

## install package remotes
install.packages("remotes")

## install phytools dev version
remotes::install_github("liamrevell/phytools")

## to force
remotes::install_github("liamrevell/phytools",
  force=TRUE)

## Bolitoglossa, Felis, Bison_bison, Iguana_iguana, Sapajus_nigritis, Scaphiopus
"((((Felis,Bison_bison),Sapajus_nigritis),Iguana_iguana),(Bolitoglossa,Scaphiopus));"->
  newick_string

## load phytools
library(phytools)

## parse my Newick string
vertebrate_tree<-read.tree(text=newick_string)
vertebrate_tree

## object of class "phylo"
class(vertebrate_tree)

## call print function on object
print(vertebrate_tree)

## plot our tree
plot(vertebrate_tree)

?plot

## if the method is documented, the documentation will live at
## method_name.class_name
?plot.phylo

## make a simpler Newick tree with edge lengths
"(shark:462,(human:319,lizard:319):143);"->newick_string2

## graph our tree with edge lengths
vertebrate_tree_edge_lengths<-read.tree(text=newick_string2)
plot(vertebrate_tree_edge_lengths)
axis(1)

## see a couple of other plot styles supported by plot.phylo
plot(vertebrate_tree,type="cladogram") ## more accurately a "slanted cladogram
## or phylogram"
plot(vertebrate_tree_edge_lengths,type="cladogram")
plot(vertebrate_tree,type="unrooted")

## take a look inside the "phylo"
plot(vertebrate_tree)

str(vertebrate_tree)

nodelabels()
tiplabels()

## let's see just the "edge" matrix
vertebrate_tree$edge

unroot(vertebrate_tree)->unrooted.vertebrate_tree

plot(unrooted.vertebrate_tree)
plot(unrooted.vertebrate_tree,type="unrooted")

## read primate tree from file
primate_tree<-read.tree(file="primates.phy")

## read directly from the website
primate_tree<-read.tree(file="https://liamrevell.github.io/biol634/data/primates.phy")

## look at our primate tree
primate_tree

## plot the tree using defaults of plot.phylo
plot(primate_tree,cex=0.5)

## different function for plotting tree -- plotTree
plotTree(primate_tree,fsize=0.6,ftype="i")

## plotting in "fan" style
plotTree(primate_tree,fsize=0.6,ftype="i",type="fan")

## write to a vector based graphical
pdf(file="example-primate-tree.pdf",width=8,height=8)
plotTree(primate_tree,fsize=0.6,ftype="i",type="fan")
dev.off()
