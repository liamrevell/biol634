## load phangorn package
library(phangorn)

## read our data from file using read.phyDat
mammal_dna<-read.phyDat(file="Laurasiatherian.nex",
  format="nexus")
mammal_dna

## look at the structure of the object
str(mammal_dna)
names(mammal_dna)->our_taxa
our_taxa

## create a random tree with these taxa
random_tree<-rtree(n=length(our_taxa),
  tip.label=our_taxa,br=NULL,rooted=FALSE)
random_tree

## plot this tree
plot(random_tree,type="unrooted",lab4ut="axial",
  cex=0.9,no.margin=TRUE)

## estimate a tree using optimization via NNI
nni_tree<-optim.parsimony(random_tree,mammal_dna,
  rearrangements="NNI",trace=2)

## estimate a tree using optimization via a combo of
## SPR and NNI
spr_tree<-optim.parsimony(random_tree,mammal_dna,
  rearrangements="SPR",trace=2)

### plot our tree
spr_tree
phytools::plotTree(spr_tree,fsize=0.8)

## root the tree using Platypus
rooted.spr_tree<-root(spr_tree,outgroup="Platypus")
phytools::plotTree(rooted.spr_tree,fsize=0.8)
