## get entire path to our morphological trait matrix file
file_path<-file.choose() ## then select file
file_path

## load packages
## update phytools
## install.packages("phytools") ## if you want to update phytools
library(phytools)
library(phangorn)

## read our input data
vert_data<-read.csv(file_path,row.names=2)
head(vert_data)

## save common names
common_names<-vert_data$Taxon
names(common_names)<-rownames(vert_data)
common_names

## remove the common names from our data frame
vert_data<-vert_data[,-1]

## convert our data frame to a phyDat object
vert_phyDat<-as.phyDat(as.matrix(vert_data),
  type="USER",levels=c(0,1))
vert_phyDat

## estimate a tree using parsimony
mp_vert<-pratchet(vert_phyDat)
mp_vert

## plot unrooted tree with radial labels
plot(mp_vert,type="unrooted",lab4ut="axial",
  cex=0.8,no.margin=TRUE)

rooted.mp_vert<-root(mp_vert,
  outgroup=names(common_names[which(common_names=="Lamprey")]),
  resolve.root=TRUE)

## plot our newly rooted tree
plotTree(rooted.mp_vert,lwd=1,ftype="i")

## get the "true" tree to compare
craniate_timetree<-read.tree(
  file="https://liamrevell.github.io/biol634/data/craniata_tree.nwk")
craniate_timetree
plotTree(craniate_timetree,ftype="i",lwd=1)

## create cophylo plot
vert_cophylo<-cophylo(rooted.mp_vert,
  craniate_timetree)
plot(vert_cophylo,link.type="curved")

## substitute names
ind<-grep("Gallus",craniate_timetree$tip.label)
craniate_timetree$tip.label[ind]<-
  rooted.mp_vert$tip.label[grep("Gallus",
    rooted.mp_vert$tip.label)]

ind<-grep("Physeter",craniate_timetree$tip.label)
craniate_timetree$tip.label[ind]<-
  rooted.mp_vert$tip.label[grep("Physeter",
    rooted.mp_vert$tip.label)]

ind<-grep("Dryophytes",craniate_timetree$tip.label)
craniate_timetree$tip.label[ind]<-
  rooted.mp_vert$tip.label[grep("Hyla",
    rooted.mp_vert$tip.label)]

## get the parsimony score
attr(rooted.mp_vert,"pscore")

## compute the parsimony score on our "true" tree
parsimony(craniate_timetree,vert_phyDat)

## estimate the tree using likelihood
rVert<-rtree(n=length(vert_phyDat),
  tip.label=names(vert_phyDat),
  rooted=FALSE)
vert_pml<-pml(rVert,data=vert_phyDat)
vert_pml

mk.vert_pml<-optim.pml(vert_pml,rearrangement="stochastic")
mk.vert_pml
plotTree(mk.vert_pml$tree,lwd=1)

ml_tree<-root(mk.vert_pml$tree,outgroup="Petromyzon_marinus",
  resolve.root=TRUE)
plotTree(compute.brlen(ml_tree),ftype="i")

vert_cophylo2<-cophylo(compute.brlen(ml_tree),
  craniate_timetree)
plot(vert_cophylo2,link.type="curved")
