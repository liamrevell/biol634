## in this first part we're just reiterating what we did at the
## end of class on Wednesday

## https://liamrevell.github.io/biol634/data/Jackman-etal.nex

## load packages
library(phangorn)
library(phytools)

## read Jackman et al. data from course website
Jackman.data<-read.phyDat(
  file="https://liamrevell.github.io/biol634/data/Jackman-etal.nex",
  format="nexus")
Jackman.data

## alternative: use download.file to create a local copy of the file
?download.file
download.file(
  url="https://liamrevell.github.io/biol634/data/Jackman-etal.nex",
  destfile="Jackman-etal.nex")
Jackman.data

## generate random starting tree
jackman_rand<-rtree(
  n=length(Jackman.data),
  tip.label=names(Jackman.data),
  rooted=FALSE)

## create an unfitted "pml" object (this is our unoptimized model)
jackman_unfitted.pml<-pml(
  tree=jackman_rand,
  data=Jackman.data)
jackman_unfitted.pml

## optimize under the Jukes-Cantor ("JC") model
jackman_jc_fitted.pml<-optim.pml(
  object=jackman_unfitted.pml,
  optEdge=TRUE,
  rearrangement="NNI")

## check class of fitted object
class(jackman_jc_fitted.pml)

## there is an S3 plot method for this object type
## but it plots the tree with midpoint root (we don't want that)
plot(jackman_jc_fitted.pml)

## let's pull of the "phylo" object to plot directly
str(jackman_jc_fitted.pml)
jackman_jc.phy<-jackman_jc_fitted.pml$tree
plot(jackman_jc.phy)

## outgroup root and re-plot
rooted.jackman_jc.phy<-root(jackman_jc.phy,
  outgroup="Diplolaemus_darwinii",resolve.root=TRUE)

## (also how to get the names off the "phyDat" object)
names(Jackman.data)

## plot rooted tree
plot(rooted.jackman_jc.phy) #,use.edge.length=FALSE)

## our fitted model
jackman_jc_fitted.pml

## now let's try to estimate with a different starting
## tree
jackman_mp<-pratchet(Jackman.data)

unclass(jackman_mp)

plot(jackman_mp)

## set arbitrary branch lengths for pml
jackman_mp$edge.length<-
  runif(n=nrow(jackman_mp$edge))

## another option
jackman_mp<-compute.brlen(jackman_mp)

plot(jackman_mp)

## create unfited "pml" object
jackman_unfitted.pml<-pml(
  tree=jackman_mp,
  data=Jackman.data)
jackman_unfitted.pml

## optimize
jackman_jc_fitted.pml<-optim.pml(
  object=jackman_unfitted.pml,
  optEdge=TRUE,
  rearrangement="NNI")

## pull tree and root
jackman_jc.phy<-jackman_jc_fitted.pml$tree
rooted.jackman_jc.phy<-root(
  jackman_jc.phy,
  outgroup="Diplolaemus_darwinii",
  resolve.root=TRUE)

## plot rooted tree
plotTree(rooted.jackman_jc.phy)

## now try with pml_bb (this should be best)
jackman_jc_fitted.pml.2<-pml_bb(Jackman.data,model="JC")

## compare to K80 model
jackman_k80_fitted.pml<-pml_bb(Jackman.data,model="K80")
