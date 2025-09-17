### Starting out by reviewing the MP homeworks assignment

## load packages
library(phangorn)
library(phytools)

## read data from file (directly from course page)
Jackman.data<-read.phyDat(
  file="https://liamrevell.github.io/biol634/data/Jackman-etal.nex",
  format="nexus")
Jackman.data

## generate starting tree for optimization
starting_tree<-rtree(n=length(Jackman.data),
  rooted=FALSE,
  tip.label=names(Jackman.data),
  br=NULL)

starting_tree

## try to find the MP solution via simple NNI
## operations from the random starting tree
## (fails to find a good solution)
mp.nni<-optim.parsimony(starting_tree,Jackman.data,
  rearrangements="NNI")

## try to find the MP solution via a mix of
## NNI and SPR moves from a random starting tree
## (finds much better tree)
mp.spr<-optim.parsimony(starting_tree,Jackman.data,
  rearrangements="SPR")

## root trees using out group & plot
mp.nni_rooted<-root(mp.nni,outgroup="Diplolaemus_darwinii",
  resolve.root=TRUE)

plotTree(mp.nni_rooted)

mp.spr_rooted<-root(mp.spr,outgroup="Diplolaemus_darwinii",
  resolve.root=TRUE)

plotTree(mp.spr_rooted)

## plot both trees in multi-panel figure
par(mfrow=c(1,2))
plotTree(mp.nni_rooted)
plotTree(mp.spr_rooted,direction="leftwards")

## plot both trees as "tanglegram"
mp_cophylo<-cophylo(mp.nni_rooted,mp.spr_rooted)

dev.off() ## reset plotting device

plot(mp_cophylo,fsize=0.6,link.type="curved")

?plot.cophylo

## estimate MP tree using parsimony ratchet
mp_ratchet<-pratchet(Jackman.data)

plotTree(mp_ratchet,ftype="i",fsize=0.9)

## root tree
mp_ratchet.rooted<-root(mp_ratchet,
  outgroup="Diplolaemus_darwinii",
  resolve.root=TRUE)

## plot MP tree from parsimony ratchet
plotTree(mp_ratchet.rooted,ftype="i",fsize=0.9)

## compare to published tree (requires image of published tree)
spp<-c("darwinii","acutirostris","occultus","luteogularis","equestris",
  "nicefori","microtus","agassizi","luciae","richardi","aeneus",
  "bartschi","vermiculatus","bahorucoensis","coelestinus","aliniger",
  "olssoni","insolitus","etheridgei","barbouri","lucius","guamuhaya",
  "chamaeleonides","cuvieri","christophei","barahonae","bimaculatus",
  "wattsi","brevirostris","distichus","krugi","cristatellus","acutus",
  "stratulus","marcanoi","strahmi","vanidicus","alutaceus","pumilis",
  "loysiana","carolinensis","maynardi","sheplani","paternus","angusticeps",
  "ahli","sagrei","ophiolepis","lineatus","limifrons","humilis",
  "valencienni","lineatopus","garmani","grahami")
spp<-sapply(spp[length(spp):1],function(x,y) y[grep(x,y)],
  y=mp_ratchet.rooted$tip.label)
tip<-setNames(1:length(spp),spp)
tree<-minRotate(minRotate(mp_ratchet.rooted,sort(tip)),
  sort(tip))

## plot comparison (requires jpeg package)
library(jpeg)
download.file(
  url="https://liamrevell.github.io/biol634/data/Jackman-etal-tree.jpg",
  destfile="Jackman-etal-tree.jpg",
  mode="wb")
img<-readJPEG(source="Jackman-etal-tree.jpg")
layout(matrix(c(1,2),1,2),widths=c(0.5,0.5))
plot.new()
par(mar=rep(0,4))
plot.window(xlim=c(0,1),ylim=c(0,1))
rasterImage(img,0,0,1,1)
plotTree(compute.brlen(tree,power=1),fsize=0.7,ftype="i",
  direction="leftwards",
  ylim=c(1,Ntip(tree))+c(-11,4.5))

### Now switch to Maximum Likelihood

## load data from web
Jackman.data<-read.phyDat(
  file="https://liamrevell.github.io/biol634/data/Jackman-etal.nex",
  format="nexus")

## random starting tree
starting_tree<-rtree(
  n=length(Jackman.data),
  tip.label=names(Jackman.data),
  rooted=FALSE)

## create "pml" object
unfitted.pml<-pml(tree=starting_tree,
  data=Jackman.data)

unfitted.pml

## optimize
Jackman_jc<-optim.pml(unfitted.pml,
  optEdge=TRUE,model="JC",
  rearrangement="NNI")

str(Jackman_jc)

## pull out & root tree
ml_tree<-Jackman_jc$tree
rooted.ml_tree<-root(ml_tree,outgroup="Diplolaemus_darwinii",
  resolve.root=TRUE)
dev.off()
plot(rooted.ml_tree)

## compare to parsimony tree
plot(cophylo(rooted.ml_tree,
  mp_ratchet.rooted),fsize=0.6)
