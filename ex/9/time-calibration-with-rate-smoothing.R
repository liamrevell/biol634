## load package
library(phangorn)
library(phytools)

## https://liamrevell.github.io/biol634/data/LaurasiatherianML.nex
## read ML tree in from file
ml_tree<-read.nexus(
  file="https://liamrevell.github.io/biol634/data/LaurasiatherianML.nex")

## plot tree (it's unrooted)
plotTree(ml_tree)

## root tree with platypus
rooted.ml_tree<-root(ml_tree,outgroup="Platypus")
plotTree(rooted.ml_tree)

## drop platypus
rooted.ml_tree<-drop.tip(rooted.ml_tree,"Platypus")
plotTree(rooted.ml_tree)

## we're going to use these upper & lower bounds on the
## common ancestors for the following species pairs, originally
## from https://timetree.org/ (but perhaps out-of-date)
## Possum & Cat: 159, 166
## Squirrel & Mouse: 66, 75
## Pig & BlueWhale: 59, 66
## Human & Baboon: 27.95, 31.35
## Horse & Donkey: 6.2, 10

## function to find MRCA (most-recent-common-ancestor) of species
## pairs
?findMRCA

## compute all node indices
nodes<-c(
  findMRCA(rooted.ml_tree,c("Possum","Cat")),
  findMRCA(rooted.ml_tree,c("Squirrel","Mouse")),
  findMRCA(rooted.ml_tree,c("Pig","BlueWhale")),
  findMRCA(rooted.ml_tree,c("Human","Baboon")),
  findMRCA(rooted.ml_tree,c("Horse","Donkey")))
nodes

## check individual nodes
nodelabels(text="is this right?",node=47)
nodelabels(text="is this correct?",node=nodes[2])

## check all nodes
plotTree(rooted.ml_tree)
nodelabels(node=nodes)

## set min & max ages
age_min<-c(159,66,59,27.95,6.2)
age_max<-c(166,75,66,31.35,10)

## function we'll use for our calibration data frame
?makeChronosCalib

## make our calibration data frame
calibration<-makeChronosCalib(
  rooted.ml_tree,
  node=nodes,
  age.min=age_min,
  age.max=age_max)

## conduct penalized-likelihood rate smoothing using
## default penalty term (lambda = 1)
pl_tree<-chronos(rooted.ml_tree,calibration=calibration)

## compare rate-smoothed / calibrated tree to original
## uncalibrated ML tree
par(mfrow=c(1,2))
plotTree(rooted.ml_tree,fsize=0.7)
plotTree(pl_tree,fsize=0.7)

## try lambda = 0.1 penalty term (lower penalty for rate
## heterogeneity)
pl_tree0.1<-chronos(rooted.ml_tree,calibration=calibration,
  lambda=0.1)

## replot with comparison
par(mfrow=c(1,2))
plotTree(pl_tree,fsize=0.7)
plotTree(pl_tree0.1,fsize=0.7)

## try the phytools::compare.chronograms function
dev.off()
compare.chronograms(pl_tree,pl_tree0.1)

## try lambda = 10
pl_tree10<-chronos(rooted.ml_tree,calibration=calibration,
  lambda=10)

## compare again
compare.chronograms(pl_tree,pl_tree10)

## re-plot our lambda = 1 calibrated tree, but with
## lines showing our calibration intervals.
dev.off()
plotTree(pl_tree,mar=c(2.1,2.1,0.1,0.1))
axis(1)

par()$usr

## this is the value we'll use in our offset
192.951687-max(nodeHeights(pl_tree))

## plot tree "backwards" but facing left
plotTree(pl_tree,xlim=c(max(nodeHeights(pl_tree)),-33),
  direction="leftwards",
  mar=c(2.1,1.1,1.1,1.1))
axis(1)

## add calibration intervals
abline(v=age_min,lty="dotted",col=palette()[1:5])
abline(v=age_max,lty="dotted",col=palette()[1:5])

nodelabels(node=nodes,cex=0.6)
