## this code uses the MP and ML trees from HW 2 to explore tree 
## comparison and other things

## MP tree
https://liamrevell.github.io/biol634/data/Goby_MP.tre

## ML tree
https://liamrevell.github.io/biol634/data/Goby_ML.tre

## load packages
library(phangorn)
library(phytools)

## read MP tree from the web
goby_mp<-read.tree(
  file="https://liamrevell.github.io/biol634/data/Goby_MP.tre")
goby_mp

## read ML tree from the web
goby_ml<-read.tree(
  file="https://liamrevell.github.io/biol634/data/Goby_ML.tre")
goby_ml

## compute co-phylogenetic plot (without rooting)
goby_cophylo<-cophylo(goby_mp,goby_ml)
plot(goby_cophylo,fsize=0.6)

## set outgroup
outgroup<-"Oreochromis_niloticus"

## root MP tree
## we need to use edgelabel=TRUE to ensure that the bootstrap
## value node labels are maintained with the correct edges during
## re-rooting
rooted.goby_mp<-root(goby_mp,outgroup=outgroup,
  edgelabel=TRUE,resolve.root=TRUE)
rooted.goby_mp

## root ML tree
rooted.goby_ml<-root(goby_ml,outgroup=outgroup,
  edgelabel=TRUE,resolve.root=TRUE)
rooted.goby_ml

## plot trees (we could also using phangorn::plotBS)
plotTree(rooted.goby_mp,fsize=0.7,type="cladogram")
plotTree(rooted.goby_ml,fsize=0.7)

## recalculate co-phylogenetic plot
goby_cophylo<-cophylo(rooted.goby_mp,rooted.goby_ml)
plot(goby_cophylo,pts=FALSE,fsize=0.6)

## calculate Robinson-Foulds distance
RF.dist(rooted.goby_mp,rooted.goby_ml)

## prune out Schindleria because it seems to be affected
## by long-branch-attraction
nn<-rooted.goby_mp$tip.label
tip_to_prune<-nn[grep("Schindleria",nn)]
tip_to_prune

## prune both trees usind drop.tip
pruned_rooted.goby_mp<-drop.tip(rooted.goby_mp,
  tip_to_prune)
pruned_rooted.goby_ml<-drop.tip(rooted.goby_ml,
  tip_to_prune)

## re-plot co-phylo to make sure we got it correct
plot(cophylo(pruned_rooted.goby_mp,
  pruned_rooted.goby_ml),pts=FALSE)

## recalculate RF distance
RF.dist(pruned_rooted.goby_mp,pruned_rooted.goby_ml)

## next we want to conduct SH test

## read original data from web
goby_dat<-read.phyDat(file=
    "https://liamrevell.github.io/biol634/data/concat.goby_75pct_JE2.nex.fas",
  format="fasta")
goby_dat

## fit "pml" object with ML tree
goby_ml.pml<-pml(goby_ml,data=goby_dat,k=4)
optimized.goby_ml.pml<-optim.pml(
  goby_ml.pml,optBf=TRUE,optQ=TRUE,
  optGamma=TRUE,optEdge=FALSE)

## do the same thing, but using MP tree
goby_mp.pml<-pml(
  compute.brlen(goby_mp),data=goby_dat,k=4)
optimized.goby_mp.pml<-optim.pml(
  goby_mp.pml,optBf=TRUE,optQ=TRUE,
  optGamma=TRUE,optEdge=TRUE)

## calculate co-phylo and plot
optimized.cophylo<-cophylo(
  root(optimized.goby_ml.pml$tree,
    outgroup=outgroup,resolve.root=TRUE),
  root(optimized.goby_mp.pml$tree,
    outgroup=outgroup,resolve.root=TRUE))
plot(optimized.cophylo,pts=FALSE,fsize=0.6)

## run SH test
goby_shtest<-SH.test(optimized.goby_ml.pml,
  optimized.goby_mp.pml)
goby_shtest
