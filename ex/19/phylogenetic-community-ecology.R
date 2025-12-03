## install the picante package
install.packages("picante")

## load picante (and phytools)
library(picante)
library(phytools)

## download these files
download.file(
  url="https://liamrevell.github.io/biol634/data/SJ_ComMatrix.csv",
  destfile="SJ_ComMatrix.csv")
download.file(
  url="https://liamrevell.github.io/biol634/data/SJtree.phy",
  destfile="SJtree.phy")

## read in the tree
SJ_tree<-read.tree(file="SJtree.phy")
SJ_tree

## graph this tree
plotTree(SJ_tree,ftype="i",fsize=0.6)

## read in my data
SJ_data<-read.csv(file="SJ_ComMatrix.csv",row.names=1)
dim(SJ_data)

## plot presence/absence of each species from the
## regional pool at the tips of the tree
plotFanTree.wTraits(SJ_tree,SJ_data,type="fan",
  ftype="off",
  colors=replicate(ncol(SJ_data),
    setNames(c("#CCCDCB","black"),0:1),
    simplify=FALSE))
legend("topleft",c("present","absent"),
  lwd=4,col=c("black","#CCCDCB"))

## let's calculate PD from all our communities
?pd
SJ_all_pd<-pd(t(SJ_data),SJ_tree,include.root=TRUE)
head(SJ_all_pd)

## just to see that this is equal to the sum of the edge
## lengths, we can do an experiment using one of our
## communities: "Aleck_Rock"
species_on_Aleck_Rock<-which(SJ_data[,"Aleck_Rock"]==1)
species_on_Aleck_Rock<-rownames(SJ_data)[species_on_Aleck_Rock]
species_on_Aleck_Rock

## first drop all but the species on Aleck Rock
Aleck_Rock.tree<-keep.tip(SJ_tree,species_on_Aleck_Rock)

## unfortunately, we lose our stem here and there's no
## easy way to retain it, so we need to add it back
Aleck_Rock.tree$root.edge<-max(nodeHeights(SJ_tree))-
  max(nodeHeights(Aleck_Rock.tree))

## convert our root edge to an unbranching node
Aleck_Rock.tree<-rootedge.to.singleton(Aleck_Rock.tree)

## graph our two trees to show that we did it correctly
## (actually, it turns out that the community tree for Aleck
## Rock includes the global root!)
par(mfrow=c(2,1))
h<-max(nodeHeights(SJ_tree))
plotTree(SJ_tree,mar=c(2.1,1.1,1.1,1.1),ftype="off",lwd=1)
axis(1,at=h-seq(0,600,by=100),seq(0,600,by=100))
plotTree(Aleck_Rock.tree,mar=c(2.1,1.1,1.1,1.1),
  ftype="off",lwd=1)
axis(1,at=h-seq(0,600,by=100),seq(0,600,by=100))

## compute PD for Aleck Rock
sum(Aleck_Rock.tree$edge.length) ## matches picante::pd

## test hypothesis about PD using ses.pd
?ses.pd
SJ_pd.test<-ses.pd(t(SJ_data),SJ_tree,
  null.model="richness")
