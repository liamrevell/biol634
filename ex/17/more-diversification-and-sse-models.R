## load phytools
library(phytools)

## load primate trees
primate_trees<-read.nexus(
  file="https://liamrevell.github.io/biol634/data/10kTrees_Primates.nex")
primate_trees

## compute one tree from this set: the MCC tree
library(phangorn)
?maxCladeCred
primate_mcc<-maxCladeCred(primate_trees)
primate_mcc

## let's plot this tree
plotTree(primate_mcc,ftype="i",fsize=0.5,type="fan",
  lwd=1)

## how many primates are there?
nlow<-376
nhigh<-524

## fit our bd model (fails)
primates_bd.low<-fit.bd(primate_mcc,
  rho=Ntip(primate_mcc)/nlow)

## is our tree ultrametric
plotTree(primate_mcc,ftype="off",lwd=1)

?force.ultrametric

## apply force.ultrametric
primate_mcc<-force.ultrametric(primate_mcc)

## now run our bd model
primates_bd.low<-fit.bd(primate_mcc,
  rho=Ntip(primate_mcc)/nlow)
primates_bd.low

## now run with high n
primates_bd.high<-fit.bd(primate_mcc,
  rho=Ntip(primate_mcc)/nhigh)
primates_bd.high

## test to see what happens with an unrealistically
## high n
n_unrealistic<-1000
fit.bd(primate_mcc,
  rho=Ntip(primate_mcc)/n_unrealistic)

## download files
download.file(
  url="https://liamrevell.github.io/biol634/data/grunts.csv",
  destfile="grunts.csv")
download.file(
  url="https://liamrevell.github.io/biol634/data/grunts.phy",
  destfile="grunts.phy")

## let's read our tree from file
grunt_tree<-read.tree(file="grunts.phy")
grunt_tree

## get our data
grunt_data<-read.csv(file="grunts.csv",row.names=1)
head(grunt_data)

## pull out habitat
grunt_habitat<-setNames(grunt_data$habitat,
  rownames(grunt_data))
head(grunt_habitat)

## plot our tree
plotTree(grunt_tree,ftype="i",fsize=0.8,
  offset=0.5)
tiplabels(
  pie=to.matrix(grunt_habitat,0:1)[grunt_tree$tip.label,],
  cex=0.4,piecol=hcl.colors(n=2))
legend("bottomleft",c("non-reef","reef"),pch=21,
  pt.bg=hcl.colors(n=2),pt.cex=2,bty="n")

## need to use the diversitree package
install.packages("diversitree") ## only if you don't have it installed!
library(diversitree)

## make our BiSSE likelihood function
grunt_bisse.fn<-make.bisse(grunt_tree,grunt_habitat)
grunt_bisse.fn

## let's get starting values for optimization
starting_p<-starting.point.bisse(grunt_tree)
starting_p

## let's get our total tree height
max(nodeHeights(grunt_tree))

## get our total tree height
h<-max(nodeHeights(grunt_tree))
h
rescaled.grunt_tree<-grunt_tree
rescaled.grunt_tree$edge.length<-rescaled.grunt_tree$edge.length/h*45

## let's check if it worked!
plotTree(rescaled.grunt_tree,mar=c(3.1,1.1,1.1,1.1))
axis(1)

## make our BiSSE likelihood function
grunt_bisse.fn<-make.bisse(rescaled.grunt_tree,
  grunt_habitat)
grunt_bisse.fn

## let's get starting values for optimization
starting_p<-starting.point.bisse(rescaled.grunt_tree)
starting_p

## compare starting points to fitted birth-death model
fit.bd(rescaled.grunt_tree)

## optimize BiSSE model
grunt_bisse.mle<-find.mle(grunt_bisse.fn,starting_p)
grunt_bisse.mle

## print the parameter values alone
round(grunt_bisse.mle$par,6)
