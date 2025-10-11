## this is showing how to run a bootstrapping analysis with neighbor-joining

## preliminaries

## load packages
library(phytools)
library(phangorn)

## read data from file
dna<-read.dna(file=
    "https://liamrevell.github.io/biol634/data/primates.dna")

## get NJ from this dataset
nj_tree<-NJ(dist.dna(dna,model="JC69"))
plot(nj_tree,type="unrooted")

## step 1: bootstrap the data 100 times
nrep<-100 ## this is my number of replicates
nsites<-ncol(dna) ## this is the number of sites in my alignment
bs_dna<-vector(mode="list",length=nrep)
for(i in 1:nrep) 
  bs_dna[[i]]<-dna[,sample(1:nsites,replace=TRUE)]

## inspect our object to make sure we're doing things correctly
## (this will mainly tell us if we have done things incorrectly)
length(bs_dna)
head(bs_dna)

## step 2: compute NJ tree on all 100 alignments
bs_trees<-vector(mode="list",length=nrep)
for(i in 1:nrep){
  D<-dist.dna(bs_dna[[i]],model="JC69")
  bs_trees[[i]]<-NJ(D)
}

## step 3 (optional): root all trees
## we can compute bootstrap proportions for unrooted trees, but
## rooting makes it more intuitive
outgroup<-"Lemur"
rooted.nj_tree<-root(nj_tree,outgroup=outgroup,
  resolve.root=TRUE)
plotTree(rooted.nj_tree)

## optional: root all the trees in our list using lapply
rooted.bs_trees<-lapply(bs_trees,FUN=root,
  outgroup=outgroup,resolve.root=TRUE)
par(mfrow=c(10,10))
nulo<-lapply(rooted.bs_trees,plotTree,fsize=0.5,lwd=1)
class(rooted.bs_trees)<-"multiPhylo" ## assign class attribute
rooted.bs_trees

## optional: root all trees using root (which is vectorized)
class(bs_trees)<-"multiPhylo"
rooted.bs_trees<-root(bs_trees,outgroup=outgroup,
  resolve.root=TRUE)
rooted.bs_trees

## try to plot all our BS trees in a grid
dev.off()
par(mfrow=c(10,10))
plotTree(rooted.bs_trees,fsize=0.6,lwd=1)

## step4 (optional): use ape::prop.clades to compute BS %
bs<-prop.clades(rooted.nj_tree,rooted.bs_trees,
  rooted=TRUE)
bs

dev.off() ## reset plotting device

## plot our bootstrap percent
plotTree(rooted.nj_tree)
nodelabels(bs)

## let's make it look a little nicer
plotTree(rooted.nj_tree)
nodelabels(bs/100,adj=c(1.2,1.2),cex=0.9,frame="none")

## we could've skipped all this....
plotBS(rooted.nj_tree,rooted.bs_trees)

## clean up our workspace
rm(list=ls())

## let's make a function for all this
njBoot<-function(dna,nrep=100,outgroup){
  nsites<-ncol(dna)
  foo<-function(X,nsites) X[,sample(1:nsites,replace=TRUE)]
  bs_dna<-replicate(nrep,foo(dna,nsites),simplify=FALSE)
  rooted.nj_tree<-root(NJ(dist.dna(dna,model="JC69")),
    outgroup=outgroup,resolve.root=TRUE)
  foo<-function(X) root(NJ(dist.dna(X,model="JC69")),
    outgroup=outgroup,resolve.root=TRUE)
  rooted.bs_trees<-lapply(bs_dna,foo)
  bs<-prop.clades(rooted.nj_tree,rooted.bs_trees)
  plotTree(rooted.nj_tree,lwd=1)
  nodelabels(bs/100,adj=c(1.2,1.2),cex=0.9,frame="none")
  rooted.nj_tree$node.label<-bs/100
  return(rooted.nj_tree)
}

## we can run a test on our original dataset
XX<-read.dna(file=
    "https://liamrevell.github.io/biol634/data/primates.dna")

njBoot(XX,outgroup="Lemur")->bootstrapped_tree

## now let's do it with the Jackman et al. data
## (we read it in using phangorn::read.phyDat to more easily convert
## to "DNAbin" in matrix)
YY<-read.phyDat(
  file="https://liamrevell.github.io/biol634/data/Jackman-etal.nex",
  format="nexus")
YY<-as.DNAbin(YY)

njBoot(YY,outgroup="Diplolaemus_darwinii")

