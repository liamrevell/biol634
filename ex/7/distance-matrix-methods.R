## start by loading a simple distance matrix (matrix from Felsenstein 2004, 
## p. 163)
D.carnivora<-read.csv(
  file="https://liamrevell.github.io/biol634/data/carnivore.distances.csv",
  row.names=1)
D.carnivora

## load phylogeny packages
library(phangorn)
library(ape)
library(phytools)

## help page of UPGMA function in phangorn
?upgma

## estimate a tree using UPGMA
upgma_carnivora<-upgma(D.carnivora) ## works
upgma_carnivora

## we could also have converted our input matrix into an object of the "dist" 
## class
upgma_carnivora<-upgma(as.dist(D.carnivora)) ## also works
upgma_carnivora

## graph tree
plotTree(upgma_carnivora,direction="upwards",fsize=1.5)

## just to see what the "dist" class looks like
as.dist(D.carnivora)

## get the NJ tree
nj_carnivora<-NJ(as.dist(D.carnivora))

## intentionally plot our tree incorrectly (NJ gives an unrooted tree)
plotTree(nj_carnivora)

## midpoint rooting
rooted.nj_carnivora<-midpoint(nj_carnivora)
plotTree(rooted.nj_carnivora)

## compute Robinson-Foulds distance between UPGMA and NJ trees
RF.dist(upgma_carnivora,rooted.nj_carnivora)

## typically we obtain distances from DNA sequences. Let's try that
primate_dna<-read.dna(
  file="https://liamrevell.github.io/biol634/data/primates.dna")
primate_dna

## use model="JC69"
D_primates<-dist.dna(primate_dna,model="JC69")
D_primates

## compute NJ tree
primates_nj<-NJ(D_primates)
primates_nj

## plot in unrooted style
plot(primates_nj,type="unrooted")

## root our tree and re-plot it
rooted.primates_nj<-root(primates_nj,outgroup="Lemur",
  resolve.root=TRUE)
plotTree(rooted.primates_nj)


## homework is going to be to program a bootstrapping analysis using NJ for the 
## primate DNA dataset

## hint number 1: how to create a SINGLE bootstrap resampled alignment
m<-ncol(primate_dna)
ii<-sample(1:m,m,replace=TRUE)
bs_dataset<-primate_dna[,ii]
bs_dataset

## comparison to our original data show that it has been resampled
primate_dna

## hint number 2: how to iterate across our bootstrap resamples
iter<-100
bs_datasets<-list()
for(i in 1:iter){
  ii<-sample(1:m,m,replace=TRUE)
  bs_datasets[[i]]<-primate_dna[,ii]
}
