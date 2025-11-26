## load phytools
library(phytools)

## read example files
## https://liamrevell.github.io/biol634/data/pagel94-example.csv
pagel94_data<-read.csv(
  file="https://liamrevell.github.io/biol634/data/pagel94-example.csv",
  row.names=1,stringsAsFactors=TRUE)
head(pagel94_data)
## https://liamrevell.github.io/biol634/data/pagel94-example.tre
pagel94_tree<-read.nexus(
  file="https://liamrevell.github.io/biol634/data/pagel94-example.tre")
plotTree(pagel94_tree,fsize=0.7)

## plot tree with data
plotTree.datamatrix(pagel94_tree,pagel94_data[,4:5])

## pull out x4 & x5
x4<-setNames(pagel94_data$x4,rownames(pagel94_data))
x5<-setNames(pagel94_data$x5,rownames(pagel94_data))

## they should look like factors
head(x4)
head(x5)

## fit Pagel '94 model
fit_x4x5<-fitPagel(pagel94_tree,x=x4,y=x5)
fit_x4x5

## let's plot the model
plot(fit_x4x5,lwd.by.rate=TRUE)

## now let's do characters x2 & x3
x2<-setNames(pagel94_data$x2,rownames(pagel94_data))
x3<-setNames(pagel94_data$x3,rownames(pagel94_data))

## fit the model
fit_x2x3<-fitPagel(pagel94_tree,x=x2,y=x3)
fit_x2x3

## plot our fitted model
plot(fit_x2x3,lwd.by.rate=TRUE)

## let's graph our data
dev.off()
plotTree.datamatrix(pagel94_tree,pagel94_data[,2:3])

## homework: analyze the two discrete characters in bonyfish.csv
## including graphing the data and the fitted model, and decide
## whether you think a hypothesis of correlated character evolution
## is warranted
## as a bonus, explore what the argument dep.var does in fitPagel

## let's do something else.
liolaemid_data<-read.csv(
  file="https://liamrevell.github.io/biol634/data/Liolaemidae.data.csv",
  row.names=1)
head(liolaemid_data)
liolaemid_tree<-read.nexus(
  file="https://liamrevell.github.io/biol634/data/Liolaemidae.MCC.nex")
liolaemid_tree

## run name.check
geiger::name.check(liolaemid_tree,liolaemid_data)

## pull out liolaemid parity mode
liol_parity<-setNames(liolaemid_data$parity_mode,
  rownames(liolaemid_data))
head(liol_parity)

## start by fitting a simple Mk (ARD) model
liol_mk<-fitMk(liolaemid_tree,liol_parity,
  model="ARD",pi="fitzjohn")
liol_mk

## let's fit a hidden-rates model to this
?fitHRM
parallel::detectCores()
liol_hrm_2<-fitHRM(liolaemid_tree,liol_parity,
  ncat=2,parallel=TRUE,niter=10,ncores=10)

## compare to our other fitted model
anova(liol_mk,liol_hrm_2)

## to fit our threshold model...
liol_thresh<-fitThresh(liolaemid_tree,
  liol_parity,parallel=TRUE)

## compare all three of our models
anova(liol_thresh,liol_mk,liol_hrm_2)
