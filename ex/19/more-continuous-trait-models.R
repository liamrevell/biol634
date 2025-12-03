## load phytools
library(phytools)

## load some data that comes with phytools
data("liolaemid.data")
data("liolaemid.tree")

## we've seen a variant of this data already
head(liolaemid.data)

## pull out temperature
liol_temp<-setNames(liolaemid.data$temperature,
  rownames(liolaemid.data))
head(liol_temp)

## prune our data
liol_temp<-liol_temp[
  -which(names(liol_temp)=="Ctenoblepharys_adspersa")]

## fit continuous character models using geiger::fitContinuous
library(geiger)
bm.liol_temp<-fitContinuous(liolaemid.tree,liol_temp,
  model="BM")
bm.liol_temp

## let's compete BM against EB 
eb.liol_temp<-fitContinuous(liolaemid.tree,liol_temp,
  model="EB")
eb.liol_temp

## new method to test the hypothesis that the rate of evolution
## of a discrete trait depends on a continuous character
?fitcontMk
parity_mode<-setNames(liolaemid.data$parity_mode,
  rownames(liolaemid.data))
head(parity_mode)

## drop mismatch taxon from data
parity_mode<-parity_mode[-1]
head(parity_mode)

## fit the model using fitcontMk
liol.contMk_fit<-fitcontMk(liolaemid.tree,parity_mode,
  liol_temp,parallel=TRUE,levs=50,ncores=10)

## graph fitted model
## (shows that that the rate of O->V is high for low environmental
## temperatures, which is consstent with our a priori hypothesis)
plot(liol.contMk_fit)
