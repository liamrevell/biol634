## this is the code for both my solution to the homework & what we did
## together today in class

## first, the homework solution

## homework details: https://liamrevell.github.io/biol634/hw/4/hw4-problem.html

## load packages
library(phytools)
library(geiger)

## in class we loaded & cleaned up the data as follows

## load data
squamate_data<-read.csv(
  file="https://liamrevell.github.io/biol634/data/brandley_table.csv",
  row.names=1)

## load tree
squamate_tree<-read.nexus(
  file="https://liamrevell.github.io/biol634/data/squamate.tre")

## gsub " " for "_" in species labels of data
rownames(squamate_data)<-gsub(" ","_",rownames(squamate_data))

## still a name mismatch
name.check(squamate_tree,squamate_data)

## subsample squamate_data
squamate_data<-squamate_data[squamate_tree$tip.label,]

## pull out hind digits & round
hind_digits<-setNames(squamate_data$Toes,
  rownames(squamate_data))
hind_digits<-round(hind_digits)

## fit ER model
squamate_er<-fitMk(squamate_tree,hind_digits,model="ER",
  pi="fitzjohn")

## fit 1-rate loss-only model
D<-matrix(c(
  0,0,0,0,0,0,
  1,0,0,0,0,0,
  0,1,0,0,0,0,
  0,0,1,0,0,0,
  0,0,0,1,0,0,
  0,0,0,0,1,0),
  6,6,byrow=TRUE,
  dimnames=list(0:5,0:5))
squamate_lossonly_1r<-fitMk(squamate_tree,hind_digits,
  model=D,pi="fitzjohn")

## fit loss-only multi-rate model
D<-matrix(c(
  0,0,0,0,0,0,
  1,0,0,0,0,0,
  0,2,0,0,0,0,
  0,0,3,0,0,0,
  0,0,0,4,0,0,
  0,0,0,0,5,0),
  6,6,byrow=TRUE,
  dimnames=list(0:5,0:5))
squamate_lossonly_5r<-fitMk(squamate_tree,hind_digits,
  model=D,pi="fitzjohn")

## fit ordered + jump model
D<-matrix(c(
  0,0,0,0,0,0,
  1,0,0,0,0,0,
  6,2,0,0,0,0,
  6,0,3,0,0,0,
  6,0,0,4,0,0,
  6,0,0,0,5,0),
  6,6,byrow=TRUE,
  dimnames=list(0:5,0:5))
squamate_lossonly_jump<-fitMk(squamate_tree,hind_digits,
  model=D,pi="fitzjohn")

## fit loss/gain two-rate model
D<-matrix(c(
  0,2,0,0,0,0,
  1,0,2,0,0,0,
  0,1,0,2,0,0,
  0,0,1,0,2,0,
  0,0,0,1,0,2,
  0,0,0,0,1,0),
  6,6,byrow=TRUE,
  dimnames=list(0:5,0:5))
squamate_lossgain_2r<-fitMk(squamate_tree,hind_digits,
  model=D,pi="fitzjohn")

## fit loss/gain multi-rate model
D<-matrix(c(
  0,6,0,0,0,0,
  1,0,7,0,0,0,
  0,2,0,8,0,0,
  0,0,3,0,9,0,
  0,0,0,4,0,10,
  0,0,0,0,5,0),
  6,6,byrow=TRUE,
  dimnames=list(0:5,0:5))
squamate_lossgain_10r<-fitMk(squamate_tree,hind_digits,
  model=D,pi="fitzjohn")

## fit loss/gain + jump model
D<-matrix(c(
  0,7,0,0,0,0,
  1,0,8,0,0,0,
  6,2,0,9,0,0,
  6,0,3,0,10,0,
  6,0,0,4,0,11,
  6,0,0,0,5,0),
  6,6,byrow=TRUE,
  dimnames=list(0:5,0:5))
squamate_lossgain_jump<-fitMk(squamate_tree,hind_digits,
  model=D,pi="fitzjohn")

## compare among all fitted models (sorted by complexity)
anova(
  squamate_er,
  squamate_lossonly_1r,
  squamate_lossgain_2r,
  squamate_lossonly_5r,
  squamate_lossonly_jump,
  squamate_lossgain_10r,
  squamate_lossgain_jump)

## from this we see evidence that we FAILED to converge to the
## true MLEs for our two most complex models. Why do I think that?

## the answer is that in all cases when we compare two NESTED models
## (in which the more complex model has the simpler model as a special
## case) the more complex model has to have an equal or higher 
## log-likelihood

## let's run multiple optimization iterations in parallel to see if
## we can find the true MLEs for our two models that indicate 
## failure to converge

## load required packages
library(foreach)
library(doParallel)

niter<-30 ## number of iterations
ncores<-min(niter,detectCores()-4) ## number of cores (leave 4 free)
mc<-makeCluster(ncores,type="PSOCK") ## make cluster
registerDoParallel(cl=mc) ## register cluster

## loss/gain multi-rate model design matrix
D<-matrix(c(
  0,6,0,0,0,0,
  1,0,7,0,0,0,
  0,2,0,8,0,0,
  0,0,3,0,9,0,
  0,0,0,4,0,10,
  0,0,0,0,5,0),
  6,6,byrow=TRUE,
  dimnames=list(0:5,0:5))

## run 10 optimization iterations across cores
lossgain10r_fits<-foreach(i=1:niter)%dopar%{
  phytools::fitMk(squamate_tree,hind_digits,
    model=D,pi="fitzjohn",rand_start=TRUE,
    lik.func="lik",logscale=sample(c(TRUE,FALSE),1))
}

## get log-likelihoods
lnL<-sapply(lossgain10r_fits,logLik)
lnL

## new best model
squamate_lossgain_10r<-lossgain10r_fits[[which.max(lnL)]]

## loss/gain + jump design matrix
D<-matrix(c(
  0,7,0,0,0,0,
  1,0,8,0,0,0,
  6,2,0,9,0,0,
  6,0,3,0,10,0,
  6,0,0,4,0,11,
  6,0,0,0,5,0),
  6,6,byrow=TRUE,
  dimnames=list(0:5,0:5))

## run 10 optimization iterations across cores
lossgainjump_fits<-foreach(i=1:niter)%dopar%{
  phytools::fitMk(squamate_tree,hind_digits,
    model=D,pi="fitzjohn",rand_start=TRUE,
    lik.func="lik",logscale=sample(c(TRUE,FALSE),1))
}

## get log-likelihoods
lnL<-sapply(lossgainjump_fits,logLik)
lnL

## new best model
squamate_lossgain_jump<-lossgainjump_fits[[which.max(lnL)]]

## stop cluster
stopCluster(cl=mc)

## re-run comparison
anova(
  squamate_er,
  squamate_lossonly_1r,
  squamate_lossgain_2r,
  squamate_lossonly_5r,
  squamate_lossonly_jump,
  squamate_lossgain_10r,
  squamate_lossgain_jump)

## now, what we did exactly in class

## load packages
library(phytools)
library(geiger)

## read data
squamate_data<-read.csv(
  file="https://liamrevell.github.io/biol634/data/brandley_table.csv",
  row.names=1)

## read tree
squamate_tree<-read.nexus(
  file="https://liamrevell.github.io/biol634/data/squamate.tre")

## address mismatch problem
rownames(squamate_data)<-gsub(" ","_",rownames(squamate_data))

## run name.check
geiger::name.check(squamate_tree,squamate_data)

## subsample my data frame
squamate_data<-squamate_data[squamate_tree$tip.label,]

## re-run name.check (should pass)
name.check(squamate_tree,squamate_data)

## extract hind digit number
hind_digits<-setNames(round(squamate_data$Toes),
  rownames(squamate_data))
head(hind_digits)

## first thing: let's compare to geiger::fitDiscrete

## fit using phytools first
squamate_er<-fitMk(squamate_tree,hind_digits,
  model="ER",pi="fitzjohn")
squamate_er

## convert to factor because geiger doesn't like 0s
as.factor(hind_digits)->hind_digits

## fit same model using geiger::fitDiscrete
squamate_er.geiger<-fitDiscrete(squamate_tree,
  hind_digits,model="ER")
squamate_er.geiger

## custom design matrix (simple loss-only one-rate model)
D<-matrix(c(
  0,0,0,0,0,0,
  1,0,0,0,0,0,
  0,1,0,0,0,0,
  0,0,1,0,0,0,
  0,0,0,1,0,0,
  0,0,0,0,1,0),
  6,6,byrow=TRUE,
  dimnames=list(levels(hind_digits),levels(hind_digits)))
D

## fit using phytools
squamate_lossonly_1r<-fitMk(squamate_tree,
  hind_digits,model=D,pi="fitzjohn")

## same model using geiger
squamate_lossonly_1r.geiger<-fitDiscrete(squamate_tree,
  hind_digits,model=D)

## graph fitted model
plot(squamate_lossonly_1r,show.zeros=FALSE)

## for some reason we seem to be getting non-convergence with
## geiger (not sure what's going on)
squamate_lossonly_1r
squamate_lossonly_1r.geiger

## all the rates are effectively zero
plot(squamate_lossonly_1r.geiger,show.zeros=FALSE)

## design matrix for custom gain/loss model
D<-matrix(c(
  0,2,0,0,0,0,
  1,0,2,0,0,0,
  0,1,0,2,0,0,
  0,0,1,0,2,0,
  0,0,0,1,0,2,
  0,0,0,0,1,0),
  6,6,byrow=TRUE,
  dimnames=list(levels(hind_digits),levels(hind_digits)))
D

## fit using phytools and geiger
squamate_lossgain_2r<-fitMk(squamate_tree,hind_digits,
  model=D,pi="fitzjohn")
squamate_lossgain_2r

squamate_lossgain_2r.geiger<-fitDiscrete(squamate_tree,
  hind_digits,model=D)
squamate_lossgain_2r.geiger

## now our results look to be identical

## design matrix for multi-rate loss/gain model
D<-matrix(c(
  0,10,0,0,0,0,
  1,0,9,0,0,0,
  0,2,0,8,0,0,
  0,0,3,0,7,0,
  0,0,0,4,0,6,
  0,0,0,0,5,0),
  6,6,byrow=TRUE,
  dimnames=list(levels(hind_digits),levels(hind_digits)))
D

## design matrix for loss/gain + jump
D<-matrix(c(
  0,10,0,0,0,0,
  1,0,9,0,0,0,
  11,2,0,8,0,0,
  11,0,3,0,7,0,
  11,0,0,4,0,6,
  11,0,0,0,5,0),
  6,6,byrow=TRUE,
  dimnames=list(levels(hind_digits),levels(hind_digits)))
D

## log-likelihoods (as obtained by students)
## loss-only no jump: -211.8014
## loss-only + jump : -207.5148
## loss/gain no jump: -220.3415
## loss/gain + jump : -209.6297
## ARD              : -192.4331

## what's wrong here (nested models don't necessarily increase in
## log-likelihood)

## let's load my own solution from earlier to compare
load("hw4-solution.RData")

anova(
  squamate_er,
  squamate_lossonly_1r,
  squamate_lossgain_2r,
  squamate_lossonly_5r,
  squamate_lossonly_jump,
  squamate_lossgain_10r,
  squamate_lossgain_jump)

## these findings underscore the importance of running multiple
## optimization iterations for difficult models

## how do we do this using phytools::fitMk?

## (1.) using a for loop
er_fits<-list()
niter<-10

for(i in 1:niter){
  cat(paste("running optimization iteration",i,"\n"))
  er_fits[[i]]<-fitMk(squamate_tree,hind_digits,
    model="ER",pi="fitzjohn",rand_start=TRUE)
}

er_fits ## all of these are the same in this simple case

## (2.) using parallelization and foreach

## "Danny's model" of three-rate loss-only, plus jump losses
D<-matrix(c(
  0,0,0,0,0,0,
  1,0,0,0,0,0,
  4,2,0,0,0,0,
  4,0,2,0,0,0,
  4,0,0,2,0,0,
  4,0,0,0,3,0),
  6,6,byrow=TRUE,
  dimnames=list(levels(hind_digits),levels(hind_digits)))
D

## load packages
library(foreach)
library(doParallel)

## set number of iterations and number of cores
niter<-10
detectCores()

ncores<-min(niter,detectCores()-2) ## two fewer than ALL cores
ncores

## make cluster
mc<-makeCluster(ncores,type="PSOCK")

## register cluster
registerDoParallel(cl=mc)

## obtain fits using foreach
dannys_model_fits<-foreach(i=1:niter)%dopar%{
  phytools::fitMk(squamate_tree,hind_digits,
    model=D,pi="fitzjohn",rand_start=TRUE,
    lik.func="lik",logscale=sample(c(TRUE,FALSE),1))
}

## here is all fits
dannys_model_fits

## don't forget to stop your cluster!
stopCluster(cl=mc)

## get all the log-likelihoods
lnL<-sapply(dannys_model_fits,logLik)
lnL

## which one is the highest?
which.max(lnL)

## pull out best-fitting model
danny_model<-dannys_model_fits[[which.max(lnL)]]
danny_model

## graph best-fitting model
plot(danny_model,width=TRUE,color=TRUE)
