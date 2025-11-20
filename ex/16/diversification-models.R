## load phytools
library(phytools)

## let's "prove" the pull-of-the-present
?pbtree

## we can figure out what net-diversification rate I need
## to get (say) 1,000 taxa after 100 ma.
## we need to solve 1000 = 2*exp(nd*100)
## ln(1000) = ln(2) + nd*100
## (ln(1000) - ln(2))/100 = nd
nd<-(log(1000)-log(2))/100
lambda<-0.1
mu<-lambda-nd
lambda
mu

## now let's simulate a tree
tree<-pbtree(b=lambda,d=mu,t=100)

## plot this tree
plotTree(tree,lwd=1,ftype="off")

## generate a lineage through time plot
object<-ltt(tree,plot=FALSE)
dev.off()
plot(object)

## prune out all the extinct tips
?getExtinct
extantonly_tree<-drop.tip(tree,getExtinct(tree))
extantonly_tree

## plot my "reconstructed tree"
plotTree(extantonly_tree,ftype="off",lwd=1)

## recalculate my LTT plot
object<-ltt(extantonly_tree,plot=FALSE)
dev.off()
plot(object)

## create a pure-birth tree of the same total 
## length & diversity
purebirth_tree<-pbtree(n=Ntip(extantonly_tree),
  scale=max(nodeHeights(extantonly_tree)))

## generate an LTT plot for this
purebirth_tree.ltt<-ltt(purebirth_tree)

## we're going to a use a very simple phytools
## function to fit the birth-death model
bd_model<-fit.bd(extantonly_tree)
bd_model

## apply the same thing to our pure-birth tree
bd_model.purebirth<-fit.bd(purebirth_tree)
bd_model.purebirth

## we can test a hypothesis about this model
yule_model<-fit.yule(extantonly_tree)
yule_model

## compare models
anova(yule_model,bd_model)

## inspect the object from ltt
object

## read tree from website
elapidae_tree<-read.nexus(
  file="https://liamrevell.github.io/biol634/data/elapidae.nex")
elapidae_tree

## approximate number of total elapids
nlow<-300
nhigh<-360

## let's fit a birth-death model accounting for sampling fraction
elapid_bd.low<-fit.bd(elapidae_tree,
  rho=Ntip(elapidae_tree)/nlow)
elapid_bd.low
elapid_bd.high<-fit.bd(elapidae_tree,
  rho=Ntip(elapidae_tree)/nhigh)
elapid_bd.high

## we discovered that, though the model of lambda + mu + rho is
## theoretically non-identifiable, it becomes identifiable in this
## case because mu hits its lower bound of 0

## we can fit the same model using diversitree
library(diversitree)
bd_lik.low<-make.bd(elapidae_tree,
  sampling.f=Ntip(elapidae_tree)/nlow)
bd_lik.low
mle.bd_low<-find.mle(bd_lik.low,x.init=c(0.2,0.05))
mle.bd_low

bd_lik.high<-make.bd(elapidae_tree,
  sampling.f=Ntip(elapidae_tree)/nhigh)
bd_lik.high
mle.bd_high<-find.mle(bd_lik.high,x.init=c(0.2,0.05))
mle.bd_high
