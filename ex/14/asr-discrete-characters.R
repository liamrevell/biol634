## let's do ancestral state reconstruction of discrete characters!!

## load phytools
library(phytools)

## we can start by keeping it simple:
data(primate.tree)
data(primate.data)

## look at data
head(primate.data)

## pull out activity pattern.
activity<-setNames(primate.data$Activity_pattern,
  rownames(primate.data))
head(activity)

## fit models:
primate_er<-fitMk(primate.tree,activity,
  model="ER")
primate_sym<-fitMk(primate.tree,activity,
  model="SYM")
primate_ard<-fitMk(primate.tree,activity,
  model="ARD")

## look at ARD model (for instance)
primate_ard
plot(primate_ard,width=TRUE,color=TRUE,
  xlim=c(-2,1),ylim=c(-1,1))

## we would normally do model comparison
anova(primate_er,
  primate_sym,
  primate_ard)

## the "best" model is the ER model
primate_er.marginal<-ancr(primate_er)

## plot the result
plot(primate_er.marginal,
  args.plotTree=list(type="arc",fsize=0.6,arc_height=0.2,part=1),
  args.nodelabels=list(piecol=c("purple","orange","black")))

## let's get model-averaged marginal states
anova(primate_er,
  primate_sym,
  primate_ard)->primate_aov
primate_aov
primate_model_averaged.marginal<-ancr(primate_aov)

## plot model averaged result
plot(primate_model_averaged.marginal,
  args.plotTree=list(type="arc",fsize=0.6,arc_height=0.2,part=1),
  args.nodelabels=list(piecol=c("purple","orange","black")))

## let's do joint reconstruction for the same data (to compare)
## we're just going to use the "best-supported" model, the ER model
primate_er.joint<-ancr(primate_er,type="joint")

## let's plot our result
plot(primate_er.joint)

## let's compare our marginal and joint reconstructions
par(mfcol=c(1,2))
plot(primate_er.marginal)
plot(primate_er.joint)

## let's do stochastic character mapping!!
?simmap
primate_simmap<-simmap(primate_aov)
primate_simmap
class(primate_simmap)

## plot our trees
par(mfrow=c(10,10))
plot(primate_simmap,lwd=1,ftype="off")

dev.off()
plot(primate_simmap[[4]],
  colors=setNames(c("purple","orange","black"),
    levels(activity)),fsize=0.7,ftype="i")

## let's get the density of changes of different types
primate_simmap.dd<-density(primate_simmap)
primate_simmap.dd

## let's graph this
dev.off()
plot(primate_simmap.dd)
