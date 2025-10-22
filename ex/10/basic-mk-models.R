## load libraries
library(phytools)
library(geiger)

## load pre-packaged data
data("primate.tree")
data("primate.data")

## let's get a sense of these data
head(primate.data)

## let's plot our tree
plotTree(primate.tree,fsize=0.7,ftype="i")

## our data
primate.data[,"Activity_pattern"]
primate.data$Activity_pattern
rownames(primate.data)

## let's get our discrete character

## one way
activity_pattern<-primate.data$Activity_pattern
names(activity_pattern)<-rownames(primate.data)
head(activity_pattern,10)

## using setNames
activity_pattern<-setNames(primate.data$Activity_pattern,
  rownames(primate.data))
head(activity_pattern,10)

## plot our tree & add an axis
plotTree(primate.tree,mar=c(3.1,1.1,1.1,1.1),
  fsize=0.7,ftype="i")
axis(1)

## fit simplest model: the "ER" model
primate_er<-fitMk(primate.tree,activity_pattern,
  model="ER")
primate_er

## how much evolution?
sum(primate.tree$edge.length)*0.002802

## plot our fitted model
plot(primate_er)

## fit most complex model: the "ARD" model
primate_ard<-fitMk(primate.tree,activity_pattern,
  model="ARD")
primate_ard

## plot this model too
plot(primate_ard,width=TRUE,color=TRUE)

## fit an intermediate model: the "SYM" model
primate_sym<-fitMk(primate.tree,activity_pattern,
  model="SYM")
primate_sym
plot(primate_sym,width=TRUE,color=TRUE)

## statistically compare model
anova(primate_er,primate_sym,primate_ard)
