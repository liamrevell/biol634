## load packages
library(phytools)

## load data
data(mammal.tree)
data(mammal.data)

## look at data
mammal.tree
plotTree(mammal.tree,fsize=0.8,ftype="i")
head(mammal.data)

## to test for an evolutionary correlation using
## contrasts we need to pull out our variables of
## interest
lnBodyMass<-setNames(log(mammal.data$bodyMass),
  rownames(mammal.data))
lnHomeRange<-setNames(log(mammal.data$homeRange),
  rownames(mammal.data))

## compute contrasts for each of my two traits
pic_bodymass<-pic(lnBodyMass,mammal.tree)
pic_homerange<-pic(lnHomeRange,mammal.tree)

## here are our contrasts
pic_bodymass

## add nodelabels
nodelabels(frame="circle",bg="white",cex=0.5)

## we're already ready to fit a regression model
## based on these contrasts
mammal_pic_regression<-lm(pic_homerange~pic_bodymass+0)
mammal_pic_regression

## let's plot our contrasts and fitted model
dev.off()
plot(pic_bodymass,pic_homerange,
  pch=16,bty="n",las=1,
  xlab="PICs for log(body mass)",
  ylab="PICs for log(home range)")
grid()
abline(mammal_pic_regression,lwd=2)

## load nlme package
library(nlme)

## we're going to build our phylogenetic error
## structure (our phylogenetic covariance matrix)
?corBrownian
mammal_spp<-rownames(mammal.data)
mammal_corBM<-corBrownian(phy=mammal.tree,
  form=~mammal_spp)
mammal_corBM

## now we're ready to fit our model!
?gls
mammal_pgls_regression<-gls(log(homeRange)~log(bodyMass),
  data=mammal.data,correlation=mammal_corBM)
mammal_pgls_regression

## now that we understand that PGLS and contrasts are
## equivalent in this simple case, let's explore situations
## in which contrasts cannot be used.
data(primate.tree)
data(primate.data)

## let's look at our primate
primate.data
head(primate.data)

## create our primate "corStruct" object for our primate
## tree
primates_spp<-rownames(primate.data)
primates_corBM<-corBrownian(phy=primate.tree,
  form=~primates_spp)
primates_corBM
colnames(primate.data)
primate_pgls_ancova<-gls(
  log(Orbit_area)~log(Skull_length)+Activity_pattern,
  data=primate.data,
  correlation=primates_corBM)
primate_pgls_ancova

## test for significance
anova(primate_pgls_ancova)

## run summary on the model
summary(primate_pgls_ancova)

## run our PGLS ANCOVA with an interaction term
primate_pgls_ancova.interaction<-gls(
  log(Orbit_area)~log(Skull_length)*Activity_pattern,
  data=primate.data,
  correlation=primates_corBM)
primate_pgls_ancova.interaction
