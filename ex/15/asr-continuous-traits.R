## demo on the covariances of related taxa.

library(phytools)

tree<-read.tree(
  text="((A:0.1,B:0.1):0.6,C:0.7);"
)
plotTree(tree)

X<-fastBM(tree,nsim=1)
X
X<-fastBM(tree,nsim=1000,internal=TRUE)
X
X<-t(X)
X<-X[,-4]
head(X)

pairs(X)

## ancestral state reconstruction of continuous traits.

library(phytools)

## load Liolaemus data.

data("liolaemid.data")
data("liolaemid.tree")

## pull out temperature data.

head(liolaemid.data)
liol_temp<-setNames(liolaemid.data$temperature,
  rownames(liolaemid.data))

## check for mismatches (there is one mismatch)
geiger::name.check(data=liolaemid.data,phy=liolaemid.tree)

## we are going to reconstruct using phytools::fastAnc
?fastAnc
liol_temp<-liol_temp[liolaemid.tree$tip.label]
liol_temp.asr<-fastAnc(liolaemid.tree,liol_temp,
  CI=TRUE)

## here is our reconstruction
liol_temp.asr
print(liol_temp.asr,printlen=10)

## compare CI at root to the range of the data
range(liol_temp)

## let's plot our tree
plotTree(liolaemid.tree,ftype="off",lwd=1)

## let's do a contMap
liol_temp.cMap<-contMap(liolaemid.tree,liol_temp,plot=FALSE)

## recolor this & regraph it
liol_temp.cMap<-setMap(liol_temp.cMap,c("darkgreen","yellow","red"))
plot(liol_temp.cMap,lwd=c(3,10),outline=FALSE,fsize=c(0.4,1),
  ftype="off")

