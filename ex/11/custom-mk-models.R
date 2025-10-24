## we can start by using download.file to pull down the data files
## we're going to use in this example
?download.file

## this is url
## https://liamrevell.github.io/biol634/data/brandley_table.csv

## download file to local destination file
download.file(
  url="https://liamrevell.github.io/biol634/data/brandley_table.csv",
  destfile="brandley_table.csv")

## read file
squamate_data<-read.csv(file="brandley_table.csv",
  row.names=1)
head(squamate_data)

## alternatively, we could've read directly from the web
squamate_data<-read.csv(
  file="https://liamrevell.github.io/biol634/data/brandley_table.csv",
  row.names=1)
head(squamate_data)

## download squamate tree
## https://liamrevell.github.io/biol634/data/squamate.tre
download.file(
  url="https://liamrevell.github.io/biol634/data/squamate.tre",
  destfile="squamate.tre")

## load our libraries
library(phytools)
library(geiger)

## read in our squamate tree using ape::read.nexus
squamate_tree<-read.nexus(file="squamate.tre")
squamate_tree

## take a look at our data and tree
head(squamate_data)
squamate_tree

## function in geiger to check names
?name.check
chk<-name.check(squamate_tree,squamate_data)
summary(chk)

## this tells us that there are mismatches between our tree
## and data. Specifically, labels are "Genus species" in data but
## "Genus_species" in tree

## fix our input data
spp<-rownames(squamate_data)
spp
?gsub
spp<-gsub(pattern=" ",replacement="_",x=spp)
spp
rownames(squamate_data)<-spp
head(squamate_data)

## we could've done this all in one step as follows
rownames(squamate_data)<-gsub(" ","_",rownames(squamate_data))

## check our data
head(squamate_data)

## check for alignment of data & tree
name.check(squamate_tree,squamate_data)

## we found a mismatch, but in this case there are just three
## taxa in our data but not the tree

## let's pull out hind digits
hind_digits<-setNames(squamate_data$Toes,
  rownames(squamate_data))
hind_digits

## let's round our character
hind_digits<-round(hind_digits)
hind_digits

## subsample my data to include only species present
## in the tree. This remedies the problem we identified using 
## geiger::name.check
hind_digits<-hind_digits[squamate_tree$tip.label]

## let's start by using fitMk to fit an ER model
squamate_er<-fitMk(squamate_tree,hind_digits,model="ER",
  pi="fitzjohn")
squamate_er
plot(squamate_er)

## this is not very biologically realistic. To fit more biologically
## realistic models, we need to define their structure via custom
## design matrices

## let's now create a design matrix
## to start, let's imagine a 1-rate loss-only model
D<-matrix(0,6,6,dimnames=list(0:5,0:5))
D
D[2,1]<-1
D
D[3,2]<-1
D
D[4,3]<-1
D
D[5,4]<-1
D
D[6,5]<-1
D

## here's another way to build the same matrix in just one step
## (this is what I normally do)
D<-matrix(c(
  0,0,0,0,0,0,
  1,0,0,0,0,0,
  0,1,0,0,0,0,
  0,0,1,0,0,0,
  0,0,0,1,0,0,
  0,0,0,0,1,0),
  6,6,byrow=TRUE,
  dimnames=list(0:5,0:5))
D

## fit the 1-rate loss-only model
squamate_lossonly_1r<-fitMk(squamate_tree,hind_digits,
  model=D,pi="fitzjohn")
squamate_lossonly_1r

## plot this model
plot(squamate_lossonly_1r,show.zeros=FALSE)

## compare to ER model
anova(squamate_er,squamate_lossonly_1r)

## now let's make loss-only multi-rate model
D<-matrix(c(
  0,0,0,0,0,0,
  1,0,0,0,0,0,
  0,2,0,0,0,0,
  0,0,3,0,0,0,
  0,0,0,4,0,0,
  0,0,0,0,5,0),
  6,6,byrow=TRUE,
  dimnames=list(0:5,0:5))
D

## alternatively
D<-matrix(0,6,6,dimnames=list(0:5,0:5))
D
D[2,1]<-1
D
D[3,2]<-2
D
D[4,3]<-3
D
D[5,4]<-4
D
D[6,5]<-5
D

## let's fit it!!
squamate_lossonly_5r<-fitMk(squamate_tree,hind_digits,
  model=D,pi="fitzjohn")
squamate_lossonly_5r

## make a model comparison
anova(squamate_er,squamate_lossonly_1r,
  squamate_lossonly_5r)

## plot model
plot(squamate_lossonly_5r,width=TRUE,color=TRUE,
  show.zeros=FALSE)

## now let's consider a two process model in which there is both
## "jump" (e.g., developmental switch) evolution to "0" and 
## sequential ordered loss
## two process model
D[3,1]<-6
D[4,1]<-6
D[5,1]<-6
D[6,1]<-6
D

## equivalently
D<-matrix(c(
  0,0,0,0,0,0,
  1,0,0,0,0,0,
  6,2,0,0,0,0,
  6,0,3,0,0,0,
  6,0,0,4,0,0,
  6,0,0,0,5,0),
  6,6,byrow=TRUE,
  dimnames=list(0:5,0:5))

## fit this model
squamate_lossonly_jump<-fitMk(squamate_tree,hind_digits,
  model=D,pi="fitzjohn")

## compare all four models
anova(squamate_er,squamate_lossonly_1r,
  squamate_lossonly_5r,squamate_lossonly_jump)

## for homework try to fit the following additional models
## 1. A loss/gain two rate model (one rate of loss, one of gain)
## 2. A loss/gain multi-rate model (5 rates of loss, 5 rates of gain)
## 3. A loss/gain + "jump" model (5 rates of loss, 5 rates of gain, 
##    1 jump loss rate)
## 4. If you get bored, just add in the ARD & SYM models
## Make a comparison among all of them.

