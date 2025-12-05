## today we're learning about how to customize our own 
## plotted phylogenies in R

## start by downloading a file
download.file(
  url="https://liamrevell.github.io/biol634/data/vertebrate_phylogeny.tre",
  destfile="vertebrate_phylogeny.tre")

## let's get our phylogeny into R
library(phytools)
vertebrate_tree<-read.tree(file="vertebrate_phylogeny.tre")

## this is our tree
vertebrate_tree

## let's see how this tree is plotted in R!
## we learn that's it plotted in a standard R graphical
## device
plotTree(vertebrate_tree,
  mar=c(3,3,1,1))
args(axis)
axis(side=1)
axis(2)
grid()

## compute the total height of our tree
nodeHeights(vertebrate_tree)
max(nodeHeights(vertebrate_tree))

## let's add the simplest of graphical features!
abline(v=100,col="blue",lwd=4)

## let's replot our tree
plotTree(vertebrate_tree,
  mar=c(4,1,1,1),ftype="i")

## let's add an axis that runs backwards in time
## from the present day
args(axis)
h<-max(nodeHeights(vertebrate_tree))
h
seq(0,450,by=50)
h-seq(0,450,by=50)
axis(side=1,at=h-seq(0,450,by=50),
  labels=seq(0,450,by=50))

## let's add a line 65 mybp
abline(v=h-65,lty="dotted",lwd=2,col="red")

## let's add an arrow
args(arrows)
arrows(x0=h-125,y0=8,x1=h-65,y1=9.5,length=0.1,
  lwd=3,col="blue")
text(x=h-125,y=8,"end of Cretaceous",pos=2)

## let's go back & plot our tree again.
plotTree(vertebrate_tree,mar=c(4,4,1,1),ftype="i")
axis(1)
axis(2)
grid()

## call the function nodelabels
nodelabels()

## don't do this at home! (but as a demo, I showed that
## node labels will still plot even after the "phylo"
## object has been deleted)
rm(list=ls())

nodelabels()

dev.off()
plot(NA,xlim=c(0,473),ylim=c(1,11))
nodelabels()

## let's go back & plot our tree again.
plotTree(vertebrate_tree,mar=c(4,4,1,1),ftype="i")
axis(1)
axis(2)
grid()

## we're going to see how a special object, "last_plot.phylo"
## has been created in a parallel "hidden" R environment
## this environment is called .PlotPhyloEnv
pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
pp

## overlay plots
points(pp$xx,pp$yy,pch=21,bg="grey",cex=1.5)

## copy something from another namespace of phytools in
## GlobalEnv
copy_of_plotTree<-phytools::plotTree
copy_of_plotTree(vertebrate_tree)

## back to the phylogeny stuff
plotTree(vertebrate_tree,mar=c(4,4,1,1),ftype="i")
axis(1)
axis(2)
grid()
pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)

## get the common ancestor of mammals
mammal_node<-getMRCA(vertebrate_tree,
  c("Homo_sapiens","Bos_taurus"))
mammal_node

## add our point
points(pp$xx[mammal_node],pp$yy[mammal_node],
  pch=21,bg="blue",cex=2)

## get all the nodes that descend from the common
## ancestor of all mammals.
dd<-getDescendants(vertebrate_tree,node=mammal_node)
dd
points(pp$xx[dd],pp$yy[dd],pch=16,col="grey")

## let's figure out the limits of a clade box
min_x<-pp$xx[mammal_node]-10
max_x<-pp$x.lim[2] ## or max(pp$xx[dd])
min_y<-min(pp$yy[dd])-0.5
max_y<-max(pp$yy[dd])+0.5

## add the clade box
polygon(
  x=c(min_x,max_x,max_x,min_x),
  y=c(min_y,min_y,max_y,max_y),
  col="lightblue",border=FALSE)
## the only problem is that it's superimposed on top
## of the tree

## to fix this, let's graph our phylogeny without plotting it
plotTree(vertebrate_tree,mar=c(4,4,1,1),ftype="i",
  plot=FALSE)
axis(1)
axis(2)
grid()
pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)

## get common ancestor
dd<-getDescendants(vertebrate_tree,node=mammal_node)

## set limits of clade box
min_x<-pp$xx[mammal_node]-10
max_x<-pp$x.lim[2] ## or max(pp$xx[dd])
min_y<-min(pp$yy[dd])-0.5
max_y<-max(pp$yy[dd])+0.5

## graph the clade box
polygon(
  x=c(min_x,max_x,max_x,min_x),
  y=c(min_y,min_y,max_y,max_y),
  col="lightblue",border=FALSE)

## add our tree
plotTree(vertebrate_tree,mar=c(4,4,1,1),ftype="i",
  plot=TRUE,add=TRUE)

## how about we try to create a custom phylogeny plotter
## using what we've learned?

## first plot the tree without plotting
plotTree(vertebrate_tree,plot=FALSE)

## get the coordinates
pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(pp$xx,pp$yy)

## add the edges
for(i in 1:nrow(pp$edge)){
  lines(x=pp$xx[pp$edge[i,]],
    y=pp$yy[pp$edge[i,]])
}

## add the labels
for(i in 1:Ntip(vertebrate_tree)){
  text(pp$xx[i],pp$yy[i],
    vertebrate_tree$tip.label[i],pos=4)
}

## put it all together as function
lior_tree_plotter<-function(phy){
  plotTree(phy,plot=FALSE)
  pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  points(pp$xx,pp$yy)
  for(i in 1:nrow(pp$edge)){
    lines(x=pp$xx[pp$edge[i,]],
      y=pp$yy[pp$edge[i,]])
  }
  for(i in 1:Ntip(phy)){
    text(pp$xx[i],pp$yy[i],
      phy$tip.label[i],pos=4)
  }
  cat("Did this work?\n") ## print some text
}

## test it
lior_tree_plotter(vertebrate_tree)

## for a truer test, let's load a different phylogeny
data("salamanders")
lior_tree_plotter(salamanders)

## other functions that work with phylogenies
## still work here!
nodelabels()
