
#-----------------------------------------------#
#####Transmission tree doe#####

#by Nick Fountain-Jones, November 2018
#-----------------------------------------------#

#install if necesary

#library(devtools)

#devtools::install_github('xavierdidelot/TransPhylo')

#-----------------------------------------------#
#####Prelims#####
#-----------------------------------------------#

library(TransPhylo)
library(ape)

rm(list=ls())
set.seed(0)

#-----------------------------------------------#
#####Load in the phylogenetic tree#####
#-----------------------------------------------#

#has to be a newick file
ptree<-ptreeFromPhylo(read.tree('FIV_WS_A_brownian_MCCb'),dateLastSample=2014.2) #read.tree is also an option.
str(ptree)


#IMPORTANT callibration step - need a reasonable distribution of generation time. We need to discuss this with team virus.


w.shape=1 # w.shape=1,w.scale=0.5 = 6 months on average between transmission events, 
w.scale=0.5
#scale=1,w.shape=1 we get an Exponential distribution with mean 1 year. 

startPi <- 0.95 # proportion sampled

dateT=2014.9 #end of sampling date

#-----------------------------------------------#
#####Calculate transmission tree#####
#-----------------------------------------------#

#update pi
recordWSA<-inferTTree(ptree,mcmcIterations=10,w.shape=w.shape,w.scale=w.scale,dateT=dateT, startPi = startPi, updatePi = F)

lastIteration<-recordWSA[[length(recordWSA)]]
plotCTree(lastIteration$ctree)

#check out mcmc diagnostics. WIll need to add some more

par(mfrow=c(2,2))
plot(sapply(recordWSA,function(x) x$pTTree+x$pPTree),ylab='Posterior probability',
     xlab='MCMC iterations',type='l')
plot(sapply(recordWSA,function(x) x$pi),ylab='Sampling proportion pi',
     xlab='MCMC iterations',type='l')
plot(sapply(recordWSA,function(x) x$neg),ylab='Within-host coalescent rate Ne*g',
     xlab='MCMC iterations',type='l')
plot(sapply(recordWSA,function(x) x$off.r),ylab='Basic reproduction number R',
     xlab='MCMC iterations',type='l')

cons=consTTree(recordWSA);str(cons)
dev.off()

plotTTree2(cons,w.shape,w.scale)
plotTTree(cons,w.shape,w.scale)

tre <-computeMatWIW(recordWSA, burnin = 0.1);str(tre)

#calculate R2

R2 <- sapply(recordWSA,function(x) x$off.r)

summary(R2)

#-----------------------------------------------#
#####Turn into a graph#####
#-----------------------------------------------#

library(igraph)

G <- as.directed(graph.adjacency(tre, weighted = T));str(G)
#calculate degree
deg <- degree(G, mode="all")

V(net)$size <- deg*3
plot(G, edge.color=c("red","green")[sign(E(G)$weight)], 
     edge.width = 10 *abs(E(G)$weight), vertex.size=1+(5*deg), vertex.color='gold', edge.curved=0.6, edge.arrow.size=10 *abs(E(G)$weight), vertex.label.dist = 2, vertex.frame.color = 'gray', vertex.label.color='black')


