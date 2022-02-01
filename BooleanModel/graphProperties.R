#Clear environment
rm(list=ls())
gc(reset = TRUE)
###FT-GRN Graph propierties, plots Supplementary Figure 2
library(BoolNet)
library(igraph)

#### Load Network, model logic rules
net <- loadNetwork("ftgrn.txt")

net <- plotNetworkWiring(net, plotIt = FALSE)


###Degree
grado<- degree(net)
####Mean degree
mean(grado)

##Outdegree
gradoS <- degree(net, mode = "out")
mean(gradoS)

##Indegree
gradoE <- degree(net, mode = "in")
mean(gradoE)
##Plot indegree and outdegree
pos_vector <- rep(3, length(gradoE))
pos_vector[names(gradoE) %in% c("CO", "FT", "AP1")] <- c(2,2,4)

plot(gradoS ~gradoE, cex= 0.8, xlab="Indegree", ylab = "Outdegree", type="p")
text(gradoS ~gradoE, labels= names(gradoE), cex=0.6, pos=pos_vector)
axis(1, at = c(1:14))
axis(2, at= c(1:12))

#Outdegree
z_ind<- scale(gradoS)
pvalue_ind<-pnorm(-abs(z_ind))
pvalue_ind<=0.1
print("outdegree top:")
names(gradoS)[which(pvalue_ind<=0.1)]

#Indegree
z_out<-scale(gradoE)
pvalue_out <- pnorm(-abs(z_out))
print("indegree top:")
names(gradoE)[which(pvalue_out<=0.1)]

#Degree
z_degree<-scale(grado)
pvalue_degree <- pnorm(-abs(z_degree))
print("degree top:")
names(grado)[which(pvalue_degree<=0.1)]

#boxplot(grado, gradoS, gradoE, names = c("degree", "outdegree","indegree")) #, plot=FALSE)
####Density## https://kateto.net/netscix2016.html
#The proportion of present edges from all possible edges in the network.
#for directed networks
netDensity <- ecount(net)/(vcount(net)*(vcount(net)-1))
print(c("NET density",netDensity))

###Reciprocity
#The proportion of reciprocated ties (for a directed network).
reciprocity <- reciprocity(net)
print(c("NET RECIPROCITY",reciprocity))

###Transitivity
#Transitivity measures the probability that the adjacent vertices of a vertex are connected. This is sometimes also called the clustering coefficient
#The global transitivity of an undirected graph (directed graphs are considered as undirected ones as well). 
#This is simply the ratio of the triangles and the connected triples in the graph. For directed graph the direction of the edges is ignored.
transitivity <- transitivity(net, type="global")
print(c("Transitivity",transitivity))
#Type of triads for directed networks
triad_census(net)
