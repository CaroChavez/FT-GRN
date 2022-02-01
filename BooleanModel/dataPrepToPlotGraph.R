###Returns a csv table with colapsed edges and information of the type of regulatory interaction between 2 nodes
#(positive regulation or negative regulation and the direction),
#if 2 nodes have feedback, the 2 directed edges are colapsed into 1  but the information of both ends of the edge is maintained.
# the csv file is then used in yED or Cytoscape for plotting the net
#Type of regulation: 0 = negative regulation, 1 positive regulation, NA no regulation in that direction of the edge
library("BoolNet")
library("BoolNetPerturb")
library("igraph")

######Create iGraph object
net<- loadNetwork(file = "ftgrn.txt")

grafo <- plotNetworkWiring(net, plotIt = FALSE)

grafo <- as.undirected(grafo, mode = "collapse")

write.graph(grafo, file="grafo.csv",format="ncol") 

grafoS <- read.csv(file = "grafo.csv", header = FALSE, sep = " ", stringsAsFactors = FALSE)
names(grafoS) <- c("g1", "g2")

###target node edges sign
interaccion <- getNetTopology(net)

###Should TRUE
length(unique(c(interaccion$Source,interaccion$Target))) == length(net$genes)

unique(c(interaccion$Source,interaccion$Target)) %in% net$genes
net$genes %in% c(interaccion$Source,interaccion$Target)

interaccion ###inpect they are all defined before using next function

####Fuction to write the type of interaction as 1 if its positive and 0 if negative
vector<- sapply(interaccion$Interaction,FUN= function(i){
  if(i=="+") {return(1)}
  if(i=="-") {return(0)}
  else {return(0)} #case of FCA to FLC wich is negative
}
)

interaccion <- cbind.data.frame(interaccion,vector)

####Functions to get the type of arrow between a pair of genes

##Type of interaction from the directed edge g1->g2
indiceTF<-function(i){which(grafoS$g1[i]==interaccion$Source & grafoS$g2[i]==interaccion$Target)}

##Type of interaction from the directed edge g2->g1
indiceSF<-function(i){which(grafoS$g2[i]==interaccion$Source & grafoS$g1[i]==interaccion$Target)}

indicesT <- unlist(sapply(c(1:length(grafoS$g1)), FUN = indiceTF))

indicesS <- unlist(sapply(c(1:length(grafoS$g2)), FUN = indiceSF))

typeT <- interaccion$vector[indicesT] #58
##Type of interaction between source node g1 to target node g2 of table interaccion

typeS <- interaccion$vector[indicesS] #70
##Type of interaction between source node g2 to target node g1 of table interaccion

############################################################################################
##typeS and TypeT vectors indexes for adding to table
indices <- function(i,x1,l1,x2,l2){which(x1[i]==l1 & x2[i]==l2)}

indicestt <- unlist(sapply(indicesT, FUN = indices, l1=grafoS$g1, x1=interaccion$Source, l2=grafoS$g2, x2=interaccion$Target))

names(typeT) <- indicestt

indicest2 <- unlist(sapply(indicesS, FUN = indices, l1=grafoS$g2, x1=interaccion$Source, l2=grafoS$g1, x2=interaccion$Target))
names(typeS) <- indicest2

##Column for grafoS table with the type of arrow that gets to node g2
typeG2 <- rep(NA,length=c(nrow(grafoS)))

for(i in c(0:length(indicestt))){
  typeG2[indicestt[i]] <- typeT[i]
}

##Column for grafoS table with the type of arrow that gets to node g1
typeG1 <- rep(NA,length=c(nrow(grafoS)))

for(i in c(0:length(indicest2))){
  typeG1[indicest2[i]]<-typeS[i]
}

(length(typeG2)==nrow(grafoS)) & length(typeG2)==length(typeG2)
###########################################################################################
#grafoS is the table for plotting the network in yED or Cytoscape
grafoS <- cbind(grafoS,typeG2,typeG1)
##
conc<-function(i){paste(as.character(grafoS$typeG1[i]),as.character(grafoS$typeG2[i]))}

g1g2 <- sapply(c(1:nrow(grafoS)), conc)

grafoS <- cbind(grafoS,g1g2)

graphX <- graph_from_data_frame(grafoS, vertices = net$genes)

###Save data in several formats

#write.graph(graphX, file = "graphX.gml", format = "gml")
#plot(graphX)

##Tabla que 
###Guarda csv
#write.csv(grafoS, file = "graphS.csv")

#write.csv(net$genes, file = "nodenames.csv")
####
#Caro Chavez UNAM March 2020