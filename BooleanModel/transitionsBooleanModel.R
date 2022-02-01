library("BoolNet")
library("BoolNetPerturb")
library("igraph")

####Load FT GRN net and calculate attractors
net<- loadNetwork(file = "ftgrn.txt")

attr <- getAttractors(net)
##############
####Labeling rules
labels <- c("JVM", "AVM","IM", "FM")

rules <- c("!SOC1 & !LFY & !AP1 & !SPL3 & !SPL9",
           "(SPL3 | SPL9) & !SOC1 & !LFY & !AP1",
           "TFL1 & SOC1 & !AP1 & !LFY",
           "(AP1 | LFY) & !TFL1")

labels.rules <- data.frame(labels, rules)

##### Dataframe attractor info
attr.df <- attractorToDataframe(attr)
TablaB <- matrix(ncol = 2)

###########################################################################
####Label attractors as meristem types
attrEt <- labelAttractors(attr, labels.rules)

attr.df$label <- sapply(attrEt, function(label) {
  paste(as.character(label), collapse='/')
})

###############################################################################
####Change bits of each gene state in the atractor vector, one bit at a time, return the 23 states one bit apart from the atractor

flipBits <- function(attractorDec) {
  lista1p <- c()
  for(i in c(1:length(net$genes))){
    vector <- int2binState(attractorDec, net$genes)
    vector[i] = !vector[i]
    z <- bin2intState(vector)
    lista1p <- c(lista1p,z)
    lista1p
  } 
  return(lista1p)
}

####Regresa una matriz con columnas del numero de atractor wt inicial y filas con el gen cuyo estado se cambio.
flip_1by1 <- sapply(attr.df$involvedStates, FUN = flipBits)

rownames(flip_1by1) <- net$genes

##Look for the states in flip_1by1 change from one basin of attraction to another one

#########################Get basin of attraction index where the perturbed state is after 1 bitflip
#
saltoDbasin <- function(estado){
  vecr <- int2binState(unlist(estado), node.names = net$genes)
  attractorN <- getStateSummary(vecr, attractorInfo = attr) #state info and index of basin of attraction
  lista_finalBA <- attractorN[[47]] ##basin of attraction  of final state. Use print(getStateSummary(0, attr)) to modify [47] if using a different net 
  return(lista_finalBA)
}

###basinF is a matrix with the basin of attraction index of the new state
basinF <- apply(flip_1by1, FUN = saltoDbasin, MARGIN = c(1,2))

##Phenotype is a matrix with the phenotype of the attractor of the new basin of attraction reached
phenotype <- apply(basinF, FUN= function(i) attr.df$label[i], MARGIN = c(1,2))

################################################################
################################################################
###Jumps of phenotype
###Dataframe with nodes's state flip-change changed the phenotype of the attractor of the new basin of attraction reached

brincos <- sapply(c(1:length(attr.df$involvedStates)), function(i) which(phenotype[,c(i)] != attr.df$label[i]))
names(brincos) <- attr.df$involvedStates

#Rearrange data for graph
##
initialAttractorN <- names(unlist(brincos))

initialAttractor <- sapply(c(1:length(initialAttractorN)), FUN = function(i) as.numeric(gsub("[a-zA-Z. ]", "", initialAttractorN[i])))
##

####NodeP is the node that caused the change in phenotype when perturbed for only one time step
nodeP <- sapply(c(1:length(initialAttractorN)), FUN = function(i) gsub("[0-9. ]", "", initialAttractorN[i]))
##

initialBasin <- match(initialAttractor, attr.df$involvedStates)
##

initialPhenotype <- sapply(initialBasin,FUN= function(i) attr.df$label[i])
##

finalBasin <- unlist(sapply(c(1:length(attr.df$involvedStates)), function(i) which(phenotype[,c(i)] != attr.df$label[i])), use.names= FALSE)

finalBasin <- sapply(c(1:length(initialAttractorN)), FUN= function(i)basinF[finalBasin[i],initialBasin[i]])
##

finalAttractor <- sapply(finalBasin,FUN= function(i) attr.df$involvedStates[i])
##

finalPhenotype <- sapply(finalBasin,FUN= function(i) attr.df$label[i])

#####Dataframe with data for the graph
datosNet <- data.frame(initialBasin,initialAttractor, initialPhenotype,nodeP,finalBasin,finalAttractor,finalPhenotype)

#####Save data for cytoscape
#write.csv(datosNet, file = "transitionsBoolean.csv") ##Open file and delete first column before using in cytoscape

###########################
######Caro Chavez 061119
###