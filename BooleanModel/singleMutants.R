#Clear environment
rm(list=ls())
gc(reset = TRUE)
#Load packages
library(BoolNet)
library(BoolNetPerturb)
library(gplots)
library(RColorBrewer)

#########Gets single mutants attractors and label them as phenotypes
###########Load network
net <- loadNetwork("ftgrn.txt")

#####Label attractors unique labels if cyclic attractor
labelAttractorsU <- function (attr, label.rules, node.names = NULL, sep = "/") 
{
  if (!is(attr, "AttractorInfo")) {
    stop("Error: non-valid attractor")
  }
  node.names <- attr$stateInfo$genes
  res <- list()
  for (i in 1:length(attr$attractors)) {
    label <- sapply(attr$attractors[[i]]$involvedStates, 
                    function(state) {
                      state <- int2binState(state, node.names)
                      l <- labelState(state, node.names, label.rules, 
                                      sep = "")
                    })
    if (!is.null(sep)) {
      label <- paste(unique(label), collapse = sep)
    }
    res <- append(res, list(label))
  }
  unlist(res)
}
##################Function to get involved states and basin size of the single mutant networks attractors

attrMUT <- function(net, typeM, genIndex, label.rules, sep= "/"){
  genM <- net$genes[genIndex]
  net <- fixGenes(network = net, fixIndices =  genM, values = typeM)
  attrsM <- getAttractors(net)
  etiqueta <- labelAttractorsU(attrsM, label.rules)
  etiqueta <- sapply(etiqueta,function(a) paste(as.character(a), collapse = "/"))
  etiq <- sapply(etiqueta, function(x) if(nchar(x) > 3) x = "Other" else x=x)
  attrsM <- attractorToDataframe(attrsM)
  attrsM <- cbind(attrsM, etiq, etiqueta)
}


##################Rules for labeling attractors as types of meristems

labels <- c("JVM", "AVM","IM", "FM", "IFM")
rules <- rules <- c("!SOC1 & !LFY & !AP1 & !SPL3 & !SPL9",
                    "(SPL3 | SPL9) & !SOC1 & !LFY & !AP1",
                    "TFL1 & SOC1 & !AP1 & !LFY",
                    "(AP1 | LFY) & !TFL1",
                    "(AP1 | LFY) & TFL1")

labelrules <- data.frame(labels, rules)
#################Get KO and OX data for all genes
dataKO <- lapply(1:length(net$genes), FUN = attrMUT, net = net, typeM = 0, label.rules=labelrules)
names(dataKO) <- t(data.frame(net$genes))

dataOX <- lapply(1:length(net$genes), FUN = attrMUT, net = net, typeM = 1, label.rules=labelrules)
names(dataOX) <- t(data.frame(net$genes))

#####################Function to get percentage of basin size that gets to each phenotype
attractors.by.label <- function(atracM) data.frame(
  states = tapply(atracM$involvedStates,atracM$etiq, paste),
  basin = tapply(atracM$basinSize, atracM$etiq, sum),
  percentage = tapply(atracM$basinSize, atracM$etiq, sum)*100/sum(atracM$basinSize))

##################Size of the basin of attraction that gets to a phenotype

percentageKO <- lapply(dataKO[],attractors.by.label)

percentageOX <- lapply(dataOX[],attractors.by.label)

#################Prepare matrix for data storage
matrizj <- matrix(nrow = length(labels)+1, ncol = length(net$genes)+1)
colnames(matrizj) <- c("WT",net$genes)
row.names(matrizj) <- c("Other", labels)

#################Function for all single mutant info
fillmatrix <- function(index, datos, matriz2F){
  values <- datos[[index]][match(rownames(matriz2F), rownames(datos[[index]])),]
  values <- unlist(values$percentage)
  values
}

#######################Calculate single mutant data

mKOpercent <- sapply(1:length(net$genes),fillmatrix, datos= percentageKO, matriz2F=matrizj)
colnames(mKOpercent) <- net$genes
rownames(mKOpercent) <- c("Other",labels)

mOXpercent <- sapply(1:length(net$genes),fillmatrix, datos= percentageOX, matriz2F=matrizj)
colnames(mOXpercent) <- net$genes
rownames(mOXpercent) <- c("Other",labels)


#####Save data for ploting with plotSingleMutants.R
#st=format(Sys.time(), "%Y-%m-%d-%a")
filename <- paste("singlemut_ftgrn.RData", sep = "")

salvar <- 0

if(salvar == 1) save.image(file = filename)
###
#Caro Chavez UNAM, Mexico Marzo 2020
