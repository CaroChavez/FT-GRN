#Clear environment
rm(list=ls())
gc(reset = TRUE)
#Load packages
library(BoolNet)
library(BoolNetPerturb)
library(gplots)
library(RColorBrewer)

#########Gets double mutant attractors and label them as meristem types
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

attrMUT <- function(net, label.rules, genM, sep= "/"){
  #genM <- genIndex[c(tabla[,x]==1)]
  net <- fixGenes(network = net, fixIndices =  genM, values = 0)
  attrsM <- getAttractors(net)
  etiqueta <- labelAttractorsU(attrsM, label.rules)
  etiqueta <- sapply(etiqueta,function(a) paste(as.character(a), collapse = "/"))
  etiq <- sapply(etiqueta, function(x) if(nchar(x) > 3) x = "Other" else x=x)
  attrsM <- attractorToDataframe(attrsM)
  attrsM <- cbind(attrsM, etiq, etiqueta)
}

##################Rules for labeling attractors as types of meristems
######

labels <- c("JVM", "AVM","IM", "FM", "IFM")
rules <- rules <- c("!SOC1 & !LFY & !AP1 & !SPL3 & !SPL9",
                    "(SPL3 | SPL9) & !SOC1 & !LFY & !AP1",
                    "TFL1 & SOC1 & !AP1 & !LFY",
                    "(AP1 | LFY) & !TFL1",
                    "(AP1 | LFY) & TFL1")

labelrules <- data.frame(labels, rules)


#################Get double mutant networks
######

flcsvp <- attrMUT(net= net, label.rules = labelrules, genM = c("SVP", "FLC"))

svpmir156 <- attrMUT(net= net, label.rules = labelrules, genM = c("SVP", "MIR156"))


#####################Function to get percentage of basin size that gets to each phenotype
attractors.by.label <- function(atracM) data.frame(
  states = tapply(atracM$involvedStates,atracM$etiq, paste),
  basin = tapply(atracM$basinSize, atracM$etiq, sum),
  percentage = tapply(atracM$basinSize, atracM$etiq, sum)*100/sum(atracM$basinSize))

##################Size of the basin of attraction that gets to a phenotype

percentageMM <- lapply(list(flcsvp, svpmir156), attractors.by.label)


#################Function to get mutant phenotype info in matrix
matrizj <- matrix(nrow = length(labels), ncol = 2)
colnames(matrizj) <- c("svp flc",
                       "svp miR156"
                      )
row.names(matrizj) <- c("JVM", "AVM", "IM", "FM", "Other")


fillmatrix <- function(index, datos, matriz2F){
  values <- datos[[index]][match(rownames(matriz2F), rownames(datos[[index]])),]
  values <- unlist(values$percentage)
  values
}

#######################Matrix with mutant data

multiM <- sapply(1:length(percentageMM),fillmatrix, datos= percentageMM, matriz2F=matrizj)
colnames(multiM) <- colnames(matrizj) 
rownames(multiM) <- c("JVM", "AVM", "IM", "FM", "Other")

# Uncomment to save miltiM data for plotting
saveRDS(multiM, file = "doubleMutants.rds")


#####Save data
st<-format(Sys.time(), "%Y-%m-%d-%a")
filename <- paste("mutantesM2_",st, ".RData", sep = "")

salvar <- 0

if(salvar == 1) save.image(file = filename)
###

#Caro Chavez UNAM, Mexico Septiembre 2021
