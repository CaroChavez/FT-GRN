#Clear environment
rm(list=ls())
gc(reset = TRUE)

library(gplots)
library(RColorBrewer)
library(BoolNetPerturb)
#library(cowplot)

####Plot single mutants of FT-GRN Model Figure 7
#################################################
####Load WT data
wtData <- readRDS(file = "wtdataframe.rds") #class data.frame from ftgrnWT.R

####Datos de mutantes
load("singlemut_ftgrn.RData", verbose = FALSE)

########Rearrange wt data
wt <- wtData$percentage
orden <- c(row.names(mKOpercent)[2:length(row.names(mKOpercent))])
wt <- wt[orden]

wt <- c(NA, wt)

names(wt) <- rownames(mKOpercent)

#########Add wt data to first column
mKOpercent <- cbind(wt, mKOpercent)

mOXpercent <- cbind(wt, mOXpercent)

#######Plot
greys <- colorRampPalette(brewer.pal(9,"Greys"))(100)

greys2 <- colorRampPalette(c(greys[5], greys[50], greys[100]))(100)

ordenLOF <- c(1,2,5,6,9,10,11,13,16,20,21,8,12,7,15,22,17,14,19,4,18,3) #ordenados por tipos de cambios, si son iguales a wt, si pierden y si ganan atractores
lista <- list(1,2)
lossMap <-
  heatmap.2(
    x = mKOpercent[,ordenLOF],
    Rowv = FALSE,
    Colv = FALSE, 
    dendrogram = "none",
    cellnote = round(mKOpercent[,ordenLOF], digits = 2),
    notecol = "black",
    notecex = .7,
    trace = "none",
    key = TRUE,
    density.info = "none",
    key.xlab = "Basin size (%)",
    lhei = c(.7,2),
    lwid=c(2,7),
    breaks = c(0:100),
    col = greys2,
    na.color = "white",
    main = "Loss of function",
    ylab = "Attractor phenotype",
    xlab = "Genotype",
    cexRow = 0.9,
    cexCol = 0.9,
    sepwidth=c(0.01,0.01),
    sepcolor="white",
    colsep=1:ncol(mKOpercent[,ordenLOF]+1),
    rowsep=1:nrow(mKOpercent[,ordenLOF]+1)
  )

ordenGOF <- c(1,2,9,11,13,21,5,6,18,8,22,12,16,17,15,14,7,10,3,4,19,20) #ordenados por tipos de cambios, si son iguales a wt, si pierden y si ganan atractores

gainMap <-  heatmap.2(
    x = mOXpercent[, ordenGOF],
    Rowv = FALSE,
    Colv = FALSE,
    dendrogram = "none",
    lhei = c(.7,2),
    lwid=c(2,7),
    cellnote = round(mOXpercent[,ordenGOF], digits = 2),
    notecol = "black",
    notecex = .7,
    trace = "none",
    key = TRUE,
    density.info = "none",
    key.xlab = "Basin size (%)",
    breaks = c(0:100),
    col = greys2,
    na.color = "white",
    main = "Constitutive Expression",
    ylab = "Attractor phenotype",
    xlab = "Genotype",
    cexRow = 0.9,
    cexCol = 0.9,
    sepwidth=c(0.01,0.01),
    sepcolor="white",
    colsep=1:ncol(mKOpercent[,ordenLOF]+1),
    rowsep=1:nrow(mKOpercent[,ordenLOF]+1)
  )

###Helpfull to examine "other" attractors
#
#others <- dataOX$FT[which(dataOX$FT$etiq == "Other"),]
#others
#sapply(dataKO$AP2$involvedStates[which(dataKO$AP2$etiq == "Other")], FUN= int2binState, node.names = net$genes)
#int2binState(5307041, net$genes)
###For fixed point attractors
#heatmap.2(x = sapply(dataKO$LFY$involvedStates[which(dataKO$LFY$etiq == "Other")], FUN= int2binState, node.names = net$genes),          
#          Rowv = FALSE, Colv = TRUE, dendrogram = "none",
#          trace = "none", key = FALSE,density.info = "none", 
#          col = colorRampPalette(c("red", "#4daf4a"))(2),
#          main = "Mutant phenotypes", cexCol = 0.7)
###For cyclic attractors
#sapply(c(4655224,4670560), FUN= int2binState, node.names = net$genes)
#######Caro Chavez JULIO 2019 UNAM
