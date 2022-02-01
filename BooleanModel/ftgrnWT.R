#Clear environment
rm(list=ls())
gc(reset = TRUE)
#######Boolean Model attractors and comparison with gene expression data Figure 4 and Supplementary Figure 3
library(BoolNet)
library(BoolNetPerturb)
library(gplots)

#### Load logic rules
net <- loadNetwork("ftgrn.txt")

####Get attractors
attr <- getAttractors(net)

####Labeling rules

labels <- c("JVM", "AVM","IM", "FM")

rules <- c("!SOC1 & !LFY & !AP1 & !SPL3 & !SPL9",
           "(SPL3 | SPL9) & !SOC1 & !LFY & !AP1",
           "TFL1 & SOC1 & !AP1 & !LFY",
           "(AP1 | LFY) & !TFL1")

labels.rules <- data.frame(labels, rules)

##### Dataframe with all attractor info
attr.df <- attractorToDataframe(attr)

#### Label attractors as phenotypes
attrEt <- labelAttractors(attr, labels.rules)

attr.df$label <- sapply(attrEt, function(label) {
  paste(as.character(label), collapse='/')
})

attractors.by.label <- data.frame(
  states = tapply(attr.df$involvedStates,attr.df$label, paste),
  basin = tapply(attr.df$basinSize, attr.df$label, sum),
  percentage = tapply(attr.df$basinSize, attr.df$label, sum)*100/sum(attr.df$basinSize)
)

attractors.by.label

#Look for unlabel attractors

unique(attr.df$label)

############################################################################
###Plot attractors by label

par(mfrow=c(1,1))
barplot(attractors.by.label$percentage, 
        xlab = "Attractor phenotype",
        ylab = "Basin of attraction size (%)", 
        ylim = c(0,round(max(attractors.by.label$percentage)+4)), 
        col = c("black"),
        las=2,
        cex.names = 0.6
)

par(mfrow= c(1,4))

plotAttractors(attr,  drawLegend = FALSE, subset = which(attr.df$label == "JVM"), title = "JVM",offColor = "white", onColor = "gray9", mode="table", grouping = list(class = c("g1", "g2"), index=list(c(2:3),list(5:23))))

plotAttractors(attr,  drawLegend = FALSE, subset = which(attr.df$label == "AVM"), title = "AVM",offColor = "white", onColor = "gray9")

plotAttractors(attr,  drawLegend = FALSE, subset = which(attr.df$label == "IM"), title = "IM",offColor = "white", onColor = "gray9")

plotAttractors(attr,  drawLegend = FALSE, subset = which(attr.df$label == "FM"), title = "FM",offColor = "white", onColor = "gray9")

par(mfrow= c(1,1))

attractors.by.label

unique(attr.df$label)

###########################################################################################
####Plot as phenotypes and compare with observed data
lista <- list(attr.df$involvedStates)
binarioAT <- sapply(attr.df$involvedStates,int2binState, net$genes)
colnames(binarioAT) <- attr.df$label

JVM <- sapply(attr.df$label, function(x) x=="JVM")
AVM <- sapply(attr.df$label, function(x) x=="AVM")
IM <- sapply(attr.df$label, function(x) x=="IM")
FM <- sapply(attr.df$label, function(x) x=="FM")

binario<-function(atractores)
{ if(class(atractores)=="matrix") 
  rowMeans(atractores) 
  else atractores}


atr1 <-binario(binarioAT[,JVM])
atr2 <- binario(binarioAT[,AVM])
atr3 <- binario(binarioAT[,IM])
atr4 <- binario(binarioAT[,FM])

matrizAt<- cbind(atr1,atr2,atr3,atr4)
colnames(matrizAt) <- c("JVM",
                        "AVM", 
                        "IM", 
                        "FM")

my_palette <- colorRampPalette(c("white", "grey40", "grey5"))(n = 3)

consenso <- heatmap.2((matrizAt), Rowv = TRUE, Colv = FALSE, #[c(1:19),]
                      col= my_palette, colsep = c(0:4), rowsep = c(0:nrow(matrizAt)), 
                      sepcolor = "black", sepwidth=c(0.01,0.02), scale="none", margins=c(5,10), 
                      dendrogram = "none", trace = "none",cexCol = 0.9,
                      key = FALSE, main = "SAM Model\nAttractors", breaks = c(0,c(0.01,0.99),1))
legend("topleft", xjust= 1, legend=c("OFF","ON","ON/OFF"),
       fill=c("white","grey5","grey40"), border=TRUE, bty="n", cex=0.8)

######################################
##Individual attractors by label all in the same plot
matrizAT1<- cbind(binarioAT[,JVM], binarioAT[,AVM], binarioAT[,IM], binarioAT[,FM])

my_palette2 <- colorRampPalette(c("white", "grey5"))(n = 2)

consenso <- heatmap.2((matrizAT1), Rowv = TRUE, Colv = FALSE, #[c(1:19),]
                      col= my_palette2, colsep = c(0:33), rowsep = c(0:nrow(matrizAT1)), 
                      sepcolor = "black", sepwidth=c(0.01,0.02), scale="none", margins=c(5,10), 
                      dendrogram = "none", trace = "none",cexCol = 0.9,
                      key = FALSE, breaks = c(0,0.1:0.99,1)) #main = "FT-GRN Model\nAttractors"
legend("topleft", xjust= 1, legend=c("OFF","ON"),
       fill=c("white","grey5"), border=TRUE, bty="n", cex=0.8)

######################################
###Compared with observed data
observado <- read.csv("expressiondata.csv", header = TRUE)
row.names(observado) <- observado$Gen
observado <- observado[net$genes,c(4,5,2,3)] #ordered as in attractor data

######Fuction for comparison plot
plotHM <- function(mdatos) {
  heatmap.2(mdatos, Rowv = FALSE, Colv = FALSE, col= my_palette , scale="none", #margins=c(6,10),
            colsep = c(0:3), rowsep = c(0:nrow(mdatos)), 
            sepcolor = "black", sepwidth=c(0.01,0.02),
            dendrogram = "none", trace = "none", key = FALSE, na.color = "grey80", cexCol = 1, breaks = c(0,c(0.01,0.99),1))
  legend("topleft", xjust= 1, legend=c("OFF","ON","ON/OFF", "No Info"),
         fill=c("white","grey5","grey40", "grey80"), border=TRUE, bty="n", cex=0.7)
}

########

########
oJVMm <- cbind(observado$JVM, binario(matrizAt[,1]))
colnames(oJVMm) <- c("JVM\n Observed", "JVM\n Model")
plotHM(oJVMm)
########
########
oAVMm <- cbind(observado$AVM, binario(matrizAt[,2]))
colnames(oAVMm) <- c("AVM\n Observed", "AVM\n Model")
plotHM(oAVMm)
########
########
oIMm <- cbind(observado$IM, binario(matrizAt[,3]))
colnames(oIMm) <- c("IM\n Observed", "IM\n Model")
plotHM(oIMm)
########
########
oFMm <- cbind(observado$FM, binario(matrizAt[,4]))
colnames(oFMm) <- c("FM\n Observed", "FM\n Model")
plotHM(oFMm)
########
###Save WT data for single mutant plot
#saveRDS(attractors.by.label, file = "wtdataframe.rds")
######### Caro Chavez Marzo 2020 UNAM
