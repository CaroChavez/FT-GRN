# Clear environment
rm(list=ls())
gc(reset = TRUE)

library(gplots)
library(RColorBrewer)
library(BoolNetPerturb)

#### Plot single and double mutants for Supplementary Figure 6
#################################################
#### Load WT attractors dataframe
wtData <- readRDS(file = "wtdataframe.rds") #class data.frame #wtdataframe.rds

#### Single mutant attractors data
mutantesS <- load("singlemut_ftgrn.RData",verbose = FALSE) #("~/Dropbox/red130320/V4/mutantesSAMv4_2020-03-13-Fri.RData", verbose = FALSE) #

### Double mutants attractors data
multiM <- readRDS("doubleMutants.rds")

multiM <- rbind(multiM[5,], multiM[1:4,],c(NA,NA))
rownames(multiM) <- rownames(mKOpercent)

######## Rearrange data
# WT attractor basin size
wt <- wtData$percentage
names(wt) <- rownames(wtData)

# Single mutants basin size
orden <- c(row.names(mKOpercent)[2:length(row.names(mKOpercent))])

wt <- c(NA, wt[rownames(mKOpercent)[2:5]], NA)

names(wt) <- rownames(mKOpercent)

######### Add wt data to first column, then single mutants and double mutants
multiplesM <- cbind( wt,mKOpercent[,c(11,7,17)], multiM)

etiqM<- c("wt", "mir156", "flc", "svp", "svp flc", "svp mir156")

colnames(multiplesM) <- c(etiqM)

####### Plot color palette
greys <- colorRampPalette(brewer.pal(9,"Greys"))(100)

greys3 <- colorRampPalette(c(greys[5], greys[50], greys[100]))(5)


multiplesMX <- multiplesM[1:5,c(1,2,3,4,6,5)] # deletes last row with NA and rearrange order of columns
multiplesMX[is.na(multiplesMX)] <- 0 # si Colv = TRUE


# Data frame
df <- data.frame(Phenotype=rep(colnames(multiplesMX), each=5),
                       Genotype=rep(rownames(multiplesMX), times=6),
                       points= unlist(as.list(multiplesMX)))
# Plot

library(ggplot2)
doubleM <- ggplot(df, aes(fill=factor(Genotype, levels = c("FM", "Other", "IM", "AVM", "JVM" )), y=points, x=factor(Phenotype, levels = colnames(multiplesMX)))) + 
  geom_bar(position="stack", stat='identity')+ 
  labs(x='', y='Basin size %') +
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold')) +
  scale_fill_manual('Phenotype', values=greys3)


doubleM

###Helpfull to examine "other" attractors
#
#others <- dataOX$FT[which(dataOX$FT$etiq == "Other"),]
#others
#sapply(dataKO$GA$involvedStates[which(dataKO$GA$etiq == "Other")], FUN= int2binState, node.names = net$genes)
#int2binState(5307041, net$genes)
###For cyclic attractors
#sapply(c(4655224,4670560), FUN= int2binState, node.names = net$genes)
#######Caro Chavez Abril 2021 UNAM