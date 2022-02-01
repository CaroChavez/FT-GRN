library("BoolNet")
library("BoolNetPerturb")
library("ggplot2")

###Robustness of each attractor to random perturbations in the logic rules of the FT-GRN###
#print(Sys.time())
net<- loadNetwork(file = "ftgrn.txt")

attr<- getAttractors(net)

wtDF <- attractorToDataframe(attr)

set.seed(63)


####Labeling rules
labels <- c("JVM", "AVM","IM", "FM")

rules <- c("!SOC1 & !LFY & !AP1 & !SPL3 & !SPL9",
           "(SPL3 | SPL9) & !SOC1 & !LFY & !AP1",
           "TFL1 & SOC1 & !AP1 & !LFY",
           "(AP1 | LFY) & !TFL1")

labels.rules <- data.frame(labels, rules)


#### Label attractors as phenotypes
attrEt <- labelAttractors(attr, labels.rules)

wtDF$label <- sapply(attrEt, function(label) {
  paste(as.character(label), collapse='/')
})

attractors.by.label <- data.frame(
  states = tapply(wtDF$involvedStates,wtDF$label, paste) #,
  #basin = tapply(wtDF$basinSize, wtDF$label, sum),
  #percentage = tapply(wtDF$basinSize, wtDF$label, sum)*100/sum(wtDF$basinSize)
)

attractors.by.label

#Look for unlabel attractors

unique(wtDF$label)

#Test attractor robustness by fliping randomly one line in the true tables of the FT-GRN
######Function to create a copy of the network with a random perturbation to 1 line of 1 Boolean table.

attrRNR <- function(network=net, wtDFrame=wtDF){
	netP <- perturbNetwork(network = net, perturb = "functions", method = "bitflip", simplify = FALSE, maxNumBits = 1)
  attr2<- getAttractors(netP, type = "synchronous", method="exhaustive", returnTable = FALSE)
  attr2DF <- attractorToDataframe(attr2)
  numAtt <- c(1:length(wtDF$involvedStates))
  resu <- c()
  resu <- sapply(numAtt, FUN=function(i){wtDF$involvedStates[i] %in% attr2DF$involvedStates})
return(resu)
}

#Uncoment next line to calculate results
#resultados <- sapply(c(1:3000),FUN=attrRNR) 

#saveRDS(resultados, file = "suppFig4Data.rds")

###Uncomment next line to read results if previously calculated, if not then uncomment 2 lines before

resultados <- readRDS(file = "suppFig4Data.rds")

####Count how many times each wt Attractor was recovered
counts <- apply(resultados,MARGIN = 1, FUN=sum)
percentageRNR <- counts*100/3000 
percentageRNR

resultados <- cbind(wtDF, percentageRNR)

resultados

####PLOT Supplementary Figure 5
barras <- ggplot(resultados, aes(x=as.character(involvedStates,label), y=round(percentageRNR,digits = 1), fill=label)) + 
  geom_bar(stat = "identity") +xlab("Attractor ID")+ylab("Recovery (%)")+
  theme_bw(base_size = 12) +
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,10))+
  scale_fill_manual( values = c("gray60","gray40", "gray90", "gray15"))+
  labs(fill="Phenotype")

barras
##Caro Chavez Marzo 2020
