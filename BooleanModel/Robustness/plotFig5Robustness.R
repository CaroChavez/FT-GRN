library("BoolNet")
library("ggplot2")
library("cowplot")

###Graphs for Figure 5, takes data generated in file Fig5_Robustness.R

##Read net
net<- loadNetwork(file = "ftgrn.txt")

#####################
ttSize <- sapply(c(1:23),FUN= function(x) length(net$interactions[[x]]$fun))

indices <- which(ttSize < 3000)

###Load all data file
resultados <- readRDS("fig5Data.rds")

####Data for histogram
indicesL=c(1:length(net$genes))

##Percentage of attractor conservation = cuenta1
cuenta1<- sapply(indicesL, FUN = function(i) (sum(resultados[[i]][,1])*100/nrow(resultados[[i]])))
names(cuenta1) <- names(resultados)

#Percentage of new attractors = cuenta2
cuenta2 <- sapply(indicesL, FUN = function(i)(sum(resultados[[i]][,2])*100/nrow(resultados[[i]])))
names(cuenta2) <- names(resultados)

##
data <- data.frame(
  gene=names(resultados),  
  attrCon=cuenta1,
  newAttr=cuenta2,
  sameAtts= 100-(cuenta1+cuenta2)
)
#plotHistogram
#barplot(t(matrix(cuenta1)), names.arg = net$genes)

AtractorConservation<- ggplot(data, aes(x=gene, y=attrCon)) + 
  geom_bar(stat = "identity") +xlab("Perturbed Function")+ylab("Attractor Conservation (%)")+
  theme_bw(base_size = 12) +
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  scale_y_continuous(limits = c(0,100))
  
#New Atractors
NewAttractors<- ggplot(data, aes(x=gene, y=newAttr)) + 
  geom_bar(stat = "identity") +xlab("Perturbed Function")+ylab("New attractors (%)")+
  theme_bw(base_size = 12) +
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  scale_y_continuous(limits = c(0,100))

###Rearrange data for stacked bars plot:

orden <- c(indices, which(net$genes=="AP1"), which(net$genes=="SOC1"))

nodos <- rep(names(resultados),3)
medida <- c(rep(c("identical"),length(orden)), rep(c("subset"),length(orden)),rep(c("new"),length(orden)))
datos <- c(cuenta1,c(100-(cuenta1+cuenta2)),cuenta2)
regulators <- sapply(c(orden),FUN=function(i) return(length(net$interactions[[i]]$input)))
names(regulators)<- net$genes[orden]
regulators <- rep(regulators,3)
data1 <- data.frame(nodos, medida, datos, regulators)

data1$medida <- factor(data1$medida, levels = c("identical", "subset", "new"))

# Stacked bars plot
barras<-ggplot(data1, aes(fill= medida, y=datos, x=nodos)) + 
  geom_bar(position=position_stack(reverse = TRUE), stat="identity") +xlab("Perturbed Function")+ylab("Attractor Recovery (%)")+
  theme_bw(base_size = 12) +
  theme(axis.text.x=element_text(angle=90,hjust=1), axis.title=element_text(size=12))+
  scale_y_continuous(limits = c(0,100))+
  scale_fill_manual(values = c("grey20", "grey50", "grey85"))+
  labs(fill="") #+

###
x=regulators[1:length(orden)] #indegree
y=cuenta1 #percentage of identical attractor recovery

data2<- data.frame(x,y)
indegree1 <- ggplot(data2, aes(x=x,y=y))+geom_point()+
  xlab("Indegree")+
  ylab("Identical Attractor Recovery (%)")+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,14), breaks = seq(0,14,2))+
  theme(axis.title=element_text(size=12))

indegree1

### Plot the number of regulators (size of true table) vs identical attractor recovery
#indegree<- plot(x=regulators[1:length(net$genes)], y=cuenta1, xlab = c("Indegree"), ylab = c("Identical attractor recovery (%)"))

#plot_grid(barras, indegree1, labels = "AUTO")

#####
#########################################################
########New attractor basin size plot, data from 3rd column of resultados
newBS <- function(idx){
  nBSp <- (resultados[[idx]][,3]*100)/(2**length(net$genes))
  return(round(nBSp, digits = 2))
}

###################
basinP <- sapply(c(1:23), FUN=newBS)
names(basinP) <- net$genes[orden]

#### order indexes
genes<- unlist(sapply(orden, FUN= function(i) rep(net$genes[i], nrow(resultados[[net$genes[i]]]))))

datos3 <- data.frame(genes, unlist(basinP),stringsAsFactors = FALSE)
names(datos3) <- c("genes", "basinP")

#####Graph
violin <- ggplot(data=datos3,aes(x=genes, y=basinP))+
  geom_count(aes(size = ..prop.., group = genes),show.legend = FALSE) +
  scale_size_area(max_size = 4)+
  #geom_count( show.legend = TRUE) +
  #geom_violin(adjust=8)+
  #geom_point(size=0.4, color="grey25")+
  xlab("Perturbed Function")+
  ylab("New Attractors Basin Size (%)")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust=1), axis.title=element_text(size=12))+
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100,10))
#geom_boxplot(width = .3)

###Promedio nuevos attractores basin size
promedio <- ggplot(data=datos3,aes(x=genes, y=basinP, group = genes))+
  stat_summary(aes(group = genes), fun.y = mean, geom = 'point', size=2, alpha=0.9) +
  xlab("Perturbed Function")+
  ylab("Mean New Attractor's Basin Size (%)")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust=1), axis.title=element_text(size=12))+
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100,10)) #+
  #geom_boxplot(width = .3)
promedio

#violin
######## Plot figure 5
#plot_grid(barras, violin, indegree1, labels = "AUTO", ncol = 3, rel_widths = c(1.2, 1, 1))
plot_grid(barras, promedio, indegree1, labels = "AUTO", ncol = 3, rel_widths = c(1.2, 1, 1))
##Caro Chavez Marzo 2020