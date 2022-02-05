#Clear environment
rm(list=ls())
gc(reset = TRUE)

library(RColorBrewer)
library(BoolNetPerturb)
library(ggplot2)

####Plots figure 7 with wt and single and double loss of function mutant data in SD dynamics (fix CO=0). Genes mutated XAL2 and SPL9.
####We have prevously calculated multiple mutants with function fixGenes(network, fixIndices, values) from BoolNet package, then named the attractors as meristem types 
## based on same labeling rules for WT attractors, then we summed the basin of attraction percentage for each type of meristem. 
### Example of how to calculate multiple mutants can be found in doubleMutants.R file
#################################################

####Load Basin of attraction data from SD mutant dynamics

multiplesMX <- readRDS(file = "multiMutSD.RDS")

#Data frame
df <- data.frame(Phenotype=rep(colnames(multiplesMX), each=5),
                       Genotype=rep(rownames(multiplesMX), times=8),
                       points= unlist(as.list(multiplesMX)))

###### Color palette
greys <- colorRampPalette(brewer.pal(9,"Greys"))(100)
greys3 <- colorRampPalette(c(greys[5], greys[50], greys[100]))(5)

#With WT rules for GA
df1 <- df[c(1:5,11:15,21:25,31:35),]


sinGA <- ggplot(df1, aes(fill=factor(Genotype, levels = c("FM", "Other", "IM", "AVM", "JVM" )), y=points, x=factor(Phenotype, levels = colnames(multiplesMX)[c(1,3,5,7)]))) + 
  geom_bar(position="stack", stat='identity')+ 
  labs(x=' ', y='Basin size %') +
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold')) +
  scale_fill_manual('Phenotype', values=greys3)

# With exogenous GA (GA=1)
df2 <- df[c(6:10,16:20,26:30,36:40),]

conGA <- ggplot(df2, aes(fill=factor(Genotype, levels = c("FM", "Other", "IM", "AVM", "JVM" )), y=points, x=factor(Phenotype, levels = colnames(multiplesMX)[c(2,4,6,8)]))) + 
  geom_bar(position="stack", stat='identity')+ 
  labs(x=' ', y='Basin size %') +
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold')) +
  scale_fill_manual('', values=greys3)

par(mfrow=c(1,2)) 

sinGA
conGA

#######Caro Chavez Abril 2021 UNAM
