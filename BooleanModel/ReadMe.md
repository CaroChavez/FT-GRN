# FT-GRN Boolean Model 
Logic rules of the Flowering Transition GRN Boolean Model are in [ftgrn.txt](BooleanModel/ftgrn.txt) file.

Data and Code for figures:
  * Figure 2.- Use [dataPrepToPlotGraph.R](BooleanModel/dataPrepToPlotGraph.R) to create [ftgrnGraph.csv](BooleanModel/ftgrnGraph.csv). This last file was used  to visualize the [FT-GRN graph](BooleanModel/ftgrnGraph.graphml) with [yED](https://www.yworks.com/products/yed).
  * Figure 4 and Supplementary Figure 3.- Use [ftgrnWT.R](ftgrnWT.R).
  * Figure 6.- Data was calculated with [singleMutants.R](BooleanModel/singleMutants.R). To plot figure 6 use [plotSingleMutants.R](BooleanModel/plotSingleMutants.R) and data files [singlemut_ftgrn.RData](BooleanModel/singlemut_ftgrn.RData) and [wtdataframe.rds](BooleanModel/wtdataframe.rds).
  * Figure 7 
  * Supplementary Figure 2 can be plotted with [graphProperties.R](BooleanModel/graphProperties.R).
  * Supplementary Figure 6.- Double mutant analysis in [doubleMutants.R](BooleanModel/doubleMutants.R) gives output file [doubleMutants.rds](BooleanModel/doubleMutants.rds). To get plot of Supplementary Figure 6 use [plotFigS6.R](BooleanModel/plotFigS6.R) and data from [doubleMutants.rds](BooleanModel/doubleMutants.rds), [wtdataframe.rds](BooleanModel/wtdataframe.rds) and [singlemut_ftgrn.RData](BooleanModel/singlemut_ftgrn.RData).
  * Supplementary Figure 7.- Data was calculated with [transitionsBooleanModel.R](BooleanModel/transitionsBooleanModel.R) and gives output file [transitionsBoolean](BooleanModel/transitionsBoolean.csv) that was used to visualize net in Supplementary Figure 7 with [Cytoscape](https://cytoscape.org/).

Data and Code for Figure 5 and Supp Figures 4 and 5 can be found on [Robustness](https://github.com/CaroChavez/FT-GRN/tree/main/BooleanModel/Robustness) folder.
