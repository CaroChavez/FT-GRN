library("BoolNet")
library("BoolNetPerturb")

###Robustness analysis of the FT-GRN, plots Supplementary Figure 5
#print(Sys.time())
net1<- loadNetwork(file = "ftgrn.txt")

wtData<- getAttractors(net1)
  
wtDF <- attractorToDataframe(wtData)
set.seed(33)

#With measure= "hamming", it returns ($value) the normalized hamming distance (fraction of different bits) between each state and the perturbed copy. 
#Hamming distance should be low for robust networks, take for example the cellcycle network is .107
#mSAMPTrajectoriesHamming <- perturbTrajectories(net1, measure = "hamming", numSamples = 10000,flipBits = 1)
#mSAMPTrajectoriesHamming$value

####### Now we compared to a set of randomly perturbed networks.

#This returns a histogram that compares the network (red line) with the distribution of 1000 topologically similar random networks

plotHamm = 0

if(plotHamm ==1){
  histHam<- testNetworkProperties(net1, testFunction="testTransitionRobustness",
                                  testFunctionParams=list(numSamples=1000),
                                  alternative="less", xlim = c(0,0.15),ylim=c(0,50))
  histHam
  #saveRDS(histHam, file = "suppFig5Data.rds")
}

###
