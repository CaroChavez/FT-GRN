#Clear environment
rm(list=ls())
gc(reset = TRUE)
###This code calculates which nodes can induce a change in meristem type when its decay rate is increased in the FT-GRN Continous model
library("ggplot2")
library("BoolNet")
library("BoolNetPerturb")
library("deSolve")
library("gplots")
library("RColorBrewer")
#####Load Boolean Network and calculates its attractors
net<- loadNetwork(file = "ftgrn.txt") #Choose path to file ftgrn.txt
attr <- getAttractors(net)
####
####Labeling rules for the clasification of attractors in meristem types
labels <- c("JVM", "AVM","IM", "FM","IFM")

rules <- c("!SOC1 & !LFY & !AP1 & !SPL3 & !SPL9",
           "(SPL3 | SPL9) & !SOC1 & !LFY & !AP1",
           "TFL1 & SOC1 & !AP1 & !LFY",
           "(AP1 | LFY) & !TFL1",
           "(AP1 | LFY) & TFL1")

labels.rules <- data.frame(labels, rules)

##### Gives table of initial conditions and the nodes that will be modified in each perturbation (Only nodes that are ON on the attractor)
attr.df <- attractorToDataframe(attr)
tabla <- matrix(ncol = 2)

for(fila in c(1:length(attr.df$involvedStates))){
  attBool <- int2binState(attr.df$involvedStates[fila], net$genes)
  TablaA <- matrix(data = NA, ncol=2, nrow = sum(attBool))
  TablaA[,1] <- as.integer(attr.df$involvedStates[fila])
  TablaA[,2] <- as.integer(which(attBool==TRUE))
  tabla <- rbind(tabla, TablaA)
  tabla
}

tabla <- na.omit(tabla)
colnames(tabla) <- c("InitialAttr", "GeneIndexON")

####################Functions##########################
########Function to give the final vector steady state as a decimal number that can be converted to binary####
Comparar <- function(Vector){
  vecbinario <- as.vector(sapply(round(Vector), FUN= as.integer))# rounds the steady state of each gene, which is a vector of decimals close to 0 or 1
  asDec <- bin2intState(vecbinario)
  #At <- labelState(state= asDec, node.names = net$genes, label.rules= labels.rules)  # label binary state
  return(asDec)
}

###############Function to transform a Boolean Model to a ODE Model
net.ode <- booleanToODE(net)


net.ode$func <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # Input Nodes
    w_AGL24 = max(min(SOC1,GA),VER)
    w_AP1 = min(max(1-TFL1,PNY),max(max(max(max(max(LFY,min(min(FT,FD),SPL3)),min(SPL9,1-GA)),min(AGL24,SVP)),AP1),XAL2))
    w_LFY = max(1-TFL1,min(PNY,max(max(max(max(max(max(GA,XAL2),min(SPL3,FD)),AP1),min(AGL24,SOC1)),FUL),SVP)))
    w_AP2 = min(min(min(max(min(SVP,FLC),AP1),1-SOC1),1-FUL),1-MIR172)
    w_AP2L = min(1-MIR172,max(max(min(FLC,SVP),min(min(min(1-SOC1,1-FUL),1-PNY),1-AP1)),AP2))
    w_GA = min(min(min(1-LFY,1-AP1),1-SVP),AGE)
    w_FLC = max(min(min(1-FCA,1-VER),max(1-FT,1-FD)),min(min(FCA,1-AGE),1-VER))
    w_FD = min(min(1-AP1,1-FLC),max(max(max(PNY,LFY),AGE),GA))
    w_FT = min(min(min(max(1-FLC,min(FLC,GA)),1-SVP),1-AP2L),CO)
    w_FUL = min(min(max(max(min(min(FT,FD),SPL3),min(GA,SPL9)),SOC1),1-AP1),max(1-AP2,min(1-SVP,1-AGL24)))
    w_MIR156 = max(1-AGE,min(AP2,1-PNY))
    w_MIR172 = min(max(max(min(1-SVP,1-FLC),1-AP2),FCA),max(max(SPL9,GA),SOC1))
    w_PNY = PNY
    w_SOC1 = min(min(min(min(min(min(1-SVP,1-FLC),1-AP2),1-AP2L),1-AP1),GA),max(max(max(max(max(CO,min(FT,FD)),SPL9),min(SOC1,AGL24)),XAL2),FUL))
    w_SPL3 = min(1-MIR156,max(max(SOC1,min(FT,FD)),GA))
    w_SPL9 = min(min(max(1-SVP,1-FLC),1-MIR156),max(1-AP1,GA))
    w_SVP = min(1-AP1,1-FCA)
    w_TFL1 = max(max(XAL2,min(LFY,1-AP1)),1-PNY)
    w_XAL2 = min(min(max(max(max(CO,GA),SPL9),AGE),1-AP1),1-SOC1)
    w_CO = CO
    w_FCA = FCA
    w_VER = VER
    w_AGE = AGE
    
    
    # Rates of Change
    dAGL24 = ((-exp(0.5*h)+exp(-h*(w_AGL24-0.5)))/((1-exp(0.5*h))*(1+exp(-h*(w_AGL24-0.5)))))-(alphaAGL24*AGL24)
    dAP1 = ((-exp(0.5*h)+exp(-h*(w_AP1-0.5)))/((1-exp(0.5*h))*(1+exp(-h*(w_AP1-0.5)))))-(alphaAP1*AP1)
    dLFY = ((-exp(0.5*h)+exp(-h*(w_LFY-0.5)))/((1-exp(0.5*h))*(1+exp(-h*(w_LFY-0.5)))))-(alphaLFY*LFY)
    dAP2 = ((-exp(0.5*h)+exp(-h*(w_AP2-0.5)))/((1-exp(0.5*h))*(1+exp(-h*(w_AP2-0.5)))))-(alphaAP2*AP2)
    dAP2L = ((-exp(0.5*h)+exp(-h*(w_AP2L-0.5)))/((1-exp(0.5*h))*(1+exp(-h*(w_AP2L-0.5)))))-(alphaAP2L*AP2L)
    dGA = ((-exp(0.5*h)+exp(-h*(w_GA-0.5)))/((1-exp(0.5*h))*(1+exp(-h*(w_GA-0.5)))))-(alphaGA*GA)
    dFLC = ((-exp(0.5*h)+exp(-h*(w_FLC-0.5)))/((1-exp(0.5*h))*(1+exp(-h*(w_FLC-0.5)))))-(alphaFLC*FLC)
    dFD = ((-exp(0.5*h)+exp(-h*(w_FD-0.5)))/((1-exp(0.5*h))*(1+exp(-h*(w_FD-0.5)))))-(alphaFD*FD)
    dFT = ((-exp(0.5*h)+exp(-h*(w_FT-0.5)))/((1-exp(0.5*h))*(1+exp(-h*(w_FT-0.5)))))-(alphaFT*FT)
    dFUL = ((-exp(0.5*h)+exp(-h*(w_FUL-0.5)))/((1-exp(0.5*h))*(1+exp(-h*(w_FUL-0.5)))))-(alphaFUL*FUL)
    dMIR156 = ((-exp(0.5*h)+exp(-h*(w_MIR156-0.5)))/((1-exp(0.5*h))*(1+exp(-h*(w_MIR156-0.5)))))-(alphaMIR156*MIR156)
    dMIR172 = ((-exp(0.5*h)+exp(-h*(w_MIR172-0.5)))/((1-exp(0.5*h))*(1+exp(-h*(w_MIR172-0.5)))))-(alphaMIR172*MIR172)
    dPNY = ((-exp(0.5*h)+exp(-h*(w_PNY-0.5)))/((1-exp(0.5*h))*(1+exp(-h*(w_PNY-0.5)))))-(alphaPNY*PNY)
    dSOC1 = ((-exp(0.5*h)+exp(-h*(w_SOC1-0.5)))/((1-exp(0.5*h))*(1+exp(-h*(w_SOC1-0.5)))))-(alphaSOC1*SOC1)
    dSPL3 = ((-exp(0.5*h)+exp(-h*(w_SPL3-0.5)))/((1-exp(0.5*h))*(1+exp(-h*(w_SPL3-0.5)))))-(alphaSPL3*SPL3)
    dSPL9 = ((-exp(0.5*h)+exp(-h*(w_SPL9-0.5)))/((1-exp(0.5*h))*(1+exp(-h*(w_SPL9-0.5)))))-(alphaSPL9*SPL9)
    dSVP = ((-exp(0.5*h)+exp(-h*(w_SVP-0.5)))/((1-exp(0.5*h))*(1+exp(-h*(w_SVP-0.5)))))-(alphaSVP*SVP)
    dTFL1 = ((-exp(0.5*h)+exp(-h*(w_TFL1-0.5)))/((1-exp(0.5*h))*(1+exp(-h*(w_TFL1-0.5)))))-(alphaTFL1*TFL1)
    dXAL2 = ((-exp(0.5*h)+exp(-h*(w_XAL2-0.5)))/((1-exp(0.5*h))*(1+exp(-h*(w_XAL2-0.5)))))-(alphaXAL2*XAL2)
    dCO = ((-exp(0.5*h)+exp(-h*(w_CO-0.5)))/((1-exp(0.5*h))*(1+exp(-h*(w_CO-0.5)))))-(alphaCO*CO)
    dFCA = ((-exp(0.5*h)+exp(-h*(w_FCA-0.5)))/((1-exp(0.5*h))*(1+exp(-h*(w_FCA-0.5)))))-(alphaFCA*FCA)
    dVER = ((-exp(0.5*h)+exp(-h*(w_VER-0.5)))/((1-exp(0.5*h))*(1+exp(-h*(w_VER-0.5)))))-(alphaVER*VER)
    dAGE = ((-exp(0.5*h)+exp(-h*(w_AGE-0.5)))/((1-exp(0.5*h))*(1+exp(-h*(w_AGE-0.5)))))-(alphaAGE*AGE)
    
    
    list(c(dAGL24,dAP1,dLFY,dAP2,dAP2L,dGA,dFLC,dFD,dFT,dFUL,dMIR156,dMIR172,dPNY,dSOC1,dSPL3,dSPL9,dSVP,dTFL1,dXAL2,dCO,dFCA,dVER,dAGE))
  })
}


######Function to calculate final steady state of the ODE model after increasing one gene decay rate in each simulation, 
###takes as parameter a table with the initial state as first column and the gene index to be perturbed as second column
####Ver como cambian en el tiempo los atractores toma el vector de InitialAttr, GeneIndexON, FinalAttr, LabelInitialAttr, LabelFinalAttr, NombresGenesTran
transitionsDecay <- function(vector){
  iniAttr <- vector[1] ##Initial state
  geneK <- as.numeric(vector[2]) ###Gene index to be perturbed
  ###ODE numerical solution
  net.ode <- net.ode #ODE model
  net.ode$parameters[geneK+2]= 5.0 ##PARAMETERVECTOR, FIRST 2 ARE  h=10, w=0.5 and THEN decay rate = 1 FOR ALL GENES.
  odeOut21 <- ode(func = net.ode$func, parms = net.ode$parameters, y = int2binState(iniAttr, net$genes), times = seq(0, 10, 0.1))
  ##Final steady state
  out <- odeOut21[101,2:(ncol(odeOut21))] #Final state reached
  outD<- Comparar(round(out)) #Final state as decimal that can be transformed to Boolean vector
  outD
  }

################################Run simulation
FinalAttr <- apply(tabla,FUN = transitionsDecay,MARGIN = 1) #Takes "tabla" and for each row (pair of attractor and gene to be perturbed) calculates the final steady state

FinalLabel <- sapply(FinalAttr, FUN = labelState, node.names = net$genes, label.rules= labels.rules) #Labels final states

LabelInitialAttr <- sapply(tabla[,1], FUN = labelState, node.names = net$genes, label.rules= labels.rules) #Labels initial states

todojunto <- cbind(tabla,LabelInitialAttr, FinalAttr, FinalLabel) #All info 

colnames(todojunto)[2] <- c("GeneIndexOff")

###Inspects wich genes can induce a change in the meristem type from the initial to the final state after been perturbed
TablaJ <- todojunto[which(FinalLabel!=LabelInitialAttr),]
#### Smaller table for the final figure in yED
TablaTr <- TablaJ[,c(3,5,2)]
TablaTr <- unique(TablaTr)

NombreGenesTran <- sapply(as.integer(TablaTr[,3]), function(gindx) net$genes[gindx]) #Name of the gene that caused a transition
TablaTr <- cbind(TablaTr, NombreGenesTran)

TablaTr <- as.data.frame(TablaTr) 


####To plot the simulations that caused phenotype transitions
plotF <- function(vector){
  iniAttr <- vector[1] ##Atractor inicial
  geneK <- as.numeric(vector[2]) ###Gene a perturbar su tasa de decaimiento
  ###Modelo ODEs y reslucion numerica
  net.ode <- net.ode #booleanToODE(net)#, logic = "Probabilistic")
  net.ode$parameters[geneK+2]= 5.0 ##PARAMETERVECTOR, FIRST 2 ARE  h=10, w=0.5 and THEN decay rate = 1 FOR ALL GENES.
  odeOut21 <- ode(func = net.ode$func, parms = net.ode$parameters, y = int2binState(iniAttr, net$genes), times = seq(0, 10, 0.1))
  ##Name odeout final state
  out <- round(odeOut21[101,2:(ncol(odeOut21))])
  fstate<- Comparar(out)
  flabel <- labelState(fstate, node.names = net$genes, label.rules= labels.rules)
  ilabel <- labelState(iniAttr, node.names = net$genes, label.rules= labels.rules)
  ###Plot
  odeOut21 <- melt(odeOut21[,2:24], id.vars=odeOut21[,1])
  names(odeOut21) <- c("Time", "Gene", "value")
  titulo <- paste("Initial ",ilabel," state: ",vector[1]," .","Final ",flabel, " state: ",fstate,". Gene perturbed: ", net$genes[geneK], sep = "")
  ggplot(data = odeOut21, aes(x=odeOut21$Time ,y = odeOut21$value, factor=Gene)) + geom_line()+xlab("Time") +ylab("Gene expression") +
    ggtitle(titulo)+
    facet_grid(Gene ~.)+
    scale_y_continuous(breaks = c(0,1)) +
    theme_gray()+theme(strip.text.y = element_text(angle = 0), legend.position = "none")
}

#Uncomment below to get the simulations graph for each pair of attractor and gene that can change phenotype when perturbed
#apply(TablaJ,FUN = plotF,MARGIN = 1) 

########Uncomment below to save data for the yED graph shown in figure 8 
#write.csv(TablaTr, file = paste("TablaTransODEunique",Sys.Date(), ".csv", sep = ""))
###Caro Chavez 2021