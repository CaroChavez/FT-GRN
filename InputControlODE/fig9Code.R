######This code changes the decay rate of controls conected to input nodes for differentiation simulation in the ODE FT-GRN model
###Data and plots for Figure 9
#Clear environment
rm(list=ls())
gc(reset = TRUE)
###LOAD PACKAGES
library("ggplot2")
library("BoolNet")
library("BoolNetPerturb")
library("deSolve")
library("cowplot")
########################

#####Boolean FT-GRN with aditional control nodes (CAGE, CVER, CFCA, CPNY) to convert to ODE
net<- loadNetwork(file = "ftgrnControl.txt")
##Boolean FT-GRN
bnet <- loadNetwork(file = "ftgrn.txt")

####################################
####Data transitions Boolean Model
####First Jump from JVM to AVM
##Initial conditions, JVM states
#        #ld#             #sd#
ci <- c(5178456,3081233,4654168,2294801)
nodeID <- c(which(net$genes=="CVER"), which(net$genes=="CAGE"),which(net$genes=="CVER"), which(net$genes=="CAGE"))
ci_vec1 <- sapply(ci, int2binState, bnet$genes)

###Ading control state to initial state
controles <- matrix(data = c(1,0,0,1), nrow = 4, ncol=4)
rownames(controles) <- net$genes[c(length(net$genes)-3):length(net$genes)]

####Initial conditions for ODE, state (from Boolean model) and control node (from control initial state)
ci_vec <- rbind(ci_vec1, controles)

initc <- apply(ci_vec, FUN = bin2intState, MARGIN = 2)

tabla <- cbind(initc,nodeID)

##REMOVE
rm(ci_vec1)
rm(bnet)

####Labeling rules
labels <- c("JVM", "AVM","IM", "FM","IFM")

rules <- c("!SOC1 & !LFY & !AP1 & !SPL3 & !SPL9",
           "(SPL3 | SPL9) & !SOC1 & !LFY & !AP1",
           "TFL1 & SOC1 & !AP1 & !LFY",
           "(AP1 | LFY) & !TFL1",
           "(AP1 | LFY) & TFL1")

labels.rules <- data.frame(labels, rules)

########Function to convert into decimal state the final state of the ODE simulation
Comparar <- function(Vector){
  vecbinario <- as.vector(sapply(round(Vector), FUN= as.integer))# rounds values of final state of the simulation
  asDec <- bin2intState(vecbinario)
  #At <- labelState(state= asDec, node.names = net$genes, label.rules= labels.rules)  # name final steady state
  return(asDec)
}
############Function to get photoperiod as Long Day or Short Day if CO = 1 or CO =0
photoperiod <-function(estado){if(int2binState(estado, net$genes)[which(net$genes=="CO")]==1)return("LD")
  else return("SD")}

#########ODE MODEL
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
    w_PNY = 1-CPNY
    w_SOC1 = min(min(min(min(min(min(1-SVP,1-FLC),1-AP2),1-AP2L),1-AP1),GA),max(max(max(max(max(CO,min(FT,FD)),SPL9),min(SOC1,AGL24)),XAL2),FUL))
    w_SPL3 = min(1-MIR156,max(max(SOC1,min(FT,FD)),GA))
    w_SPL9 = min(min(max(1-SVP,1-FLC),1-MIR156),max(1-AP1,GA))
    w_SVP = min(1-AP1,1-FCA)
    w_TFL1 = max(max(XAL2,min(LFY,1-AP1)),1-PNY)
    w_XAL2 = min(min(max(max(max(CO,GA),SPL9),AGE),1-AP1),1-SOC1)
    w_CO = CO
    w_FCA = 1-CFCA
    w_VER = 1-CVER
    w_AGE = 1-CAGE
    w_CFCA = CFCA
    w_CVER = CVER
    w_CAGE = CAGE
    w_CPNY = CPNY
    
    
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
    dCFCA = ((-exp(0.5*h)+exp(-h*(w_CFCA-0.5)))/((1-exp(0.5*h))*(1+exp(-h*(w_CFCA-0.5)))))-(alphaCFCA*CFCA)
    dCVER = ((-exp(0.5*h)+exp(-h*(w_CVER-0.5)))/((1-exp(0.5*h))*(1+exp(-h*(w_CVER-0.5)))))-(alphaCVER*CVER)
    dCAGE = ((-exp(0.5*h)+exp(-h*(w_CAGE-0.5)))/((1-exp(0.5*h))*(1+exp(-h*(w_CAGE-0.5)))))-(alphaCAGE*CAGE)
    dCPNY = ((-exp(0.5*h)+exp(-h*(w_CPNY-0.5)))/((1-exp(0.5*h))*(1+exp(-h*(w_CPNY-0.5)))))-(alphaCPNY*CPNY)
    
    
    list(c(dAGL24,dAP1,dLFY,dAP2,dAP2L,dGA,dFLC,dFD,dFT,dFUL,dMIR156,dMIR172,dPNY,dSOC1,dSPL3,dSPL9,dSVP,dTFL1,dXAL2,dCO,dFCA,dVER,dAGE,dCFCA,dCVER,dCAGE,dCPNY))
  })
}


#########################################################
####Second jump from AVM to IM

####Data from jumps boolean
##Initial conditions, AVM states
#      #LD          #SD
ci2 <- c(82806913, 82282625)
ci_vec2 <- sapply(ci2, int2binState, net$genes)

###Ading control state initial state
controles2 <- matrix(data = c(0,0,0,1), nrow = 4, ncol=2)
rownames(controles2) <- net$genes[c(length(net$genes)-3):length(net$genes)]

####Initial conditions
ci_vec2 <- rbind(ci_vec2[1:23,], controles2)

initc2 <- apply(ci_vec2, FUN = bin2intState, MARGIN = 2)

tabla2 <- cbind(initc2,which(net$genes=="CFCA"))
##################################################################
####Third jump from IM to FM

####Data from jumps boolean
##Initial conditions, IM states
#      #LD          #SD
ci3 <- c(75164577, 74640033)
ci_vec3 <- sapply(ci3, int2binState, net$genes)

###Ading control state initial state
controles3 <- matrix(data = c(0,0,0,0), nrow = 4, ncol=2)
rownames(controles3) <- net$genes[c(length(net$genes)-3):length(net$genes)]

####Initial conditions
ci_vec3 <- rbind(ci_vec3[1:23,], controles3)

initc3 <- apply(ci_vec3, FUN = bin2intState, MARGIN = 2)

tabla3 <- cbind(initc3,which(net$genes=="CPNY"))

######################################################
#plot all genes separate in b&w 

plotFS <- function(vector){
  iniAttr <- vector[1] ##Atractor inicial
  ilabel<- labelState(iniAttr,node.names = net$genes, label.rules = labels.rules)
  geneK <- as.numeric(vector[2]) ###Gene a perturbar su tasa de decaimiento
  geneM <- function(geneK=geneK){
    if(net$genes[geneK]=="CVER") return("VER") 
    if(net$genes[geneK]=="CAGE") return("AGE") 
    if(net$genes[geneK]=="CPNY") return("PNY")
    if(net$genes[geneK]=="CFCA") return("FCA")
    else(return("not input"))}
  ###Modelo ODEs y reslucion numerica
  net.ode <- net.ode #booleanToODE(net)#, logic = "Probabilistic")
  net.ode$parameters[geneK+2]= 1.0 ##PARAMETERVECTOR, FIRST 2 ARE  h=10, w=0.5 and THEN decay rate = 1 FOR ALL GENES.
  odeOut21 <- ode(func = net.ode$func, parms = net.ode$parameters, y = int2binState(iniAttr, net$genes), times = seq(0, 10, 0.1))
  ##Name odeout final state
  out <- round(odeOut21[101,2:(ncol(odeOut21))])
  fstate<- Comparar(out)
  flabel <- labelState(fstate, node.names = net$genes, label.rules= labels.rules)
  ###Plot
  odeOut21A <- melt(odeOut21[,2:24], id.vars=odeOut21[,1])
  names(odeOut21A) <- c("Time", "Gene", "value")
  #titulo1 <- paste("Initial ",ilabel," state: ",vector[1],".","\nFinal ",flabel, " state: ",fstate,".", sep = "")
  titulo <- paste("From ",ilabel," to ",flabel,"\nPhotoperiod: ",photoperiod(fstate),"\nInput modified: ",geneM(geneK), sep = "")
  p1<-ggplot(data = odeOut21A, aes(x=odeOut21A$Time ,y = odeOut21A$value, factor=Gene))+ 
    geom_line(position=position_dodge(width=1))+xlab("Time") +ylab("Node State") +
    ggtitle(titulo)+
    facet_grid(Gene ~.)+
    scale_y_continuous(breaks = c(0,1)) +
    scale_x_continuous(breaks = c(0,50,100)) +
    theme(axis.text=element_text(size = rel(.5)),legend.position = "none",axis.title=element_text(size=rel(.7)), plot.title = element_text(size = rel(.7)), strip.text.y = element_text(angle = 0, size = rel(.7)))
  p1
}

#######Image b/w all genes separated, ld and sd 
##all initial conditions data for plotFS
tabla4<- rbind(tabla, tabla2, tabla3)

allS <- apply(tabla4,FUN = plotFS,MARGIN = 1)
##ld
plot_grid(allS[[1]],allS[[2]],allS[[5]],  allS[[7]], ncol = 4)
##sd
plot_grid(allS[[3]],allS[[4]],allS[[6]],allS[[8]], ncol= 4)

#Caro Chavez 2020 Abril
