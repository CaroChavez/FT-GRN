library("BoolNet")
library("BoolNetPerturb")
library("methods")

###Robustness analysis for FT-GRN, code to create data for Fig 5###
#######If only plotting data use plotFig5.R code #####

print(Sys.time())
net<- loadNetwork(file = "ftgrn.txt")
#net2 <- fixGenes(net1, fixIndices = c("VER"), values = c(1))
attr<- getAttractors(net)

wtDF <- attractorToDataframe(attr)

set.seed(38)

######Comparing the number of attractors that are the same in the network and their pertubed copies.
#Test attractor robustness by fliping one line in the true tables of each gene
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
  states = tapply(wtDF$involvedStates,wtDF$label, paste),
  basin = tapply(wtDF$basinSize, wtDF$label, sum),
  percentage = tapply(wtDF$basinSize, wtDF$label, sum)*100/sum(wtDF$basinSize)
)

attractors.by.label

#Look for unlabel attractors

unique(wtDF$label)

# Generate a list of all assignments of n variables with N possible values
allcombn <- function(N,n)
{
  rownum = N^n
  sapply(n:1,function(i)
  {
    rep(seq_len(N),each=N^(i-1),len=rownum)
  })
}

# Get DNF representation of a truth table <truthTable>
# using the gene names in <genes>. 
# If <mode> is "canonical", build a canonical DNF.
# If <mode> is "short", join terms to reduce the DNF
getDNF <- function(truthTable, genes, mode = c("short","canonical"))
{
  if (mode[1] == TRUE)
    mode <- (if (length(genes) <= 12) "short" else "canonical")
  
  mode <- match.arg(mode, c("short","canonical"))
  # check for constant functions
  if (isTRUE(all.equal(truthTable,rep(0,length(truthTable)))))
    return("0")
  else
    if (isTRUE(all.equal(truthTable,rep(1,length(truthTable)))))
      return("1")
  
  # generate truth table
  entries <- allcombn(2,length(genes)) - 1
  colnames(entries) <- genes
  
  
  if (mode == "short")
  {
    # heuristic minimization
    
    # the 1 terms that need to be covered
    uncoveredEntries <- which(truthTable == 1)
    
    # current conjunction list
    conjunctions <- list()  
    
    while (length(uncoveredEntries) > 0)
    {
      # take an uncovered entry and expand it
      currentEntry <- entries[uncoveredEntries[1],]
      
      for (gene in genes)
        # test for each gene whether it can be eliminated from the term
      {
        geneIdx <- which(names(currentEntry) == gene)
        candidate <- currentEntry[-geneIdx]
        condition <- rep(TRUE,length(truthTable))
        for (i in seq_along(candidate))
        {
          condition <- condition & (entries[,names(candidate)[i]] == candidate[i])
        }
        
        if (length(unique(truthTable[condition])) == 1)
          # eliminate gene
          currentEntry <- currentEntry[-geneIdx]
      }
      
      # determine which truth table result entries are now covered
      eliminatedEntries <- rep(TRUE,length(truthTable))
      for (i in seq_along(currentEntry))
      {
        eliminatedEntries <- eliminatedEntries & 
          (entries[,names(currentEntry)[i]] == currentEntry[i])
      }
      uncoveredEntries <- setdiff(uncoveredEntries, which(eliminatedEntries))
      
      # remember conjunction
      conjunctions <- c(conjunctions, list(currentEntry))
    }
    return(paste(paste("(",sapply(conjunctions, function(conj)
    {
      paste(mapply(function(gene, val)
      {
        if (val == 1)
          return(gene)
        else
          return(paste("!",gene,sep=""))
      }, names(conj), conj), collapse=" & ")
    }), ")", sep=""), collapse=" | "))
  }
  else
  {
    # canonical DNF
    conjunctions <- apply(entries[truthTable==1,,drop=FALSE],1,function(conj)
    {
      paste("(",paste(sapply(seq_along(conj),function(lit)
      {
        if (conj[lit])
          genes[lit]
        else
          paste("!",genes[lit],sep="")
      }),collapse=" & "),")",sep="")
    })
    return(paste(conjunctions[conjunctions != ""],collapse = " | "))
  }
}

# Retrieves a string representation of an interaction function by either
# building a DNF (if <readableFunction> is false)
# or returning an unspecific function description 
# (if <readableFunction> is true or "canonical" or "short").
# <truthTable> contains the result column of the interaction's truth table, and
# <genes> contains the names of the involved genes.
getInteractionString <- function(readableFunctions,truthTable,genes)
{
  if (readableFunctions != FALSE)
    getDNF(truthTable,genes, readableFunctions)
  else
  {
    if (isTRUE(all.equal(truthTable,rep(0,length(truthTable)))))
      return("0")
    else
      if (isTRUE(all.equal(truthTable,rep(1,length(truthTable)))))
        return("1")
    else
    {
      truthTable <- sapply(truthTable,function(x)
      {
        if (x == 0)
          "0"
        else
          if (x == 1)
            "1"
        else  
          "*"
      })
      paste("<f(",
            paste(genes,collapse=","),"){",
            paste(truthTable,collapse=""),"}>", sep="")
    }
  }    
}

####################Perturb the function of a single gene, and return a new network
perturbNetworkIdx <- function (network, perturb = c("functions"), method = c("bitflip", 
                                                                             "shuffle"), simplify = (perturb[1] != "functions"), readableFunctions = FALSE, 
                               maxNumBits = 1, numStates = max(1, 2^length(network$genes)/100), perturbIndex = perturbIndex, flipIndices) 
{
  stopifnot(inherits(network, "BooleanNetwork") | inherits(network, 
                                                           "ProbabilisticBooleanNetwork"))
  
  perturbIndex = perturbIndex
  if (length(perturb) == 1 && perturb == "states") {
    warning("perturb=\"states\" is deprecated. Use perturb=\"transitions\" instead!")
    perturb <- "transitions"
  }
  if (inherits(network, "BooleanNetwork")) {
    switch(match.arg(perturb, c("functions")), 
           functions = switch(match.arg(method, c("bitflip", 
                                                  "shuffle")), bitflip = {
                                                    functionIdx = perturbIndex
                                                    flipIndices <- flipIndices
                                                    network$interactions[[functionIdx]]$func[flipIndices] <- as.integer(!network$interactions[[functionIdx]]$func[flipIndices])
                                                    network$interactions[[functionIdx]]$expression <- getInteractionString(readableFunctions, 
                                                                                                                           network$interactions[[functionIdx]]$func, network$genes[network$interactions[[functionIdx]]$input])
                                                  }, shuffle = {
                                                    functionIdx = perturbIndex      
                                                    flipIndices <- flipIndices
                                                    network$interactions[[functionIdx]]$func <- network$interactions[[functionIdx]]$func[flipIndices]
                                                    network$interactions[[functionIdx]]$expression <- getInteractionString(readableFunctions, 
                                                                                                                           network$interactions[[functionIdx]]$func, network$genes[network$interactions[[functionIdx]]$input])
                                                  }, stop("'method' must be one of \"bitflip\",\"shuffle\"")), 
           stop("'perturb' must be one of \"functions\",\"transitions\""))
  }
  if (simplify) 
    network <- simplifyNetwork(network, readableFunctions)
  return(network)
}


#######Create one perturbed network with a change in one of the outputs of the truth tables in the function of each gene (indices)
indices <- c(1,6)

identicosWT <- c()
newAttr <- c()
mismosAttr <- function(indice){
  vector= c(1:length(net$interactions[[indice]]$func))
  print(indice)
  newBS <- c()
  for(i in vector){
    print(c(i))
    netp <- perturbNetworkIdx(net, perturb = "functions", method = "bitflip", maxNumBits = 1, simplify = TRUE, perturbIndex= indice, flipIndices=i)
    attrP <- getAttractors(netp)
    netDF <- attractorToDataframe(attrP)
    rm(attrP)
    idxN <- which((netDF$involvedStates %in% wtDF$involvedStates)==FALSE)
    sumaBS <- sum(netDF$basinSize[idxN])
    newBS <- append(newBS, sumaBS)
    #print(netDF$involvedStates)
    iguales <- identical(netDF$involvedStates,wtDF$involvedStates)
    if(iguales == FALSE){
      #print("not identical")
      allWT <-netDF$involvedStates %in% wtDF$involvedStates
      nuevosAt <- FALSE %in% allWT
      #if(nuevosAt == TRUE) print("new attractors")
      #if(nuevosAt==FALSE) print("only wt but some missing")
      newAttr <- append(newAttr,nuevosAt)
    }
    else {
      #print("identical attractors")
      newAttr <- append(newAttr,!iguales)}
    identicosWT <- append(identicosWT,iguales)
  } 
  return(cbind(identicosWT,newAttr,newBS))
}


###100 RANDOM PERTURBED FUNCTIONS: identicosWT means this perturbed network recovered the same wt attractors; 
#newAttr means this perturbed network recovered aditional new attractors not found in WT
cienRNg<-mismosAttr(c(1)) #dummy example for first node of the network
print(cienRNg)

########################ATENTION  RUN TIME WARNING###################
###Take into consideration the truthTable size for each gene, select only those genes with truth tables smaller than 3000 rows; which leaves AP1 and SOC1 excluded;
##SOC1 and AP1 will be evaluated separetly because it takes to much time (4000 nets ~ 9h computing time)

###################ATENTION Uncoment 5 lines below to run code for all nodes with small truth tables (all except AP1 and SOC1) and save data. This may take a few hours, better if run from command line
#ttSize <- sapply(c(1:23),FUN= function(x) length(net$interactions[[x]]$fun))

#indices <- which(ttSize < 3000) 

#funcionTodos <- lapply(indices,FUN = mismosAttr)

#names(funcionTodos)<- net$genes[indices]

#saveRDS(cienRN, file = "~/Desktop/milRN_SAM130320.rds")


###################ATENTION Uncoment 2 lines below to run code for AP1. This may take a few hours, use command line
#funcionAP1 <- lapply(which(net$genes=="AP1"),FUN = mismosAttr)
#names(funcionAP1) <- c("AP1")

#######################ATENTION Uncoment 2 lines below to run code for SOC1. This may take a few hours, use command line
#funcionSOC1 <- lapply(which(net$genes=="SOC1"),FUN = mismosAttr)
#names(funcionSOC1) <- c("SOC1d") 

#saveRDS(funcionSOC1, file = "~/Desktop/RE_SOC1dTTsam130320.rds") ###1:4000

print(Sys.time())

#Caro Chavez UNAM ABRIL 2020
