########################################################################################################################################
########################################################################################################################################
# FUNCTIONS
########################################################################################################################################
########################################################################################################################################


########################################################################################################################################
##Generic 
########################################################################################################################################
Generate.Binary.State.Space <- function(Ngenes) {
  VecMat <- matrix(0, 2^Ngenes, Ngenes)
  for(i in 1:ncol(VecMat)) VecMat[,i] <- rep(c(0,1), each=2^(Ngenes-i))
  return(VecMat)
}
########################################################################################################################################
Calculate.Hamming.Distances <- function(VectorsMat) {
  DISTS <- matrix(0, nrow(VectorsMat), nrow(VectorsMat))
  for(Si in 1:nrow(VectorsMat)) DISTS[Si,] <- sapply(1:nrow(VectorsMat), function(i) sum(stringdist(VectorsMat[Si,], VectorsMat[i,])))
   return(DISTS) 
}
########################################################################################################################################

Discretize <- function(Data) {
  DataDisc <- Data*0
  for(i in 1:nrow(Data)) {
    Model <- kmeans(as.numeric(Data[i,]), 2)
    DataDisc[i,Model$cluster==which.min(as.numeric(Model$centers))] <- 0
    DataDisc[i,Model$cluster==which.max(as.numeric(Model$centers))] <- 1
  }
  return(DataDisc)
}


########################################################################################################################################
##FUNCTIONS READ NETWORK
########################################################################################################################################
GenerateRandomNet <- function(AdjMatrx) {
# Generates a random Network in matrix form using the proportions of characters in the original net as sampling probabilities
  return(matrix(sample(c(-1,0,1), nrow(AdjMatrx)*ncol(AdjMatrx), rep=T, prob=table(as.numeric(AdjMatrx))/length(as.numeric(AdjMatrx))), nrow(AdjMatrx),ncol(AdjMatrx)))
}
########################################################################################################################################
GetGeneNames <- function(Interacs) {
#Extracts a list of Genenames from the network/interactions file
  InteracsList <- lapply(1:length(Interacs), function(i) unlist(strsplit(Interacs[i]," ")))
  GenesIn <- sapply(1:length(InteracsList), function(i) InteracsList[[i]][1])
  GenesTarget <- sapply(1:length(InteracsList), function(i) InteracsList[[i]][3])
  GeneNames <- unique(c(unique(GenesIn), unique(GenesTarget)))
  return(GeneNames)
}
########################################################################################################################################
GenerateAdjMatrix <- function(Interacs) {
  GeneNames <- GetGeneNames(Interacs)
  InteracsList <- lapply(1:length(Interacs), function(i) unlist(strsplit(Interacs[i]," ")))
  AdjMatrx <- matrix(0, length(GeneNames), length(GeneNames))
  
  rownames(AdjMatrx) <- GeneNames
  colnames(AdjMatrx) <- GeneNames

  for(i in 1:length(InteracsList)) {
    if(sum(InteracsList[[i]]=="->")>0) AdjMatrx[which(GeneNames==InteracsList[[i]][1]), which(GeneNames==InteracsList[[i]][3])] <- 1
    if(sum(InteracsList[[i]]=="-|")>0) AdjMatrx[which(GeneNames==InteracsList[[i]][1]), which(GeneNames==InteracsList[[i]][3])] <- -1
  }
  return(AdjMatrx)
}
########################################################################################################################################
########################################################################################################################################

##############################################################################################################
## Boolean Boolnet Functions
##############################################################################################################
AttrTable <- function(Net) {
  #Extracts the attractors in a table form from the attractor object created by BoolNet function "getAttractors"
  attrs <- getAttractors(Net)  
  Attrs <- print(plotAttractors(attrs))
  Attrs <- Attrs[[1]]
  return(Attrs)
}

########################################################
# Calculate KO Mutants
########################################################
Get.KOs.States.Attractor.Cluster <- function(Network) {

  States.Attractor.Cluster <- list()
  
  for(i in 1:length(Network$genes)) {
    KOi <- fixGenes(Network, Network$genes[i], 0)  
    States.Attractor.Cluster[[i]] <- Get.States.Attractor.Cluster(KOi)
  }
  names(States.Attractor.Cluster) <- Network$genes 
  return(States.Attractor.Cluster)
}  
##############################################################################################################
# Stochastic
##############################################################################################################
# Calculate Attractor's Transition Probability Matrix
##############################################################################################################
AttrTable <- function(Net) {
  #Extracts the attractors in a table form from the attractor object created by BoolNet function "getAttractors"
  attrs <- getAttractors(Net)  
  Attrs <- print(plotAttractors(attrs))
  Attrs <- Attrs[[1]]
  return(Attrs)
}

Get.States.Attractor.Cluster <- function(Network) {
  attrs <- getAttractors(Network)
  TransTable <- getTransitionTable(attrs)
  
  Inits <- TransTable[[1]]
  for(i in 2:length(Network$genes)) Inits <- cbind(Inits, TransTable[[i]])

  AttractsClusterVector <- TransTable$attractorAssignment

  return(list(Inits, AttractsClusterVector))
}


Generate.Transition.Matrix <- function(Network, Error.Prob, Reps) {
  AttrsN <- as.character(unlist(apply(AttrTable(Network),2, function(i) paste(i, collapse=""))))

  Attr.Lands <- Get.States.Attractor.Cluster(Network)
  State.Space <- Attr.Lands[[1]]
  Attr.Vector <- Attr.Lands[[2]]
  State.Space.Collapsed <- apply(State.Space, 1,  function(i) paste(i, collapse=""))

  Transition.Matrix <- matrix(0, length(table(Attr.Vector)), length(table(Attr.Vector)))
  
  #colnames(Transition.Matrix) <- AttrsN
  #rownames(Transition.Matrix) <- AttrsN
  
#i <- 1
#
  for(i in 1:dim(State.Space)[1]) {

    TempTrans <- stateTransition(Network, State.Space[i,])

    for(Rep in 1:Reps) {
    #print("############################################################")
    #print(Rep)
    #print("############################################################")
    Luckies <- sample(c(0,1), dim(State.Space)[2], replace=TRUE, prob=c(1-Error.Prob, Error.Prob))
      if(sum(Luckies)>0) {
	Temp <- TempTrans
	if(sum(TempTrans[which(Luckies==1)]==0)>0) Temp[which(TempTrans[which(Luckies==1)]==0)] <- 1  
	if(sum(TempTrans[which(Luckies==1)]==1)>0) Temp[which(TempTrans[which(Luckies==1)]==1)] <- 0
	Temp <- paste(Temp, collapse="")
	A.cual <- which(State.Space.Collapsed == Temp)
	Transition.Matrix[Attr.Vector[i], Attr.Vector[A.cual]] <- Transition.Matrix[Attr.Vector[i], Attr.Vector[A.cual]] + 1
      }
    
      if(sum(Luckies)==0) Transition.Matrix[Attr.Vector[i], Attr.Vector[i]] <- Transition.Matrix[Attr.Vector[i], Attr.Vector[i]] + 1  
    }
   # print("############################################################")
    #print("############################################################")
    print("############################################################")
    print(paste("State",i))
    print("############################################################")
    #print("############################################################")
    #print("############################################################")
  }
  
#   rownames(Transition.Matrix) <- names(table(Attr.Vector))
#   colnames(Transition.Matrix) <- names(table(Attr.Vector))
  return(Transition.Matrix)
}


Generate.Prob.Matrix <- function(TransProbs) {
  ProbMat <- TransProbs*0
  for(i in 1:nrow(TransProbs)) ProbMat[i,] <- TransProbs[i,]/sum(TransProbs[i,])
  return(as.matrix(ProbMat))
}

# function for sampling a finite Markov chain

rfmc<-function(n,P,pi0) {
        v=vector("numeric",n)
        r=length(pi0)
        v[1]=sample(r,1,prob=pi0)
        for (i in 2:n) {
                v[i]=sample(r,1,prob=P[v[i-1],])
        }
        ts(v)
}

##############################################################################################################
##############################################################################################################

##############################################################################################################
## Boolean NO LOGIC Functions
##############################################################################################################


One.Step.State.Transition <- function(AdjMat, State.t, Bs) {
  Total.Inputs.Vec <- sapply(1:ncol(AdjMat), function(i) sum(AdjMat[,i]*State.t) + Bs[i])  
  
  State <- rep(0,ncol(AdjMat))
  names(State) <- colnames(AdjMat)

  if(sum(Total.Inputs.Vec>0)!=0) State[Total.Inputs.Vec>0] <- 1
  if(sum(Total.Inputs.Vec<0)!=0) State[Total.Inputs.Vec<0] <- 0
  if(sum(Total.Inputs.Vec==0)!=0) State[Total.Inputs.Vec==0] <- State.Vector[Total.Inputs.Vec==0]
  
  return(State)
}

Generate.OneStepTransitions <- function(State.Space, AdjMat, Bs) {
  Space.Transitions.Table <- list(character(nrow(State.Space)), character(nrow(State.Space)))

  for(i in 1:nrow(State.Space)) {
    #i <- 1
    Space.Transitions.Table[[1]][i] <- paste(State.Space[i,], collapse="")
    Space.Transitions.Table[[2]][i] <- paste(One.Step.State.Transition(AdjMat, State.Space[i,], Bs), collapse="")
  }
  return(Space.Transitions.Table)
}

Forward.Set <- function(State0, Transitions) {
  State0 
  FSet <- State0

  Next <- Transitions[[2]][which(State0==Transitions[[1]])]
  while(sum(Next %in% FSet)==0){
    FSet <- c(FSet, Next)
    State0 <- Next
    Next <- Transitions[[2]][which(State0==Transitions[[1]])]
  }
  return(FSet)
}


Calcula.Attractors <- function(Transitions.Table) {
  F.SETS <- lapply(Transitions.Table[[1]], function(i) Forward.Set(i, Transitions.Table))
  B.SETS <- lapply(1:length(Transitions.Table[[1]]), function(j) Transitions.Table[[1]][sapply(F.SETS, function(i) sum(i %in% Transitions.Table[[1]][j])==1)])

  Cual <- sapply(1:length(F.SETS), function(i) sum(F.SETS[[i]]%in%B.SETS[[i]])==length(F.SETS[[i]]))
  AttrElements <- Transitions.Table[[1]][which(Cual)]
  Attractors <- lapply(which(Cual), function(i)  F.SETS[[i]])
  OUT <- list()
  OUT[[1]] <- F.SETS
  OUT[[2]] <- B.SETS
  OUT[[3]] <- Attractors
  return(OUT)
}

##########################
## Probabilistic Landscape



## Transition Probabilities

Generate.OneStepTransitionsP <- function(State.Space, AdjMat, Bs) {
  Space.Transitions.Table <- list(character(nrow(State.Space)), character(nrow(State.Space)))

  for(i in 1:nrow(State.Space)) {
    #i <- 1
    Space.Transitions.Table[[1]][i] <- paste(State.Space[i,], collapse="")
    Space.Transitions.Table[[2]][i] <- paste(One.Step.State.Transition(AdjMat, State.Space[i,], Bs), collapse="")
  }
  return(Space.Transitions.Table)
}


One.Step.State.TransitionP <- function(AdjMat, State.t, State.t1, Bs) {
  Total.Inputs.Vec <- sapply(1:ncol(AdjMat), function(i) sum(AdjMat[,i]*State.t) + Bs[i])  
  
  State <- rep(0,ncol(AdjMat))
  names(State) <- colnames(AdjMat)

  if(sum(Total.Inputs.Vec>0)!=0) State[Total.Inputs.Vec>0] <- 0.5 + 0.5*tanh(Miu*Total.Inputs.Vec)
  if(sum(Total.Inputs.Vec<0)!=0) State[Total.Inputs.Vec<0] <- 0.5 - 0.5*tanh(Miu*Total.Inputs.Vec)
  if(sum(Total.Inputs.Vec==0)!=0) State[Total.Inputs.Vec==0] <- 1-d
 
  return(prod(State))
}

#   Current <- State.Space[6824,]
#   Next  <- State.Space[2535,]
  
CalculateTransitionProb <- function(Current, Next, d, Miu, AdjMat) {

  Total.Inputs.Vec <- sapply(1:ncol(AdjMat), function(i) sum(AdjMat[,i]*Current) + Bs[i])
  
#   ProbperNode <- rep(0,ncol(AdjMat))
#   names(ProbperNode) <- colnames(AdjMat)
#   
#   if(Total.Inputs.Vec[i]==0){
#     if(Current[i]==Next[i]) ProbperNode[i] <- 1-d
#     if(Current[i]!=Next[i]) ProbperNode[i] <- d  
#   }    
#   
#   if(Next[i]==0) ProbperNode[i] <- 0.5 - 0.5*tanh(Miu*Total.Inputs.Vec[i])
#   if(Next[i]==1) ProbperNode[i] <- 0.5 + 0.5*tanh(Miu*Total.Inputs.Vec[i])

 ProbperNode <- sapply(1:ncol(AdjMat),
  
    function(i) {
      if(Total.Inputs.Vec[i]==0){
	return(1-d)
	#if(Current[i]==Next[i]) return(1-d)
	#if(Current[i]!=Next[i]) return(d)
      }    
  
      if(Next[i]==0) return(0.5 - 0.5*tanh(Miu*Total.Inputs.Vec[i]))
      if(Next[i]==1) return(0.5 + 0.5*tanh(Miu*Total.Inputs.Vec[i]))
    }
  )
  
  return(prod(ProbperNode))
  
}

Simulate.Master.Equation <- function(RED, TransProbMatrix, Timef, dt) {
  ti <- 0

  Times <- ti  

  ProbVectorT0 <- rep(1/2^ncol(RED), ncol(TransProbMatrix))
  names(ProbVectorT0) <- colnames(TransProbMatrix)

  ProbVectorTi <- ProbVectorT0

  for(Ts in 1:Timef) {
    ProbVectorTemp <- ProbVectorTi
    for(i in 1:ncol(TransProbMatrix)) {
      Diff <- sum(ProbVectorTi*TransProbMatrix[,i])*dt - sum(ProbVectorTi[i]*TransProbMatrix[i,])*dt
      ProbVectorTemp[i] <- ProbVectorTi[i] + Diff
    }
    ProbVectorTi <- ProbVectorTemp
  
    ProbVectorT0 <- rbind(ProbVectorT0, ProbVectorTi)
    ti <- ti+dt
    Times <- c(Times, ti)
  }

  return(list(Times, ProbVectorT0))
}


##############################################################################################################
## ODEs Hill Model 
##############################################################################################################

HillActivation <- function(xi, alfa, s, n) {
  return((alfa*xi^n)/(s^n + xi^n))
}

HillInhibit <- function(xi, beta, s, n) {
  return((beta*s^n)/(s^n + xi^n))
}

Generate.Fi <- function(AdjMat, ind) {
  if(sum(AdjMat[,ind])!=0) {
    if(sum(AdjMat[,ind]==1)!=0) {
      PosNames <- names(which(AdjMat[,ind]==1))
      Activs <- paste("HillActivation(c(", paste(PosNames, collapse=","), ")", sep="", ",", paste(c("alfa", "s", "n"), collapse=", "),")" )
    }
    if(sum(AdjMat[,ind]==1)==0) Activs <- "0"
    
    if(sum(AdjMat[,ind]==-1)!=0) {
      NegNames <- names(which(AdjMat[,ind]==-1))
      Negs <- paste("HillInhibit(c(", paste(NegNames, collapse=","), ")", sep="", ",", paste(c("beta", "s", "n"), collapse=", "),")" )
    }
    if(sum(AdjMat[,ind]==-1)==0) Negs <- "0"

    Term1<- paste("sum(",Activs,")", sep="")
    Term2 <- paste("sum(",Negs,")", sep="")
    Term3 <- paste("k*", colnames(AdjMat)[ind], sep="")

    F.i <- paste(Term1, "+", Term2, "-", Term3)
  }

  if(sum(AdjMat[,ind])==0) {
    Term1 <- colnames(AdjMat)[ind]
    Term2 <- paste("k*", colnames(AdjMat)[ind], sep="")
    F.i <- paste(Term1, "-", Term2)
  }
  return(F.i)
}

GeneraModeloODEsHill <- function(AdjMat, Path, FileName) { 
  Head <- paste(paste("d", colnames(AdjMat), sep=""), "=")
  Body <- sapply(1:ncol(AdjMat), function(i) paste(Head[i], Generate.Fi(AdjMat, ind=i)))
  Tail <- paste("dX <- c(", paste("d", colnames(AdjMat), sep="", collapse=","), ")")
  ODEs <- c("RedODEs <- function(t, x, parms) {" , "with(as.list(c(x, parms)),{", Body, Tail, "list(dX)", "})", "}")
  write(ODEs, paste(Path, FileName, sep=""))
}


##############################################################################################################
##############################################################################################################
##############################################################################################################