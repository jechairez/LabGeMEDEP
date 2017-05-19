
source("GRN_ODEs.r")
source("FUNCIONES.R")

library(BoolNet)
library(deSolve) # If odesolve is not already loaded
library(rootSolve)

RED.BoolNet <- "RedCC.txt"

# READ FOS-GRN Boolean/Logic  model
BoolNet <- loadNetwork(RED.BoolNet)

Attractors.Table <- AttrTable(BoolNet)

write.table(t(Attractors.Table), "Attractor_Table.txt", row.names = FALSE)
Attractors <- read.table("Attractor_Table.txt", header=TRUE)

AttrsLandscape <- Get.States.Attractor.Cluster(BoolNet)

state <- c(APC=1,KRP1=0,CYCA23=0,CDKB11=1,CYCB11=0,MYB3R=1,MYB77=0,E2Fe=0,E2Fc=1,E2Fb=0,E2Fa=0,RBR=1,SCF=0,CYCD31=0)
colnames(Attractors) <- names(state)

########################################################################

BinStateSpace <- AttrsLandscape[[1]]
colnames(BinStateSpace) <- names(state)

AttractorsVec <- apply(Attractors, 1,function(i)  paste(i, collapse=""))
AttrsBasins <- numeric(nrow(BinStateSpace))

#####################################################################################################################
#####################################################################################################################
# Analisis de Bifurcaciones
#####################################################################################################################

# Attri <- 1

AttractorsMat <- Attractors

for(Attri in 1:nrow(AttractorsMat)) {
  InitialState <- as.numeric(AttractorsMat[Attri,])
  names(InitialState) <- colnames(AttractorsMat)
  ActGen <- names(InitialState)[InitialState==1]

  for(i in 1:length(ActGen)) {

	  Kpars <- rep(1,length(state))
	  names(Kpars) <- paste("k", names(state), sep="")
	  Parms <- c( h = 15, Kpars)
	  ParIndex <- which(names(Parms)==paste("k",ActGen[i],sep=""))

	  Ks <- seq(1,5, length=100)
	  AttrSums <- numeric(length(Ks))
  
	  for(j in 1:length(Ks)) {
	
	  Parms[ParIndex] <- Ks[j]

	  AttrSums[j] <- sum(runsteady(y = InitialState, fun = network, parms = Parms, times = c(0, 1e5))$y)
	 
	}
	
	svg(paste(rownames(AttractorsMat)[Attri], "_", ActGen[i],"k", ".svg", sep=""), bg="transparent")  
	plot(Ks[2:length(Ks)], AttrSums[2:length(AttrSums)], pch=20, xlab=paste(ActGen[i], "k"), ylab="Sum(Attractor)", ylim=c(0,10))
	Parms[ParIndex] <- 1
	InitS <- Discretize(matrix(runsteady(y = InitialState, fun = network, parms = Parms, times = c(0, 1e5))$y, 1 ))
	abline(h=sum(InitS), col=4, lwd=2, lty=2)
	Parms[ParIndex] <- Ks[length(Ks)]
	FinalS <- Discretize(matrix(runsteady(y = InitialState, fun = network, parms = Parms, times = c(0, 1e5))$y, 1 ))
	abline(h=sum(FinalS), col=2, lwd=2, lty=2)
	if(paste(FinalS, collapse="")%in%AttractorsVec) { 
	  text(4, 10, paste("Initial Attractor:", rownames(AttractorsMat)[Attri]), col=4)
	  text(4, 9, paste("Final Attractor:", names(AttractorsVec)[which(paste(FinalS, collapse="")==AttractorsVec)]), col=2)
	}
	if(!paste(FinalS, collapse="")%in%AttractorsVec) {
	  text(4, 10, paste("Initial Attractor:", rownames(AttractorsMat)[Attri]), col=4)
	  text(4, 9, paste("Final Attractor:", "Other"), col=2)
	}
	dev.off()
  }		

}

