############################
#ESTE PROGRAMA RESUELVE EL MODELO CONTINUO DESDE N CONDICIONES ALEATORIAS (VARIABLE REPETICIONES) Y EVALUA SI LOS ESTADOS ESTACIONARIOS CORRESPONDEN A LOS ATRACTORES FIJOS DEL MODELO BOOLEANO. LOS PARAMETROS H (50) Y LAMBDA (1) SON IGUALES PARA TODOS LOS NODOS.
#este programa se hizo para GHRN
#UNKNOWN1 Y 2 SON C. PROV TD2 Y C. PROV TD3 EN García-Gómez, Mónica L., Eugenio Azpeitia, and Elena R. Álvarez-Buylla. "A dynamic genetic-hormonal regulatory network model explains multiple cellular behaviors of the root apical meristem of Arabidopsis thalian." PLOS Computational Biology 13.4 (2017): e1005488.
############################
#ESTA PARTE DECLARA SISTEMA CONTINUO
library(deSolve)
parameters <- c(h = 50,lambda = 1)
network<-function(t, state, parameters) {
	with(as.list(c(state, parameters)),{
	     #node inputs
		 w_CK = max(min(PHB,1-ARFR),1-SHR)
		 w_ARR1 = min(CK,1-SCR)
		 w_SHY2 = min(1-AUX,ARR1)
		 w_AUXIAAR = 1-AUX
		 w_ARFR = 1-AUXIAAR
		 w_ARF10 = min(1-AUXIAAR,1-min(JKD,SHR))
		 w_ARF5 = min(min(1-SHY2,1-AUXIAAR),1-min(SHR,MGP))
		 w_AUX = max(AUX,WOX5)
		 w_SCR = min(SHR,SCR,JKD)
		 w_SHR = max(SHR,min(SCR,JKD))
		 w_MIR165 = max(1-PHB,min(SCR,SHR,1-CK))
		 w_PHB = 1-MIR165
		 w_JKD = min(1-PHB,SCR,SHR)
		 w_MGP = min(1-WOX5,SCR,SHR)
		 w_WOX5 = min(1-ARF10,ARF5,1-CLE40)
		 w_CLE40 = 1-SHR
						
	     #rates of change
		 dCK <- ((-exp(0.5*h)+exp(-h*(w_CK)))/((1-exp(0.5*h))*(1+exp(-h*(w_CK-0.5)))))-(lambda*CK)
		 dARR1 <- ((-exp(0.5*h)+exp(-h*(w_ARR1)))/((1-exp(0.5*h))*(1+exp(-h*(w_ARR1-0.5)))))-(lambda*ARR1)
		 dSHY2 <- ((-exp(0.5*h)+exp(-h*(w_SHY2)))/((1-exp(0.5*h))*(1+exp(-h*(w_SHY2-0.5)))))-(lambda*SHY2)
		 dAUXIAAR <- ((-exp(0.5*h)+exp(-h*(w_AUXIAAR)))/((1-exp(0.5*h))*(1+exp(-h*(w_AUXIAAR-0.5)))))-(lambda*AUXIAAR)
		 dARFR <- ((-exp(0.5*h)+exp(-h*(w_ARFR)))/((1-exp(0.5*h))*(1+exp(-h*(w_ARFR-0.5)))))-(lambda*ARFR)
		 dARF10 <- ((-exp(0.5*h)+exp(-h*(w_ARF10)))/((1-exp(0.5*h))*(1+exp(-h*(w_ARF10-0.5)))))-(lambda*ARF10)
		 dARF5 <- ((-exp(0.5*h)+exp(-h*(w_ARF5)))/((1-exp(0.5*h))*(1+exp(-h*(w_ARF5-0.5)))))-(lambda*ARF5)
		 dAUX <- ((-exp(0.5*h)+exp(-h*(w_AUX)))/((1-exp(0.5*h))*(1+exp(-h*(w_AUX-0.5)))))-(lambda*AUX)
		 dSCR <- ((-exp(0.5*h)+exp(-h*(w_SCR)))/((1-exp(0.5*h))*(1+exp(-h*(w_SCR-0.5)))))-(lambda*SCR)
		 dSHR <- ((-exp(0.5*h)+exp(-h*(w_SHR)))/((1-exp(0.5*h))*(1+exp(-h*(w_SHR-0.5)))))-(lambda*SHR)
		 dMIR165 <- ((-exp(0.5*h)+exp(-h*(w_MIR165)))/((1-exp(0.5*h))*(1+exp(-h*(w_MIR165-0.5)))))-(lambda*MIR165)
		 dPHB <- ((-exp(0.5*h)+exp(-h*(w_PHB)))/((1-exp(0.5*h))*(1+exp(-h*(w_PHB-0.5)))))-(lambda*PHB)
		 dJKD <- ((-exp(0.5*h)+exp(-h*(w_JKD)))/((1-exp(0.5*h))*(1+exp(-h*(w_JKD-0.5)))))-(lambda*JKD)
		 dMGP <- ((-exp(0.5*h)+exp(-h*(w_MGP)))/((1-exp(0.5*h))*(1+exp(-h*(w_MGP-0.5)))))-(lambda*MGP)
		 dWOX5 <- ((-exp(0.5*h)+exp(-h*(w_WOX5)))/((1-exp(0.5*h))*(1+exp(-h*(w_WOX5-0.5)))))-(lambda*WOX5)
		 dCLE40 <- ((-exp(0.5*h)+exp(-h*(w_CLE40)))/((1-exp(0.5*h))*(1+exp(-h*(w_CLE40-0.5)))))-(lambda*CLE40)
	      
		 # return the rate of change
		 list(c(dCK, dARR1, dSHY2, dAUXIAAR, dARFR, dARF10, dARF5, dAUX, dSCR, dSHR, dMIR165, dPHB, dJKD, dMGP, dWOX5, dCLE40))
	})
}
############################
#ESTA PARTE DECLARA MATRIZ CON CONFIGURACIONES DE ACTIVIDAD DE LOS ATRACTORES DEL MODELO BOOLEANO
ATTRACTORS=matrix(nrow=11,ncol=16)
#ORDEN DE LOS NODOS: CK ARR1 SHY2 AUXIAA ARF ARF10 ARF5 AUX SCR SHR MIR PHB JKG MGP WOX CLE
QC<-c(0,0,0,0,1,0,1,1,1,1,1,0,1,0,1,0)
ENDODERMISPD<-c(0,0,0,0,1,0,0,1,1,1,1,0,1,1,0,0)
ENDODERMISTD<-c(0,0,0,1,0,0,0,0,1,1,1,0,1,1,0,0)
PERIPHERALPD<-c(0,0,0,0,1,1,1,1,0,1,1,0,0,0,0,0)
PERIPHERALTD<-c(0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0)
CENTRALPD<-c(0,0,0,0,1,1,1,1,0,1,0,1,0,0,0,0)
CENTRALTD<-c(1,1,1,1,0,0,0,0,0,1,0,1,0,0,0,0)
ROOTCAP1<-c(1,1,1,1,0,0,0,0,0,0,1,0,0,0,0,1)
ROOTCAP2<-c(1,1,0,0,1,1,1,1,0,0,1,0,0,0,0,1)
UNKNOWN2<-c(1,1,0,0,1,1,1,1,0,0,0,1,0,0,0,1)
UNKNOWN1<-c(1,1,1,1,0,0,0,0,0,0,0,1,0,0,0,1)
ATTRACTORS[1,]<-QC
ATTRACTORS[2,]<-ENDODERMISPD
ATTRACTORS[3,]<-ENDODERMISTD
ATTRACTORS[4,]<-PERIPHERALPD
ATTRACTORS[5,]<-PERIPHERALTD
ATTRACTORS[6,]<-CENTRALPD
ATTRACTORS[7,]<-CENTRALTD
ATTRACTORS[8,]<-ROOTCAP1
ATTRACTORS[9,]<-ROOTCAP2
ATTRACTORS[10,]<-UNKNOWN2
ATTRACTORS[11,]<-UNKNOWN1

############################
#ESTE CODIGO  RESUELVE MODELO CONTINUO DESDE N CONDICIONES ALEATORIAS (VARIABLE REPETICIONES), Y GUARDA VALORES EN MATRIZ B (>0.9=1, <0.1=0). 
repeticiones=100000
B=matrix(nrow=repeticiones,ncol=16)

for (i in 1:repeticiones){
	state<-c(CK=runif(1,0,1),ARR1=runif(1,0,1),SHY2=runif(1,0,1),AUXIAAR=runif(1,0,1),ARFR=runif(1,0,1),ARF10=runif(1,0,1),ARF5=runif(1,0,1),AUX=runif(1,0,1),SCR=runif(1,0,1),SHR=runif(1,0,1),MIR165=runif(1,0,1),PHB=runif(1,0,1),JKD=runif(1,0,1),MGP=runif(1,0,1),WOX5=runif(1,0,1),CLE40=runif(1,0,1)) 

	times <-seq(0,100,by=0.1)         #100/0.1=1000. De ahi salen los numeros de este y los sig renglones. 
   	o <- ode(y=state,times=times,func=network,parms=parameters)
   	for (j in 1:16){
   		if(o[1000,j+1]>0.9){B[i,j]<-1}
   		if(o[1000,j+1]<0.1){B[i,j]<-0}
   	}
}
 
############################  	
#ESTE CODIGO COMPARA CADA CONFIGURACIONES GUARDADA EN B CON LAS CONFIGURACIONES GUARDADAS EN ATTRACTORS. SI SON LA MISMA CONFIGURACIONES SE VA CONTANDO EN MATRIZ COUNTER. attractores en counter tienen mismo orden que ATTRACTORS
counter=matrix(c(0),nrow=1,ncol=11)
#newattractors=matrix(nrwo=100,ncol=11) ** PARA QUE EL CODIGO QUEDARA CHIDO FALTARIA GUARDAR Y CONTAR LOS NUEVOS ATRACTORES QUE SALGAN. PERO EN GHRH NO ES NECESARIO.


for (x in 1:repeticiones){
	for(i in 1:11){ #comparo cada attractor
		cohen<-0
		for(j in 1:16){   #comparo cada nodo
			if(ATTRACTORS[i,j]!=B[x,j]){cohen<-cohen+1}	
		}
			if(cohen==0){counter[1,i]=counter[1,i]+1} #si variable es 0 entonces son iguales, podría poner break...  
	}
}
counter #si sum(counter)==repeticiones entonces todas las condiciones aleatorias convergen a uno de los 11 atractores fijos declarados en ATTRACTORS