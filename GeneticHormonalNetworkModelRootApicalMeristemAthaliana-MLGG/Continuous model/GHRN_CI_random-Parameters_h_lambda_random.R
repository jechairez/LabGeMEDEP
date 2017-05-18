#ESTE CODIGO RESUELVE LA RED GHRN. SE EXPLORAN N (variable repet_param) SETS DE PARAMETROS ALEATORIOS (donde h es un valor de 10-50, y lambda es un valor de 0.5-1) DESDE N CONDICIONES ALEATORIAS (variable repeticiones). SE ANALIZA SI ATRACTORES RECUPERADOS SON LOS 11 PUNTOS FIJOS DEL SISTEMA BOOLEANO. EL ARCHIVO OUTPUT.TXT ARROJA LOS SETS DE PARAMETROS Y LAS CI QUE CONVERGEN A LOS 11 ATRACTORES. SI ESTE NUMERO ES IGUAL QUE repeticiones ENTONCES TODAS LAS CONDICIONES ALCANZAN UNO DE LOS 11 ATRACTORS. 
#García-Gómez, Mónica L., Eugenio Azpeitia, and Elena R. Álvarez-Buylla. "A dynamic genetic-hormonal regulatory network model explains multiple cellular behaviors of the root apical meristem of Arabidopsis thalian." PLOS Computational Biology 13.4 (2017): e1005488.

library(deSolve)
ATTRACTORS=matrix(nrow=11,ncol=16)
#CK ARR1 SHY2 AUXIAA ARF ARF10 ARF5 AUX SCR SHR MIR PHB JKG MGP WOX CLE
QC<-c(0,0,0,0,1,0,1,1,1,1,1,0,1,0,1,0)
ENDODERMISPD<-c(0,0,0,0,1,0,0,1,1,1,1,0,1,1,0,0)
ENDODERMISTD<-c(0,0,0,1,0,0,0,0,1,1,1,0,1,1,0,0)
PERIPHERALPD<-c(0,0,0,0,1,1,1,1,0,1,1,0,0,0,0,0)
PERIPHERALTD<-c(0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0)
CENTRALPD<-c(0,0,0,0,1,1,1,1,0,1,0,1,0,0,0,0)
CENTRALTD<-c(1,1,1,1,0,0,0,0,0,1,0,1,0,0,0,0)
ROOTCAP1<-c(1,1,1,1,0,0,0,0,0,0,1,0,0,0,0,1)
ROOTCAP2<-c(1,1,0,0,1,1,1,1,0,0,1,0,0,0,0,1)
UNKNOWN1<-c(1,1,0,0,1,1,1,1,0,0,0,1,0,0,0,1)
UNKNOWN2<-c(1,1,1,1,0,0,0,0,0,0,0,1,0,0,0,1)
ATTRACTORS[1,]<-QC
ATTRACTORS[2,]<-ENDODERMISPD
ATTRACTORS[3,]<-ENDODERMISTD
ATTRACTORS[4,]<-PERIPHERALPD
ATTRACTORS[5,]<-PERIPHERALTD
ATTRACTORS[6,]<-CENTRALPD
ATTRACTORS[7,]<-CENTRALTD
ATTRACTORS[8,]<-ROOTCAP1
ATTRACTORS[9,]<-ROOTCAP2
ATTRACTORS[10,]<-UNKNOWN1
ATTRACTORS[11,]<-UNKNOWN2

repet_param=100
for (i in 1:repet_param){
counter=matrix(c(0),nrow=1,ncol=11)
parameters <- c(h_CK=runif(1,10,50),h_ARR1=runif(1,10,50),h_SHY2=runif(1,10,50),h_AUXIAAR=runif(1,10,50),h_ARFR=runif(1,10,50),h_ARF10=runif(1,10,50),h_ARF5=runif(1,10,50),h_AUX=runif(1,10,50),h_SCR=runif(1,10,50),h_SHR=runif(1,10,50),h_MIR165=runif(1,10,50),h_PHB=runif(1,10,50),h_JKD=runif(1,10,50),h_MGP=runif(1,10,50),h_WOX5=runif(1,10,50),h_CLE40=runif(1,10,50),lambda_CK=runif(1,0.5,1),lambda_ARR1=runif(1,0.5,1),lambda_SHY2=runif(1,0.5,1),lambda_AUXIAAR=runif(1,0.5,1),lambda_ARFR=runif(1,0.5,1),lambda_ARF10=runif(1,0.5,1),lambda_ARF5=runif(1,0.5,1),lambda_AUX=runif(1,0.5,1),lambda_SCR=runif(1,0.5,1),lambda_SHR=runif(1,0.5,1),lambda_MIR165=runif(1,0.5,1),lambda_PHB=runif(1,0.5,1),lambda_JKD=runif(1,0.5,1),lambda_MGP=runif(1,0.5,1),lambda_WOX5=runif(1,0.5,1),lambda_CLE40=runif(1,0.5,1))
	write("parameters", file = "OUTPUT.txt",append = TRUE, sep = " ")
	write(parameters, file = "OUTPUT.txt",ncolumns = if(is.character(x)) 1 else 11,append = TRUE, sep = " ")

	network<-function(t, state, parameters) {
	with(as.list(c(state, parameters)),{
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
		 
				
		 dCK <- ((-exp(0.5*h_CK)+exp(-h_CK*(w_CK)))/((1-exp(0.5*h_CK))*(1+exp(-h_CK*(w_CK-0.5)))))-(lambda_CK*CK)
		 dARR1 <- ((-exp(0.5*h_ARR1)+exp(-h_ARR1*(w_ARR1)))/((1-exp(0.5*h_ARR1))*(1+exp(-h_ARR1*(w_ARR1-0.5)))))-(lambda_ARR1*ARR1)
		 dSHY2 <- ((-exp(0.5*h_SHY2)+exp(-h_SHY2*(w_SHY2)))/((1-exp(0.5*h_SHY2))*(1+exp(-h_SHY2*(w_SHY2-0.5)))))-(lambda_SHY2*SHY2)
		 dAUXIAAR <- ((-exp(0.5* h_AUXIAAR)+exp(- h_AUXIAAR*(w_AUXIAAR)))/((1-exp(0.5* h_AUXIAAR))*(1+exp(- h_AUXIAAR*(w_AUXIAAR-0.5)))))-(lambda_AUXIAAR*AUXIAAR)
		 dARFR <- ((-exp(0.5*h_ARFR)+exp(-h_ARFR*(w_ARFR)))/((1-exp(0.5*h_ARFR))*(1+exp(-h_ARFR*(w_ARFR-0.5)))))-(lambda_ARFR*ARFR)
		 dARF10 <- ((-exp(0.5*h_ARF10)+exp(-h_ARF10*(w_ARF10)))/((1-exp(0.5*h_ARF10))*(1+exp(-h_ARF10*(w_ARF10-0.5)))))-(lambda_ARF10*ARF10)
		 dARF5 <- ((-exp(0.5*h_ARF5)+exp(-h_ARF5*(w_ARF5)))/((1-exp(0.5*h_ARF5))*(1+exp(-h_ARF5*(w_ARF5-0.5)))))-(lambda_ARF5*ARF5)
		 dAUX <- ((-exp(0.5*h_AUX)+exp(-h_AUX*(w_AUX)))/((1-exp(0.5*h_AUX))*(1+exp(-h_AUX*(w_AUX-0.5)))))-(lambda_AUX*AUX)
		 dSCR <- ((-exp(0.5*h_SCR)+exp(-h_SCR*(w_SCR)))/((1-exp(0.5*h_SCR))*(1+exp(-h_SCR*(w_SCR-0.5)))))-(lambda_SCR*SCR)
		 dSHR <- ((-exp(0.5*h_SHR)+exp(-h_SHR*(w_SHR)))/((1-exp(0.5*h_SHR))*(1+exp(-h_SHR*(w_SHR-0.5)))))-(lambda_SHR*SHR)
		 dMIR165 <- ((-exp(0.5* h_MIR165)+exp(-h_MIR165*(w_MIR165)))/((1-exp(0.5* h_MIR165))*(1+exp(-h_MIR165*(w_MIR165-0.5)))))-(lambda_MIR165*MIR165)
		 dPHB <- ((-exp(0.5*h_PHB)+exp(-h_PHB*(w_PHB)))/((1-exp(0.5*h_PHB))*(1+exp(-h_PHB*(w_PHB-0.5)))))-(lambda_PHB*PHB)
		 dJKD <- ((-exp(0.5*h_JKD)+exp(-h_JKD*(w_JKD)))/((1-exp(0.5*h_JKD))*(1+exp(-h_JKD*(w_JKD-0.5)))))-(lambda_JKD*JKD)
		 dMGP <- ((-exp(0.5*h_MGP)+exp(-h_MGP*(w_MGP)))/((1-exp(0.5*h_MGP))*(1+exp(-h_MGP*(w_MGP-0.5)))))-(lambda_MGP*MGP)
		 dWOX5 <- ((-exp(0.5*h_WOX5)+exp(-h_WOX5*(w_WOX5)))/((1-exp(0.5*h_WOX5))*(1+exp(-h_WOX5*(w_WOX5-0.5)))))-(lambda_WOX5*WOX5)
		 dCLE40 <- ((-exp(0.5*h_CLE40)+exp(-h_CLE40*(w_CLE40)))/((1-exp(0.5*h_CLE40))*(1+exp(-h_CLE40*(w_CLE40-0.5)))))-(lambda_CLE40*CLE40)
	      
		 list(c(dCK, dARR1, dSHY2, dAUXIAAR, dARFR, dARF10, dARF5, dAUX, dSCR, dSHR, dMIR165, dPHB, dJKD, dMGP, dWOX5, dCLE40))
	})
}

repeticiones=1000
B=matrix(nrow=repeticiones,ncol=16)

for (i in 1:repeticiones){
	state<-c(CK=runif(1,0,1),ARR1=runif(1,0,1),SHY2=runif(1,0,1),AUXIAAR=runif(1,0,1),ARFR=runif(1,0,1),ARF10=runif(1,0,1),ARF5=runif(1,0,1),AUX=runif(1,0,1),SCR=runif(1,0,1),SHR=runif(1,0,1),MIR165=runif(1,0,1),PHB=runif(1,0,1),JKD=runif(1,0,1),MGP=runif(1,0,1),WOX5=runif(1,0,1),CLE40=runif(1,0,1)) 
	#write("CI", file = "OUTPUT.txt",append = TRUE, sep = " ")
	#write(state, file = "OUTPUT.txt",ncolumns = if(is.character(x)) 1 else 11,append = TRUE, sep = " ")


	times <-seq(0,100,by=0.1)          
   	o <- ode(y=state,times=times,func=network,parms=parameters)
   	for (j in 1:16){
   		if(o[1000,j+1]>0.9){B[i,j]<-1}
   		if(o[1000,j+1]<0.1){B[i,j]<-0}
   	}
}
   	

for (x in 1:repeticiones){
	for(i in 1:11){ 
		cohen<-0
		for(j in 1:16){   
			if(ATTRACTORS[i,j]!=B[x,j]){cohen<-cohen+1}	
		}
			if(cohen==0){counter[1,i]=counter[1,i]+1}
	}
}
vari=sum(counter)
	write(vari, file = "OUTPUT.txt",append = TRUE, sep = " ")

}
