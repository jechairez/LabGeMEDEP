#Este codigo sirve para ver qué pasa con los atractores del sistema booleano en la versión continua del modelo 
#Los parámetros h (50) y lambda (1) son los mismos para todos los nodos
#El atractor a analizar se guarda en la variable 'state' , y se puede explorar cada uno descomentándolo en las sig líneas
#Los resultados son que los 11 atractores fijos son estables también en el sistema continuo. 
#Los atractores cíclicos no lo son, y convergen a alguno de los fijos. 
#García-Gómez, Mónica L., Eugenio Azpeitia, and Elena R. Álvarez-Buylla. "A dynamic genetic-hormonal regulatory network model explains multiple cellular behaviors of the root apical meristem of Arabidopsis thalian." PLOS Computational Biology 13.4 (2017): e1005488.

library(deSolve)
parameters <- c(h = 50,lambda = 1)
#INITIAL CONDITIONS
#state<-c(CK=0,ARR1=0,SHY2=0,AUXIAAR=0,ARFR=1,ARF10=0,ARF5=1,AUX=1,SCR=1,SHR=1,MIR165=1,PHB=0,JKD=1,MGP=0,WOX5=1,CLE40=0) #1 QC
#state<-c(CK=0,ARR1=0,SHY2=0,AUXIAAR=0,ARFR=1,ARF10=0,ARF5=0,AUX=1,SCR=1,SHR=1,MIR165=1,PHB=0,JKD=1,MGP=1,WOX5=0,CLE40=0) #2 END-PD
#state<-c(CK=0,ARR1=0,SHY2=0,AUXIAAR=0,ARFR=1,ARF10=1,ARF5=1,AUX=1,SCR=0,SHR=1,MIR165=1,PHB=0,JKD=0,MGP=0,WOX5=0,CLE40=0) #3 P.P.PD
#state<-c(CK=0,ARR1=0,SHY2=0,AUXIAAR=0,ARFR=1,ARF10=1,ARF5=1,AUX=1,SCR=0,SHR=1,MIR165=0,PHB=1,JKD=0,MGP=0,WOX5=0,CLE40=0) #4 C.P.PD
#state<-c(CK=0,ARR1=0,SHY2=0,AUXIAAR=1,ARFR=0,ARF10=0,ARF5=0,AUX=0,SCR=1,SHR=1,MIR165=1,PHB=0,JKD=1,MGP=1,WOX5=0,CLE40=0) #5 END-TD
#state<-c(CK=0,ARR1=0,SHY2=0,AUXIAAR=1,ARFR=0,ARF10=0,ARF5=0,AUX=0,SCR=0,SHR=1,MIR165=1,PHB=0,JKD=0,MGP=0,WOX5=0,CLE40=0) #6 P.P.TD
#state<-c(CK=1,ARR1=1,SHY2=1,AUXIAAR=1,ARFR=0,ARF10=0,ARF5=0,AUX=0,SCR=0,SHR=1,MIR165=0,PHB=1,JKD=0,MGP=0,WOX5=0,CLE40=0) #7 C.P.TD1
#state<-c(CK=1,ARR1=1,SHY2=1,AUXIAAR=1,ARFR=0,ARF10=0,ARF5=0,AUX=0,SCR=0,SHR=0,MIR165=0,PHB=1,JKD=0,MGP=0,WOX5=0,CLE40=1) #8 C.P.TD2
#state<-c(CK=1,ARR1=1,SHY2=0,AUXIAAR=0,ARFR=1,ARF10=1,ARF5=1,AUX=1,SCR=0,SHR=0,MIR165=0,PHB=1,JKD=0,MGP=0,WOX5=0,CLE40=1) #9 C.P.TD3
#state<-c(CK=1,ARR1=1,SHY2=1,AUXIAAR=1,ARFR=0,ARF10=0,ARF5=0,AUX=0,SCR=0,SHR=0,MIR165=1,PHB=0,JKD=0,MGP=0,WOX5=0,CLE40=1) #10 RC-1
#state<-c(CK=1,ARR1=1,SHY2=0,AUXIAAR=0,ARFR=1,ARF10=1,ARF5=1,AUX=1,SCR=0,SHR=0,MIR165=1,PHB=0,JKD=0,MGP=0,WOX5=0,CLE40=1) #11 RC-2
#state<-c(CK=1,ARR1=1,SHY2=1,AUXIAAR=1,ARFR=0,ARF10=0,ARF5=0,AUX=0,SCR=0,SHR=0,MIR165=0,PHB=0,JKD=0,MGP=0,WOX5=0,CLE40=1) #CYC1.1
#state<-c(CK=1,ARR1=1,SHY2=1,AUXIAAR=1,ARFR=0,ARF10=0,ARF5=0,AUX=0,SCR=0,SHR=0,MIR165=1,PHB=1,JKD=0,MGP=0,WOX5=0,CLE40=1) #CYC1.2
#state<-c(CK=1,ARR1=1,SHY2=0,AUXIAAR=0,ARFR=1,ARF10=1,ARF5=1,AUX=1,SCR=0,SHR=0,MIR165=0,PHB=0,JKD=0,MGP=0,WOX5=0,CLE40=1) #CYC2.1
#state<-c(CK=1,ARR1=1,SHY2=0,AUXIAAR=0,ARFR=1,ARF10=1,ARF5=1,AUX=1,SCR=0,SHR=0,MIR165=1,PHB=1,JKD=0,MGP=0,WOX5=0,CLE40=1) #CYC2.2
#state<-c(CK=1,ARR1=0,SHY2=1,AUXIAAR=1,ARFR=0,ARF10=0,ARF5=0,AUX=0,SCR=0,SHR=1,MIR165=0,PHB=0,JKD=0,MGP=0,WOX5=0,CLE40=0) #CYC3.1
#state<-c(CK=0,ARR1=1,SHY2=0,AUXIAAR=1,ARFR=0,ARF10=0,ARF5=0,AUX=0,SCR=0,SHR=1,MIR165=1,PHB=1,JKD=0,MGP=0,WOX5=0,CLE40=0) #CYC3.2
#state<-c(CK=0,ARR1=0,SHY2=0,AUXIAAR=0,ARFR=1,ARF10=1,ARF5=1,AUX=1,SCR=0,SHR=1,MIR165=0,PHB=0,JKD=0,MGP=0,WOX5=0,CLE40=0) #CYC4.1
#state<-c(CK=0,ARR1=0,SHY2=0,AUXIAAR=0,ARFR=1,ARF10=1,ARF5=1,AUX=1,SCR=0,SHR=1,MIR165=1,PHB=1,JKD=0,MGP=0,WOX5=0,CLE40=0) #CYC4.2
#state<-c(CK=0,ARR1=0,SHY2=0,AUXIAAR=0,ARFR=1,ARF10=0,ARF5=0,AUX=1,SCR=1,SHR=1,MIR165=1,PHB=0,JKD=1,MGP=0,WOX5=0,CLE40=0) #CYC5.1
#state<-c(CK=0,ARR1=0,SHY2=0,AUXIAAR=0,ARFR=1,ARF10=0,ARF5=1,AUX=1,SCR=1,SHR=1,MIR165=1,PHB=0,JKD=1,MGP=1,WOX5=0,CLE40=0) #CYC5.2
#state<-c(CK=0,ARR1=0,SHY2=0,AUXIAAR=0,ARFR=1,ARF10=0,ARF5=0,AUX=1,SCR=1,SHR=1,MIR165=1,PHB=0,JKD=1,MGP=1,WOX5=1,CLE40=0) #CYC5.3
#state<-c(CK=0,ARR1=0,SHY2=0,AUXIAAR=0,ARFR=1,ARF10=0,ARF5=1,AUX=1,SCR=1,SHR=1,MIR165=1,PHB=0,JKD=1,MGP=0,WOX5=0,CLE40=0) #CYC6.1
#state<-c(CK=0,ARR1=0,SHY2=0,AUXIAAR=0,ARFR=1,ARF10=0,ARF5=1,AUX=1,SCR=1,SHR=1,MIR165=1,PHB=0,JKD=1,MGP=1,WOX5=1,CLE40=0) #CYC6.2
state<-c(CK=0,ARR1=0,SHY2=0,AUXIAAR=0,ARFR=1,ARF10=0,ARF5=0,AUX=1,SCR=1,SHR=1,MIR165=1,PHB=0,JKD=1,MGP=0,WOX5=1,CLE40=0) #CYC6.3
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

times <-seq(0,20,by=0.1)
   o <- ode(y=state,times=times,func=network,parms=parameters)
	colores <- rainbow(length(state),start=0.1,end=1) 
	salida.pdf <- paste("attractor",".pdf",sep="")
	pdf(salida.pdf, height=4)
	ejeY <- 0
		#lines(o[,1],o[,(16+1)], col=colores[16])
	for(m in 1:length(state)){ ejeY <- max(ejeY,o[,(m+1)])} #Calcute the maximum y axis value
	plot(o[,1],o[,2], ylim=c(0,ejeY), type = "l", col=colores[1], xlab="Timesteps", ylab="Relative activity")
        for(n2 in 2:length(state)){ lines(o[,1],o[,(n2+1)], col=colores[n2])}
        legend("topright", lwd=2, names(state), col=colores, cex=0.6, box.lwd = 0, box.col = "white", bg = "white") #legend
dev.off()
o[200,]

