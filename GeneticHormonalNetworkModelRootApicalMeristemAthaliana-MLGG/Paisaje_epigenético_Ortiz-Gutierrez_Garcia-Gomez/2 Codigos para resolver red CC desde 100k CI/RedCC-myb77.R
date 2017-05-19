#!/usr/var/env Rscript

library(deSolve)

parameters <- c(h = 15,lambda = 1)

for(s in 1:100000){
    state <- c(APC=sample(0:1,1),KRP1=sample(0:1,1),CYCA23=sample(0:1,1),CDKB11=sample(0:1,1),CYCB11=sample(0:1,1),MYB3R=sample(0:1,1),MYB77=sample(0:1,1),E2Fe=sample(0:1,1),E2Fc=sample(0:1,1),E2Fb=sample(0:1,1),E2Fa=sample(0:1,1),RBR=sample(0:1,1),SCF=sample(0:1,1),CYCD31=sample(0:1,1))

network<-function(t, state, parameters) {
	with(as.list(c(state, parameters)),{
	     #node inputs
		w_APC = min(1-E2Fe, max(min(E2Fa, 1-RBR), MYB3R, MYB77))
		w_KRP1 = min(max(MYB77, MYB3R), 1-min(CDKB11, CYCA23, SCF))
		w_CYCA23 = min(1-APC, max(MYB3R, MYB77))
		w_CDKB11 = max(min(max(1-RBR, min(1-KRP1, CYCD31)), E2Fb, 1-E2Fc), MYB3R, MYB77)
		w_CYCB11 = min(1-APC, max(MYB3R, MYB77, min(max(1-RBR, min(1-KRP1, CYCD31)), E2Fb, 1-E2Fc)))
		w_MYB3R = max(MYB77, min(MYB3R, CYCB11, 1-KRP1))
		w_MYB77 = min(E2Fb, max(1-RBR, min(1-KRP1, CYCD31)))
		w_E2Fe = max(1-E2Fc, min(E2Fb, max(1-RBR, min(1-KRP1, CYCD31))), MYB77)
		w_E2Fc = min(1-min(SCF, 1-KRP1, CYCD31), max(min(E2Fa, 1-RBR), MYB3R))
		w_E2Fb = min(E2Fa, 1-RBR)
		w_E2Fa = min(max(E2Fa, 1-E2Fc), 1-min(CDKB11, CYCA23))
		w_RBR = min(max(KRP1 , 1-CYCD31), max(min(E2Fa, 1-RBR), MYB3R))
		w_SCF = min(1-APC, max(min(max(1-RBR, min(1-KRP1, CYCD31)), E2Fb), MYB3R))
		w_CYCD31 = 1-SCF
						
	     #rates of change
		dAPC <- ((-exp(0.5*h)+exp(-h*(w_APC)))/((1-exp(0.5*h))*(1+exp(-h*(w_APC-0.5)))))-(lambda*APC)
		dKRP1 <- ((-exp(0.5*h)+exp(-h*(w_KRP1)))/((1-exp(0.5*h))*(1+exp(-h*(w_KRP1-0.5)))))-(lambda*KRP1)
		dCYCA23 <- ((-exp(0.5*h)+exp(-h*(w_CYCA23)))/((1-exp(0.5*h))*(1+exp(-h*(w_CYCA23-0.5)))))-(lambda*CYCA23)
		dCDKB11 <- ((-exp(0.5*h)+exp(-h*(w_CDKB11)))/((1-exp(0.5*h))*(1+exp(-h*(w_CDKB11-0.5)))))-(lambda*CDKB11)
		dCYCB11 <- ((-exp(0.5*h)+exp(-h*(w_CYCB11)))/((1-exp(0.5*h))*(1+exp(-h*(w_CYCB11-0.5)))))-(lambda*CYCB11)
		dMYB3R <- ((-exp(0.5*h)+exp(-h*(w_MYB3R)))/((1-exp(0.5*h))*(1+exp(-h*(w_MYB3R-0.5)))))-(lambda*MYB3R)
		dMYB77 <- ((-exp(0.5*h)+exp(-h*(w_MYB77)))/((1-exp(0.5*h))*(1+exp(-h*(w_MYB77-0.5)))))-(5*MYB77)
		dE2Fe<- ((-exp(0.5*h)+exp(-h*(w_E2Fe)))/((1-exp(0.5*h))*(1+exp(-h*(w_E2Fe-0.5)))))-(lambda*E2Fe)
		dE2Fc <- ((-exp(0.5*h)+exp(-h*(w_E2Fc)))/((1-exp(0.5*h))*(1+exp(-h*(w_E2Fc-0.5)))))-(lambda*E2Fc)
		dE2Fb <- ((-exp(0.5*h)+exp(-h*(w_E2Fb)))/((1-exp(0.5*h))*(1+exp(-h*(w_E2Fb-0.5)))))-(lambda*E2Fb)
		dE2Fa <- ((-exp(0.5*h)+exp(-h*(w_E2Fa)))/((1-exp(0.5*h))*(1+exp(-h*(w_E2Fa-0.5)))))-(lambda*E2Fa)
		dRBR <- ((-exp(0.5*h)+exp(-h*(w_RBR)))/((1-exp(0.5*h))*(1+exp(-h*(w_RBR-0.5)))))-(lambda*RBR)
		dSCF <- ((-exp(0.5*h)+exp(-h*(w_SCF)))/((1-exp(0.5*h))*(1+exp(-h*(w_SCF-0.5)))))-(lambda*SCF)
		dCYCD31 <- ((-exp(0.5*h)+exp(-h*(w_CYCD31)))/((1-exp(0.5*h))*(1+exp(-h*(w_CYCD31-0.5)))))-(lambda*CYCD31)
		 
		 # return the rate of change
		 list(c(dAPC, dKRP1, dCYCA23, dCDKB11, dCYCB11, dMYB3R, dMYB77, dE2Fe, dE2Fc, dE2Fb, dE2Fa, dRBR, dSCF, dCYCD31))
	})
}

times <-seq(0,20,by=0.01)
        o <- ode(y=state,times=times,func=network,parms=parameters)
        ifelse(s==1,suma<-(sum(o[2001,])-20),suma<-c(suma,(sum(o[2001,])-20))) 
        }

      pdf("histSuma-myb77.pdf")
              boxplot(suma,ylim=c(0,14))
       dev.off()
