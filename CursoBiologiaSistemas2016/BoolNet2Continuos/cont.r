network = function(t, state, parameters) {
	with(as.list(c(state, parameters)),{
		# Input Nodes
		w_CYCD31 = (1-SCF)
		w_SCF = (1-APCC)*(((((E2Fb)*(((1-RBR)+(((1-KRP1)*(CYCD31)))-(1-RBR)*(((1-KRP1)*(CYCD31)))))))+(MYB3R14)-(((E2Fb)*(((1-RBR)+(((1-KRP1)*(CYCD31)))-(1-RBR)*(((1-KRP1)*(CYCD31)))))))*(MYB3R14)))
		w_RBR = (((KRP1)+(1-CYCD31)-(KRP1)*(1-CYCD31)))*(((((E2Fa)*(1-RBR)))+(MYB3R14)-(((E2Fa)*(1-RBR)))*(MYB3R14)))
		w_E2Fa = (((E2Fa)+(1-E2Fc)-(E2Fa)*(1-E2Fc)))*(1-((CDKB11)*(CYCA23)))
		w_E2Fb = (((E2Fa)*(1-RBR)))
		w_E2Fc = (1-(((SCF)*(1-KRP1))*(CYCD31)))*(((((E2Fa)*(1-RBR)))+(MYB3R14)-(((E2Fa)*(1-RBR)))*(MYB3R14)))
		w_E2Fe = (((1-E2Fc)+(((E2Fb)*(((1-RBR)+(((1-KRP1)*(CYCD31)))-(1-RBR)*(((1-KRP1)*(CYCD31)))))))-(1-E2Fc)*(((E2Fb)*(((1-RBR)+(((1-KRP1)*(CYCD31)))-(1-RBR)*(((1-KRP1)*(CYCD31)))))))))+(MYB77)-(((1-E2Fc)+(((E2Fb)*(((1-RBR)+(((1-KRP1)*(CYCD31)))-(1-RBR)*(((1-KRP1)*(CYCD31)))))))-(1-E2Fc)*(((E2Fb)*(((1-RBR)+(((1-KRP1)*(CYCD31)))-(1-RBR)*(((1-KRP1)*(CYCD31)))))))))*(MYB77)
		w_MYB77 = (E2Fb)*(((1-RBR)+(((1-KRP1)*(CYCD31)))-(1-RBR)*(((1-KRP1)*(CYCD31)))))
		w_MYB3R14 = (MYB77)+((((MYB3R14)*(CYCB11))*(1-KRP1)))-(MYB77)*((((MYB3R14)*(CYCB11))*(1-KRP1)))
		w_CYCB11 = (1-APCC)*((((MYB3R14)+(MYB77)-(MYB3R14)*(MYB77))+((((((1-RBR)+(((1-KRP1)*(CYCD31)))-(1-RBR)*(((1-KRP1)*(CYCD31)))))*(E2Fb))*(1-E2Fc)))-((MYB3R14)+(MYB77)-(MYB3R14)*(MYB77))*((((((1-RBR)+(((1-KRP1)*(CYCD31)))-(1-RBR)*(((1-KRP1)*(CYCD31)))))*(E2Fb))*(1-E2Fc)))))
		w_CDKB11 = (((((((1-RBR)+(((1-KRP1)*(CYCD31)))-(1-RBR)*(((1-KRP1)*(CYCD31)))))*(E2Fb))*(1-E2Fc)))+(MYB3R14)-((((((1-RBR)+(((1-KRP1)*(CYCD31)))-(1-RBR)*(((1-KRP1)*(CYCD31)))))*(E2Fb))*(1-E2Fc)))*(MYB3R14))+(MYB77)-(((((((1-RBR)+(((1-KRP1)*(CYCD31)))-(1-RBR)*(((1-KRP1)*(CYCD31)))))*(E2Fb))*(1-E2Fc)))+(MYB3R14)-((((((1-RBR)+(((1-KRP1)*(CYCD31)))-(1-RBR)*(((1-KRP1)*(CYCD31)))))*(E2Fb))*(1-E2Fc)))*(MYB3R14))*(MYB77)
		w_CYCA23 = (1-APCC)*(((MYB3R14)+(MYB77)-(MYB3R14)*(MYB77)))
		w_KRP1 = (((MYB77)+(MYB3R14)-(MYB77)*(MYB3R14)))*(1-(((CDKB11)*(CYCA23))*(SCF)))
		w_APCC = (1-E2Fe)*((((((E2Fa)*(1-RBR)))+(MYB3R14)-(((E2Fa)*(1-RBR)))*(MYB3R14))+(MYB77)-((((E2Fa)*(1-RBR)))+(MYB3R14)-(((E2Fa)*(1-RBR)))*(MYB3R14))*(MYB77)))


		# Rates of Change
		dCYCD31 = 1/(1+(exp(-2*b*(w_CYCD31-h)))) - (alphaCYCD31*CYCD31)
		dSCF = 1/(1+(exp(-2*b*(w_SCF-h)))) - (alphaSCF*SCF)
		dRBR = 1/(1+(exp(-2*b*(w_RBR-h)))) - (alphaRBR*RBR)
		dE2Fa = 1/(1+(exp(-2*b*(w_E2Fa-h)))) - (alphaE2Fa*E2Fa)
		dE2Fb = 1/(1+(exp(-2*b*(w_E2Fb-h)))) - (alphaE2Fb*E2Fb)
		dE2Fc = 1/(1+(exp(-2*b*(w_E2Fc-h)))) - (alphaE2Fc*E2Fc)
		dE2Fe = 1/(1+(exp(-2*b*(w_E2Fe-h)))) - (alphaE2Fe*E2Fe)
		dMYB77 = 1/(1+(exp(-2*b*(w_MYB77-h)))) - (alphaMYB77*MYB77)
		dMYB3R14 = 1/(1+(exp(-2*b*(w_MYB3R14-h)))) - (alphaMYB3R14*MYB3R14)
		dCYCB11 = 1/(1+(exp(-2*b*(w_CYCB11-h)))) - (alphaCYCB11*CYCB11)
		dCDKB11 = 1/(1+(exp(-2*b*(w_CDKB11-h)))) - (alphaCDKB11*CDKB11)
		dCYCA23 = 1/(1+(exp(-2*b*(w_CYCA23-h)))) - (alphaCYCA23*CYCA23)
		dKRP1 = 1/(1+(exp(-2*b*(w_KRP1-h)))) - (alphaKRP1*KRP1)
		dAPCC = 1/(1+(exp(-2*b*(w_APCC-h)))) - (alphaAPCC*APCC)


		list(c(dCYCD31,dSCF,dRBR,dE2Fa,dE2Fb,dE2Fc,dE2Fe,dMYB77,dMYB3R14,dCYCB11,dCDKB11,dCYCA23,dKRP1,dAPCC))
	})
}
