#################################################################################################################################
# This file contains R functions intented to be useful while performing Epigenetic Landscape analysis over BoolNet networks.
# The functions have been tested over limited examples so any correction or suggestion is welcome. 
# The functions are intented to be as generic as possible with the purpouse that limited manual intervention is required and reproducibility might be granted.
#
# Developed by: Juan Antonio Arias Del Angel (PhD student) - Universidad Nacional Autónoma de Mexico - Instituto de Ecología
# Contact: jariasdelangel@gmail.com
# Advisor: Mariana Benitez Keinrad (PhD) - Universidad Nacional Autónoma de Mexico - Instituto de Ecología
# Contact: marianabk@gmail.com 
#
# This script was develop as part of the requirenments for the graduate course "Biologia de sistemas, complejidad, desarrollo biologico y medicina: tratamientos matematicos y computacionales" (August - December 2016) by Dr. Juan Carlos Martinez Garcia and Dr. Elisa Dominguez Huttinger in the PhD Program for Biomedical Sciences. Acknowledgement are given to them for discussion and teaching regarding the methods employed here. Also, acknowledgements to Elena Alvarez-Buylla for introductory seminars to the methods and Carlos Villarreal for explain deeper details regarding to stochastic analysis and fuzzy logic. Also, to Monica Garcia, Mariana Esther Martinez and Joel Herrera for sharing code and useful discussions during the development of this script.  
# 
# References:
# Alvarez-Buylla et al (2008). Floral morphogenesis: stochastic explorations of a gene network epigenetic landscape. Plos one, 3(11), e3626.
# Sanchez-Corrales et al (2010). The Arabidopsis thaliana flower organ specification gene regulatory network determines a robust differentiation process. Journal of theoretical biology, 264(3), 971-983.
# Villarreal et al (2012). General theory of genotype to phenotype mapping: derivation of epigenetic landscapes from N-node complex gene regulatory networks. Physical review letters, 109(11), 118102.
# Davila-Velderrain (2015). Reshaping the epigenetic landscape during early flower development: induction of attractor transitions by relative differences in gene decay rates. BMC systems biology, 9(1), 1.
#################################################################################################################################

# Translate from Boolean Network to Continuous Network.
# Description: The function translate the Boolean Network to a Determinist Continuous Network composed by two principal parts. The first one containes the translation from classic to fuzzy logic using one of two methods. The second one contains an ODE-based description of the rate of changes for each of the nodes. 
# IMPORTANT: This function call the perl script "ConvertTocontinuos.pl". For now, such script must be located in the same directory that the R working directory. 
# Arguments:
# (1) net a file containing a BoolNet network.
# (2) logic is a character string specificating either if Zadeh or Probabilistic logic is employed during the translation. Default: Zadeh.
# (3) eq is a character string specificating either if the equation proposed by Sanchez-Corrales et al., 2010 or Villarreal et al., 2012 is employed by describe the rate of change for the nodes in the network. 

BoolNet2ToContinuous = function(net, logic = "Zadeh", eq = "SQUAD", sep = ","){
    if(logic != "Zadeh" && logic != "Probabilistic"){
        stop(paste("Invalid logic value: logic must be either 'Fuzzy' or 'Probabilistic' "))
    }
    if(eq != "SQUAD" && eq != "Villarreal"){
        stop(paste("Invalid eq value: eq must be either 'LuisMendoza' or 'CarlosVillarreal' "))
    }

    cmd = paste("perl", "ConvertToContinuos.pl", net, logic, eq, sep, "> cont.r")
    system(cmd)
    source("cont.r")
    system("rm cont.r")
    return(network)
    
}

