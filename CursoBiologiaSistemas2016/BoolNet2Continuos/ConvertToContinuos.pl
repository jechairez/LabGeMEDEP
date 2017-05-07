#!/usr/bin/perl;

#################################################################################################################################
# Developed by: Juan Antonio Arias Del Angel 
# Last version: 1.0 (November 26th 2016)
# Contact: jariasdelangel@gmail.com
#
#
# Description: This program is created with the purpouse to automatize the continuous mapping of boolean network over which further analysis might be performed using R. 
#
# Input: The program takes 3 arguments as input. 
#
# (1) The first argument is a plain text file containing a boolean network encoded according to R-library BoolNet specifications. 
# (2) A string specificating the type of fuzzy logic which will be employed to generate a continuous model. 'Fuzzy' specifies Zadeh's logic based in max and min functions and 'Probabilistic' is based in sums and multiplications. 
# (3) A string specificating the type of equation employed to describe the rate of change. "CarlosVillarreal" and "LuisMendoza" specifies for the equations proposed by Villarreal et al. 2012 and Ortíz-Gutierrez et al., 2015, respecitvely. 
################################################################################################################################

use strict;

open(net, $ARGV[0]);

my @BooleanNetwork;
my $BoolFunc;

######################################################
# In this variables it is stored the information regarding to the logic and equations employed to generate the continuous model.

my $logic = $ARGV[1];
my $equation = $ARGV[2];
my $mode = $ARGV[3];
######################################################

# Get the information about the Boolean network stored in the input file. 
while(<net>){
    
    if(!/targets/ && !/\#/){
        $BoolFunc = $_;
        chomp $BoolFunc;
        
        $BoolFunc .= '\\';
        
        $BoolFunc =~ s/ //gi;
        
        push(@BooleanNetwork, $BoolFunc);
    }
    
    
} close(net);

######################################################
my @BoolFunc;
my $LuisMendoza = "dX = ((-exp(0.5*h)+exp(-h*(w_X)))/((1-exp(0.5*h))*(1+exp(-h*(w_X-0.5)))))-(alphaX*X)";
my $CarlosVillarreal = "dX = 1/(1+(exp(-2*b*(w_X-h)))) - (alphaX*X)";
my $eq = ();
my @nodes = ();
######################################################
# Print the beginning of the R function to simulate the continuous model as employed in R-library deSolve.

if($mode eq 'Deterministic'){
    print("network = function(t, state, parameters) {\n");
    print("\twith(as.list(c(state, parameters)),{\n");
}


############################################################################################################
############################################################################################################
print("\t\t\# Input Nodes\n");

# According to the parameter chosen. The inputs nodes are specified according to fuzzy or probabilistic logic.


if($logic eq 'Fuzzy'){


    for($BoolFunc = 0; $BoolFunc < scalar(@BooleanNetwork); $BoolFunc++){
        @BoolFunc = split(',', $BooleanNetwork[$BoolFunc]);
        my $characters = length($BoolFunc[1]);

        my $FuzzyExpression = &getFuzzyContent($BoolFunc[1], 0);
        $FuzzyExpression =~ s/!/1-/gi;

        print("\t\tw\_$BoolFunc[0] = $FuzzyExpression\n");

    }
} 

if($logic eq 'Probabilistic'){

    for($BoolFunc = 0; $BoolFunc < scalar(@BooleanNetwork); $BoolFunc++){
        @BoolFunc = split(',', $BooleanNetwork[$BoolFunc]);
        my $characters = length($BoolFunc[1]);

        my $FuzzyExpression = &getProbabilisticContent($BoolFunc[1], 0);
        $FuzzyExpression =~ s/!/1-/gi;

        print("\t\tw\_$BoolFunc[0] = $FuzzyExpression\n");

    }
}
############################################################################################################
############################################################################################################
print "\n\n";

if($mode eq 'Deterministic'){
    # According to the parameter chosen. The rates of changes are specified according to Villarreal or Mendoza. 
    print("\t\t\# Rates of Change\n");
    if($equation eq 'LuisMendoza'){
        for($BoolFunc = 0; $BoolFunc < scalar(@BooleanNetwork); $BoolFunc++){
            $eq = $LuisMendoza;
            @BoolFunc = split(',', $BooleanNetwork[$BoolFunc]);
            $eq =~ s/X/$BoolFunc[0]/g;
            print "\t\t$eq\n";

            push(@nodes, 'd' . $BoolFunc[0]);
        }
    }

    if($equation eq 'CarlosVillarreal'){
        for($BoolFunc = 0; $BoolFunc < scalar(@BooleanNetwork); $BoolFunc++){
            $eq = $CarlosVillarreal;
            @BoolFunc = split(',', $BooleanNetwork[$BoolFunc]);
            $eq =~ s/X/$BoolFunc[0]/g;
            print "\t\t$eq\n";

            push(@nodes, 'd' . $BoolFunc[0]);
        }
    }
}

if($mode eq 'Stochastic'){
    print "\n\n\t\t\# Add noise\n";
    foreach my $node (@nodes){
        print "\t\t$node = $node + intensity*sqrt(dt)*rnorm(1)\n";
        print "\t\tif($node < 0) $node = 0\n";
        print "\t\tif($node > 1) $node = 1\n";
    }
}

############################################################################################################
############################################################################################################

# Print the last of the R function to simulate the continuous model as employed in R-library deSolve

my $nodes = join(',', @nodes);

print("\n\n");
print("\t\tlist(c($nodes))\n");
print("\t})\n}\n");



#################################################################################################################################
#######################                                                                                   #######################
#######################          Subroutines employed to transform the boolean into fuzzy network         #######################
#######################                                                                                   #######################
#################################################################################################################################

sub getFuzzyContent{

    my $string = $_[0];
    my $pos = $_[1];
    my @variables = ();
    my @characters = ();
    my $position = 0;
    my $fuzzy = ();
    my $nextfuzzy = ();
    my $i = ();
    my $c = ();

    
    @characters = split('', $BoolFunc[1]);
    for($c = $pos; $c < scalar(@characters); $c++){
        
        
######################################################
######################################################
# Actions to be taken when a closing bracket is found. 

        if($characters[$c] eq ')'){
            if(scalar(@variables) == 1){
                $fuzzy = $variables[0];
            }
            elsif(scalar(@variables) == 3){
                if($variables[1] eq '&'){
                    $fuzzy = "min(" . $variables[0] . "," . $variables[2] . ")";
                } 
                
                if($variables[1] eq '|'){
                    $fuzzy = "max(" . $variables[0] . "," . $variables[2] . ")";
                }
            }
            elsif(scalar(@variables) >= 5){
                if($variables[1] eq '&'){
                    $fuzzy = "min(" . $variables[0] . "," . $variables[2] . ")";
                } 
                
                if($variables[1] eq '|'){
                    $fuzzy = "max(" . $variables[0] . "," . $variables[2] . ")";
                }
                
                
                for($i = 3; $i < scalar(@variables); $i += 2){
                    if($variables[$i] eq '&'){
                        $nextfuzzy = "min(" . $fuzzy . "," . $variables[$i + 1] . ")";
                    } 
                
                    if($variables[$i] eq '|'){
                        $nextfuzzy = "max(" . $fuzzy . "," . $variables[$i + 1] . ")";
                    }
                    
                    $fuzzy = $nextfuzzy;
                }
                
               
                
            }
            my $result = $fuzzy . "\t" . $c;
            return($result);
        } 
        
######################################################
######################################################
# Actions to be taken when an openning bracket is found. 
# This is the recursive part of the function.


        elsif($characters[$c] eq '('){
            my $data = &getFuzzyContent($string, $c + 1);
            my @data = split("\t", $data);
            $variables[$position] .= $data[0];
            $c = $data[1]; 
        } 
        
######################################################
######################################################
# Actions to be taken when binary logic operators are found. 

        elsif($characters[$c] eq '&' || $characters[$c] eq '|'){
            $position++;
            $variables[$position] .= $characters[$c];
            $position++;
        } 
        
######################################################
######################################################
# Actions to be taken when the end of the line is reaached. 

        elsif ($characters[$c] eq '\\'){
           # $variables[$position] .= $characters[$c];
            
            if(scalar(@variables) == 1){
                $fuzzy = $variables[0];
            }
            elsif(scalar(@variables) == 3){
                if($variables[1] eq '&'){
                    $fuzzy = "min(" . $variables[0] . "," . $variables[2] . ")";
                } 
                
                if($variables[1] eq '|'){
                    $fuzzy = "max(" . $variables[0] . "," . $variables[2] . ")";
                }
            }
            elsif(scalar(@variables) >= 5){
                if($variables[1] eq '&'){
                    $fuzzy = "min(" . $variables[0] . "," . $variables[2] . ")";
                } 
                
                if($variables[1] eq '|'){
                    $fuzzy = "max(" . $variables[0] . "," . $variables[2] . ")";
                }
                
                for($i = 3; $i < scalar(@variables); $i += 2){
                    if($variables[$i] eq '&'){
                        $nextfuzzy = "min(" . $fuzzy . "," . $variables[$i + 1] . ")";
                    } 
                
                    if($variables[$i] eq '|'){
                        $nextfuzzy = "max(" . $fuzzy . "," . $variables[$i + 1] . ")";
                    }
                    
                    $fuzzy = $nextfuzzy;
                }
                
               
                
            }
            return($fuzzy);
            
        }
        
######################################################
######################################################
# Actions to be taken when other characters are found.

        else {
            $variables[$position] .= $characters[$c];
        }
    }
    
    
    
}

############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################

sub getProbabilisticContent{

    my $string = $_[0];
    my $pos = $_[1];
    my @variables = ();
    my @characters = ();
    my $position = 0;
    my $fuzzy = ();
    my $nextfuzzy = ();
    my $i = ();
    my $c = ();
    my $term = ();

    
    @characters = split('', $BoolFunc[1]);
    for($c = $pos; $c < scalar(@characters); $c++){
        
        
######################################################
######################################################
# Actions to be taken when a closing bracket is found.
# An important difference with the other method is that parenthesis are set between variables. 

        if($characters[$c] eq ')'){
        
            for($term = 0; $term < scalar(@variables); $term++){
                $variables[$term] =~ s/!/1-/gi;
            }
            
            if(scalar(@variables) == 1){
                $fuzzy = "(" . $variables[0] . ")";
            }
            elsif(scalar(@variables) == 3){
                if($variables[1] eq '&'){
                    $fuzzy = "(" . $variables[0] . ")" . "*" . "(" . $variables[2] . ")";
                } 
                
                if($variables[1] eq '|'){
                    $fuzzy = "(" . $variables[0] . ")" . "+" . "(" . $variables[2] . ")" . "-"  . "(" . $variables[0] . ")" . "*" . "(" . $variables[2] . ")";
                }
            }
            elsif(scalar(@variables) >= 5){
                if($variables[1] eq '&'){
                    $fuzzy = "(" . $variables[0] . ")" . "*" . "(" . $variables[2] . ")";
                } 
                
                if($variables[1] eq '|'){
                    $fuzzy = "(" . $variables[0] . ")" . "+" . "(" . $variables[2] . ")" . "-"  . "(" . $variables[0] . ")" . "*" . "(" . $variables[2] . ")";
                }
                
                
                for($i = 3; $i < scalar(@variables); $i += 2){
                    if($variables[$i] eq '&'){
                        $nextfuzzy = "(" . $fuzzy . ")" . "*" . "(" .  $variables[$i + 1] . ")";
                    } 
                
                    if($variables[$i] eq '|'){
                        $nextfuzzy = "(" . $fuzzy . ")" . "+" . "(" . $variables[$i + 1] . ")" . "-" . "(" . $fuzzy . ")" . "*" . "(" . $variables[$i + 1] . ")";
                    }
                    
                    $fuzzy = $nextfuzzy;
                }
                
               
                
            }
            my $result = "(" . $fuzzy . ")" . "\t" . $c;
            return($result);
        } 
        
######################################################
######################################################
# Actions to be taken when an openning bracket is found.
# This is the recursive part of the function.

        elsif($characters[$c] eq '('){
            my $data = &getProbabilisticContent($string, $c + 1);
            my @data = split("\t", $data);
            $variables[$position] .= $data[0];
            $c = $data[1]; 
        } 
        
######################################################
######################################################
# Actions to be taken when binary logic operators are found.

        elsif($characters[$c] eq '&' || $characters[$c] eq '|'){
            $position++;
            $variables[$position] .= $characters[$c];
            $position++;
        } 
        
######################################################
######################################################
# Actions to be taken when the end of the line is reaached. 
# An important difference with the other method is that parenthesis are set between variables. 


        elsif ($characters[$c] eq '\\'){
           # $variables[$position] .= $characters[$c];
            
            for($term = 0; $term < scalar(@variables); $term++){
                $variables[$term] =~ s/!/1-/gi;
            }
            
            if(scalar(@variables) == 1){
                $fuzzy = "(" . $variables[0] . ")";
            }
            elsif(scalar(@variables) == 3){
                if($variables[1] eq '&'){
                    $fuzzy = "(" . $variables[0] . ")" . "*" . "(" . $variables[2] . ")";
                } 
                
                if($variables[1] eq '|'){
                    $fuzzy = "(" . $variables[0] . ")" . "+" . "(" . $variables[2] . ")" . "-"  . "(" . $variables[0] . ")" . "*" . "(" . $variables[2] . ")";
                }
            }
            elsif(scalar(@variables) >= 5){
                if($variables[1] eq '&'){
                    $fuzzy = "(" . $variables[0] . ")" . "*" . "(" . $variables[2] . ")";
                } 
                
                if($variables[1] eq '|'){
                    $fuzzy = "(" . $variables[0] . ")" . "+" . "(" . $variables[2] . ")" . "-"  . "(" . $variables[0] . ")" . "*" . "(" . $variables[2] . ")";
                }
                
                
                for($i = 3; $i < scalar(@variables); $i += 2){
                    if($variables[$i] eq '&'){
                        $nextfuzzy = "(" . $fuzzy . ")" . "*" . "(" .  $variables[$i + 1] . ")";
                    } 
                
                    if($variables[$i] eq '|'){
                        $nextfuzzy = "(" . $fuzzy . ")" . "+" . "(" . $variables[$i + 1] . ")" . "-" . "(" . $fuzzy . ")" . "*" . "(" . $variables[$i + 1] . ")";
                    }
                    
                    $fuzzy = $nextfuzzy;
                }
                
               
                
            }
            return($fuzzy);
            
        }
######################################################
######################################################
######################################################
######################################################
# Actions to be taken when other characters are found.

        else {
            $variables[$position] .= $characters[$c];
        }
    }
    
    
    
}


#################################################################################################################################
#################################################################################################################################
