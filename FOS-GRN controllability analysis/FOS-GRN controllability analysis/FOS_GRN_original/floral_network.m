
% ALGEBRAIC TRANSFORMATION OF THE BOOLEAN NETWORK DESCRIBING FLORAL
% MORPHOGENESIS
%
% JD & JCMG - Berkeley, CA - April 11, 2013
%
% THE TRANSFOREMATION USES THE STP TOOLBOX DEVELOPPED BY Daiz Cheng (http://lsc.amss.ac.cn/~dcheng/STP/stp.zip)

% CLEARING ALL THE VARIABLES AND ERASING THE SCREEN

clear all
clc

% DEFINITION OF THE NETWORK
%
% From: 
%                   General Theory of Genotype to Phenotype Mapping: Derivation of Epigenetic Landscapes
%                   from N-Node Complex Gene Regulatory Networks, Carlos Villarreal, Pablo Padilla-Longoria, 
%                   and Elena R. Alvarez-Buylla. PRL 109, 118102
%                   (2012), - Supplementary Material -,  DOI: 10.1103/PhysRevLett.109.118102 
%

logic_Ag = '(!EMF1 & !TFL1 & !AP2)|(!EMF1 & LFY & !AP1)|(!EMF1 & !AP2 & LFY)|(!EMF1 & !TFL1 & LFY & (AG & SEP))|(!EMF1 & (LFY & WUS))';
logic_Ap1 = '(!AG & !TFL1)|(FT & LFY & !AG)|(FT & !AG & !PI)|(LFY & !AG & !PI)|(FT & !AG & !AP3)|(LFY & !AG & !AP3)';
logic_Ap2 = '!TFL1';
logic_Ap3 = '(LFY & UFO)|(PI & SEP & AP3 & (AG|AP1))';
logic_Emf1 = '!LFY';
logic_Ft = '!EMF1';
logic_Ful = '!AP1 & !TFL1';
logic_Lfy = '(!EMF1|!TFL1) | U';
logic_Pi = '(LFY & (AG|AP3))|(PI & SEP & AP3 & (AG|AP1))';
logic_Sep = 'LFY';
logic_Tfl1 = '!AP1 & (EMF1 & !LFY)';
logic_Ufo = 'UFO';
logic_Wus = 'WUS & (!AG|!SEP)' ;


% MATRIX VERSION OF THE LOGICAL RULES FOR EACH NODE IN THE BOOLEAN NETWORK:

matrix_expr_Ag = lmparser(logic_Ag);
matrix_expr_Ap1 = lmparser(logic_Ap1);
matrix_expr_Ap2 = lmparser(logic_Ap2);
matrix_expr_Ap3 = lmparser(logic_Ap3);
matrix_expr_Emf1 = lmparser(logic_Emf1);
matrix_expr_Ft = lmparser(logic_Ft);
matrix_expr_Ful = lmparser(logic_Ful);
matrix_expr_Lfy = lmparser(logic_Lfy);
matrix_expr_Pi = lmparser(logic_Pi);
matrix_expr_Sep = lmparser(logic_Sep);
matrix_expr_Tfl1 = lmparser(logic_Tfl1);
matrix_expr_Ufo = lmparser(logic_Ufo);
matrix_expr_Wus = lmparser(logic_Wus);


% INITIALIZING:

k = 2; % Fixing logical variables as tWO values ones.
ME = lme(k);
MI = lmi(k);
MD = lmd(k);
MN = lmn(k);
MR = lmr(k);
MC = lmc(k);
MU = lmu(k); % dummy matrix
MX = lm([2 1 1 2], 2); % xor


% DEFINING THE VECTOR LOGICAL VARIABLES:

options = lmset('vars',{'AG','AP1','AP2','AP3','EMF1','FT','FUL','LFY','PI','SEP','TFL1','UFO','WUS'});

% BOOLEAN NETWORK CODED 

% write the equations in a cell array as follows:

eqn = {'MD MD MD MD MC MC MN EMF1 MN TFL1 MN AP2 MC MC MN EMF1 LFY MN AP1 MC MC MN EMF1 MN AP2 LFY MC MC MC MN EMF1 MN TFL1 LFY MC AG SEP MC MN EMF1 MC LFY WUS',...
    'MD MD MD MD MD MC MN AG MN TFL1 MC MC FT LFY MN AG MC MC FT MN AG MN PI MC MC LFY MN AG MN PI MC MC FT MN AG MN AP3 MC MC LFY MN AG MN AP3',...
    'MN TFL1',...
    'MD MC LFY UFO MC MC MC PI SEP AP3 MD AG AP1',...
    'MN LFY',...
    'MN EMF1',...
    'MC MN AP1 MN TFL1',...
    'MD MN EMF1 MN TFL1',...
    'MD MC LFY MD AG AP3 MC MC MC PI SEP AP3 MD AG AP1',...
    'LFY',...
    'MC MN AP1 MC EMF1 MN LFY',...
    'UFO',...
    'MC WUS MD MN AG MN SEP'};

% CONVERT THE EQUATION TO A CANONICAL FORM:

eqn = completeeqn(eqn,options); % make each equation have all variables
[expr,vars] = stdform(eqn,options,k);  % get the canonical form for each variable
userview=memory
for i=1:length(expr)            % calculate the structure matrix for each variable
   expr{i} = eval(expr{i});
end

% ALGEBRAIC VERSION OF THE BOOLEAN NETWORK:

L = ctimes(expr{:});

% COMPUTING THE ATRACTOR'S LANDSCAPE OF THE BOOLEAN NETWORK:

[n,l,c,r0,T] = bn(L,k);

fprintf('Number of attractors: %d\n\n',n);
fprintf('Lengths of attractors:\n');
disp(l);
for i=1:length(c)
    fprintf('\nAll attractors are displayed as follows:\n\n');
    fprintf('No. %d (length %d)\n\n',i,l(i));
    disp(c{i});
end