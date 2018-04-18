%This code computes trajectory tracking and control design
%% Initialization
clear all; clc;
k = 2; % Fixing logical variables as two values ones.
ME = lme(k); % equivalence matrix
MI = lmi(k); % inequivalence matrix 
MD = lmd(k); % OR matrix
MN = lmn(k); % NOT matrix
MR = lmr(k); % reduction power matrix
MC = lmc(k); % AND matrix
MU = lmu(k); % dummy matrix
MX = lm([2 1 1 2], 2); % xor
n=2; % Genes number
m=2;  % Control inputs
%s=10; % Number of steps
zz=1; %Index of matrix res
x0='1 0'; %Initial condition in binary
xd='1 1'; %Final condition in binary
s_max=10; %Maximun number of steps
options = lmset('vars',{'u1','u2','x1','x2'}); % Defining the vector of logical variables
A=regexp(fileread('Example16_3.txt'),'\n','split');% Open archive for reading
%% Attractors and transition matrix of non-controlled network
[expr,vars]=GetAttractors(A,options,k,ME,MI,MD,MN,MR,MC,MU,MX); % Get attractors from uploaded txt
L=ctimes(expr{:}); % Get transition matrix from uploaded matrix from uploaded txt
%% Calculating controllability matrix
Mc=Mcontrol(L,n,m);
%% Finding trajectory
x0=lm((2^n)-bin2dec(x0),2^n)
xd=lm((2^n)-bin2dec(xd),2^n)
find1=0;
find2=0;
if Mc(xd.v,x0.v)>0
   for s=1:s_max
       if find1==0
           gpower=Gpower(L,m,s);
           for alpha=0:2^m-1
               if alpha>0
                   if gpower(xd.v,x0.v+alpha*4)>0
                       s
                       alpha+1
                       break;
                   end          
               else
                   if gpower(xd.v,x0.v)>0
                       s
                       alpha+1
                       break;
                   end
               end
           end
       else
       end 
   end
else
    display('No possible')
end