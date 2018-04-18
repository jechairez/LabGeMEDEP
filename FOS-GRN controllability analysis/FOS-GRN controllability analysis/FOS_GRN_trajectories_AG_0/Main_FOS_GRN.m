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
n=13; % Genes number
m=1;  % Control inputs
s=10; % Number of steps
zz=1;
options = lmset('vars',{'AG','AP1','AP2','AP3','EMF1','FT','FUL','LFY','PI','SEP','TFL1','UFO','WUS'}); % Defining the vector of logical variables
A=regexp(fileread('FOS_GRN.txt'),'\n','split');% Open archive for reading
%% Attractors and transition matrix of non-controlled network
[expr,vars]=GetAttractors(A,options,k,ME,MI,MD,MN,MR,MC,MU,MX); % Get attractors from uploaded txt
L=ctimes(expr{:}); % Get transition matrix from uploaded matrix from uploaded txt
[Att_Land,Att_Land_size]=ShowAttractors(L,k); % Get and show attractors //Check ShowAttractors function
%% Transition matrix of controlled-network 
%This matrix L changes every loop according to the input controlled
options = lmset('vars',{'u1','AG','AP1','AP2','AP3','EMF1','FT','FUL','LFY','PI','SEP','TFL1','UFO','WUS'}); % Defining the vector logical variables
for i=1:length(A) % Loop 1
    A{i}=strcat(A{i},[blanks(1) 'u1']); % Add control input on first gene , and then the next one and so on
    A{i}=sprintf('%s',['MD' blanks(1) A{i}]); % Add logic operator on first gene, and then the next one and so on // can be OR, AND, etc operators
    [expr,vars]=GetAttractors(A,options,k,ME,MI,MD,MN,MR,MC,MU,MX); % Get attractors from modified txt (controlled network)
    L = ctimes(expr{:}); % Get transition matrix from modified txt (controlled network)
    %% Compute Ltilde
    % Needed to convert the system from x(t+1)=Lu(t)x(t) into x(t+1)=Ltx(t)u(t) in order to find u(t) sequence  
    % Mc=Mcontrol(L,n,m); % Get the controllability matrix of the system (calculate only for systems where n<5)
    Lt=L*lwij(2^n,2^m); % Calculate Lt
    %% Get reachability set
    % Depending on initial condition (x0) and final condition (xd) taking from attractors landcape , calculate gene per gene if it possible to reach that finaal state 
    for i1=1:Att_Land_size %Loop 2
        x0=lm(((2^n)-bin2dec(num2str(Att_Land(i1,:)))),2^n);% Initial condition (x0)
        Col_check=(Lt^s)*x0; % According to Theorem 9.3 (Deizhang Cheng et al.2011), the system is reachable from x0 to xd if xd is 
        for i2=1:Att_Land_size %Loop 3
            xd=lm(((2^n)-bin2dec(num2str(Att_Land(i2,:)))),2^n); % Final condition (xd)
            jj=[]; % Matrix jj / It contains all the reachable set 
            uu=1; % Dimensional index of matrix jj
            tt=[]; % Matrix tt / It contains all the reachable xd for each x0 
            kk=1; % Dimensional index of matrix tt
            for nn=1:length(Col_check.v) %Loop 4
                if (Col_check.v(nn)~=jj) % Conditional 1
                    jj(uu,:)=Col_check.v(nn);
                    uu=uu+1; %Increment index uu
                end 
                if (Col_check.v(nn)==xd.v) %Conditional 2
                    ff=lm(nn,2^(m*s));
                    tt(kk,:)=dec2bin(ff.n-ff.v,m*s);
                    kk=kk+1; %Increment index kk
                end
            end 
     %% Get Boolean control sequences 
            if (isempty(tt)==1) %Conditional 3
            else       
                res(zz,:)=[zz i i1 i2]; % Sort information
                U{zz,1}=tt; %Matrix U/ It contains all the possible control sequence for every possible transition
                zz=zz+1; %Increment index zz
            end
        end
    end
    A=regexp(fileread('FOS_GRN.txt'),'\n','split');
end
[res2]=SortInfo(res);
[Trajectories_a,Trajectories_b]=ShowTrajectories(res2,Att_Land_size);
