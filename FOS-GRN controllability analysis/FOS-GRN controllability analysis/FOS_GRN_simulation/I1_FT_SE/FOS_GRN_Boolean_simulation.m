%This code computes the FOS-GRN under the effect of the controller

clear all
clc
%% Initial conditions
%I1=[0 0 0 0 1 0 0 0 0 0 1 0 0] size of basin=136, De=260
%I2=[0 0 0 0 1 0 0 0 0 0 1 1 0] size of basin=136, De=262
%I3=[0 0 0 0 1 0 0 0 0 0 1 1 1] size of basin=72, De=263
%I4=[0 0 0 0 1 0 0 0 0 0 1 0 1] size of basin=72, De=261
%SE=[0 1 1 0 0 1 0 1 0 1 0 0 0] size of basin=812
%PE1=[0 1 1 1 0 1 0 1 1 1 0 0 0] size of basin=12
%PE2=[0 1 1 1 0 1 0 1 1 1 0 1 0] size of basin=824
%ST1=[1 0 1 1 0 1 1 1 1 1 0 0 0] size of basin=94
%ST2=[1 0 1 1 0 1 1 1 1 1 0 1 0] size of basin=3064
%CAR=[1 0 1 0 0 1 1 1 1 1 0 0 0] size of basin=2970
x0=[0 0 0 0 1 0 0 0 0 0 1 1 0];
AG=x0(:,1);
AP1=x0(:,2);
AP2=x0(:,3);
AP3=x0(:,4);
EMF1=x0(:,5);
FT=x0(:,6);
FUL=x0(:,7);
LFY=x0(:,8);
PI=x0(:,9);
SEP=x0(:,10);
TFL1=x0(:,11);
UFO=x0(:,12);
WUS=x0(:,13);
%% Simulation on sth steps
i=1;
s=6;
AG(:,i)=AG;
AP1(:,i)=AP1;
AP2(:,i)=AP2;
AP3(:,i)=AP3;
EMF1(:,i)=EMF1;
FT(:,i)=FT;
FUL(:,i)=FUL;
LFY(:,i)=LFY;
PI(:,i)=PI;
SEP(:,i)=SEP;
TFL1(:,i)=TFL1;
LFY(:,i)=LFY;
UFO(:,i)=UFO;
WUS(:,i)=WUS;
u1=[1 1 0 0 0 0 0];
% u2=[1 0 1 1 1];
% u3=[1 1 0 1 1];
% u4=[1 1 1 1 1 1 1 1 1 0 1];
x0=[AG(:,1) AP1(:,1) AP2(:,1) AP3(:,1) EMF1(:,1) FT(:,1) FUL(:,1) LFY(:,1) PI(:,1) SEP(:,1) TFL1(:,1) UFO(:,1) WUS(:,1)]


for i=1:s
    AG(:,i+1)=(not(EMF1(:,i)) & not(TFL1(:,i)) & not(AP2(:,i)))|(not(EMF1(:,i)) & LFY(:,i) & not(AP1(:,i)))|(not(EMF1(:,i)) & not(AP2(:,i)) & LFY(:,i))|(not(EMF1(:,i)) & not(TFL1(:,i)) & LFY(:,i) & (AG(:,i) & SEP(:,i)))|(not(EMF1(:,i)) & (LFY(:,i) & WUS(:,i)));
    AP1(:,i+1)=(not(AG(:,i)) & not(TFL1(:,i)))|(FT(:,i) & LFY(:,i) & not(AG(:,i)))|(FT(:,i) & not(AG(:,i)) & not(PI(:,i)))|(LFY(:,i) & not(AG(:,i)) & not(PI(:,i)))|(FT(:,i) & not(AG(:,i)) & not(AP3(:,i)))|(LFY(:,i) & not(AG(:,i)) & not(AP3(:,i)));
    AP2(:,i+1)=not(TFL1(:,i));
    AP3(:,i+1)=(LFY(:,i) & UFO(:,i))|(PI(:,i) & SEP(:,i) & AP3(:,i) & (AG(:,i)|AP1(:,i)));
    EMF1(:,i+1)=not(LFY(:,i));
    FT(:,i+1)=(not(EMF1(:,i)))|u1(:,i);
    FUL(:,i+1)=not(AP1(:,i)) & not(TFL1(:,i));
    LFY(:,i+1)=not(EMF1(:,i))|not(TFL1(:,i));
    PI(:,i+1)=(LFY(:,i) & (AG(:,i)|AP3(:,i)))|(PI(:,i) & SEP(:,i) & AP3(:,i) & (AG(:,i)|AP1(:,i)));
    SEP(:,i+1)=LFY(:,i);
    TFL1(:,i+1)=not(AP1(:,i)) & (EMF1(:,i) & not(LFY(:,i)));
    UFO(:,i+1)=UFO(:,i);
    WUS(:,i+1)=WUS(:,i) & (not(AG(:,i))|not(SEP(:,i)));
    xt(i,:)=[AG(:,i+1) AP1(:,i+1) AP2(:,i+1) AP3(:,i+1) EMF1(:,i+1) FT(:,i+1) FUL(:,i+1) LFY(:,i+1) PI(:,i+1) SEP(:,i+1) TFL1(:,i+1) UFO(:,i+1) WUS(:,i+1)];
end
xt
xd=[AG(:,s+1) AP1(:,s+1) AP2(:,s+1) AP3(:,s+1) EMF1(:,s+1) FT(:,s+1) FUL(:,s+1) LFY(:,s+1) PI(:,s+1) SEP(:,s+1) TFL1(:,s+1) UFO(:,s+1) WUS(:,s+1)]

r1=bin2dec(num2str([x0; xt]));
%  stairs(r1,'-k')

subplot(2,1,1)
stairs(r1,'-ok')
grid on
ylabel('Attractor')
refline(0,260);%I1
% refline(0,262);%I2
% refline(0,263);%I3
% refline(0,261);%I4
refline(0,3240);%SE
% refline(0,3768);%PE1
% refline(0,5880);%ST1
% refline(0,5368);%CAR


subplot(2,1,2)
stairs(u1,'-ok')
grid on
ylabel('u_1(s)')

% subplot(4,1,3)
% stairs(u2,'-ok')
% grid on
% ylabel('u_2(s)')
% 
% subplot(4,1,4)
% stairs(u3,'-ok')
% grid on
% xlabel('s')
% ylabel('u_3(s)')

% subplot(5,1,5)
% stairs(u4,'-ok')
% grid on
% xlabel('s')
% ylabel('u_4(s)')


