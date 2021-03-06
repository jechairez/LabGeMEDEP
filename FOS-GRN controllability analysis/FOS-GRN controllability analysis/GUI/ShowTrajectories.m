function [Trajectories_a,Trajectories_b]=ShowTrajectories(res2,Att_Land_size)
Trajectories_a{Att_Land_size,Att_Land_size}=[];
Trajectories_b{Att_Land_size,Att_Land_size}=[];
Att_names_x=cell(1,length(Att_Land_size));
for i=1:length(res2)
    x=res2(i,3);
    y=res2(i,4);
    switch res2(i,2)
        case 1
            Z='AG';
        case 2
            Z='AP1';
        case 3
            Z='AP2';
        case 4
            Z='AP3';
        case 5
            Z='EMF1';
        case 6
            Z='FT';
        case 7
            Z='FUL';
        case 8
            Z='LFY';
        case 9
            Z='PI';
        case 10
            Z='SEP';
        case 11
            Z='TFL1';
        case 12
            Z='UFO';
        case 13
            Z='WUS';
    end
    Trajectories_a{x,y}=strjoin({Trajectories_a{x,y},Z},'\');
    Trajectories_b{x,y}=strjoin({Trajectories_b{x,y},num2str(res2(i,1))},'\');
end
Trajectories_a
Trajectories_b