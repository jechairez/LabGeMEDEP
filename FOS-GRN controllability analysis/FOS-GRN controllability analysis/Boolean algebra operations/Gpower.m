function gpower=Gpower(L,m,s)
%This program computes a g0^s matrix
for i=1:s
    if (i>1)
        gpowern=anysp(ones(2^m,1),L);
        gpowerm=anysp(L,Bpower(gpowern,i-1));
        gpower=gpowerm;
    else
        gpower=double(L);
    end
end