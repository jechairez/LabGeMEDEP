function Mc=Mcontrol(L,n,m)
%This program computes the controllability matrix
M=anysp(double(L),ones(2^m,1));
for s=1:2^(m+n)
    if (s>1)
        Mc=Bsum(Mc,Bpower(M,s));
    else
        Mc=double(M);
    end
end