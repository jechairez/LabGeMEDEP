function Cpower=Bpower(A,n)
%This program computes a n-power of a matrix
for i=1:n
    if (i>1)
        Cpowern=Bprod(Cpower,A');
        Cpower=Cpowern;
    else
        Cpower=A;
    end
end
