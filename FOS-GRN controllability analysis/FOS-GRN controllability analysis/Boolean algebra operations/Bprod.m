function [C]=Bprod(A,B)
%This program calculates a Boolean product of two matrices mxn and nxp.
sz1=size(A);
sz2=size(B);

for k=1:sz1(2)
    for i=1:sz1(1)
        for j=1:sz2(2)
            if (k>1)
               Cm(i,j)=A(i,k)&B(k,j);
               Cn(i,j)=Bsum(C(i,j),Cm(i,j));
               C(i,j)=Cn(i,j);
            else
                C(i,j)=A(i,k)&B(k,j); 
            end
        end
    end
end
