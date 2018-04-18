function [C]=Bsum(A,B)
%This program calculates a Boolean n-sum of two matrices mxn.
%Chapter 11. Stability and Stabilization. Definition 11.2
sz=size(A);

for i=1:sz(1)
    for j=1:sz(2)
        C(i,j)=A(i,j)|B(i,j);
    end
end
