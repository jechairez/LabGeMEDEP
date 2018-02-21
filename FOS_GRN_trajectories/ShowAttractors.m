% This function 
function [Att_Land,Att_Land_size]=ShowAttractors(L,k)
[n,l,c,r0,T] = bn(L,k);

fprintf('Number of attractors: %d\n\n',n);
fprintf('Lengths of attractors:\n');
disp(l);
k=1;
Att_Land=[];
for i=1:length(c)
    fprintf('\nAll attractors are displayed as follows:\n\n');
    fprintf('No. %d (length %d)\n\n',i,l(i));
    disp(c{i});
    Attr_check=c{i};
    size_ACh=size(Attr_check);
    size_AttCheck=size_ACh(1);
    for j=1:size_AttCheck
        if j==1
            Att_Land(k,:)=Attr_check(j,:);
            k=k+1;
        elseif Attr_check(j,:)~=Attr_check(j-1,:)
            Att_Land(k,:)=Attr_check(j,:);
            k=k+1;
        end
    end
end
Att_Land_sizexy=size(Att_Land);
Att_Land_size=Att_Land_sizexy(1);