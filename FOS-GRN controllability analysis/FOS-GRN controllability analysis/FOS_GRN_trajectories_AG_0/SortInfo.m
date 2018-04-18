function [res2]=SortInfo(res)
gg=1;
for fg=1:length(res)
    if (res(fg,3)~=res(fg,4))
        res2(gg,:)=res(fg,:);
        gg=gg+1;
    end
end