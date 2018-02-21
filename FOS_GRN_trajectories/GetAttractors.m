function [expr,vars]=GetAttractors(A,options,k,ME,Mi,MD,MN,MR,MC,MU,MX)
A = completeeqn(A,options); % make each equation have all variables
[expr,vars] = stdform(A,options,k);  % get the canonical form for each variable
for i=1:length(expr)            % calculate the structure matrix for each variable
    expr{i} = eval(expr{i});
end 