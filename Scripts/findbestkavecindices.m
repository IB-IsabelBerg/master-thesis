function indexnumber_in_kavec1 = findbestkavecindices(kavec1,kavec2)
% This function gives the index numbers in kavec1 such that:
%
% the selected indices in kavec1 are as close as
% possible to the kavec2-values.
%
% Example:
%    kavec1 = [1   2   3   4 5];
%    kavec2 = [1 2.8 4.2 5.6 7.0];
%
% For this, we want:
% indexnumber_in_kavec1 = [1 3 4 5 5]; 
% 
% indexnumber_in_kavec1 = findbestkavecindices(kavec1,kavec2);

nvals = length(kavec2);

indexnumber_in_kavec1 = zeros(1,nvals);

for ii = 1:nvals
   absdistance = abs(kavec1-kavec2(ii));
   [~,minloc] = min(absdistance);
   indexnumber_in_kavec1(ii) = minloc(1);
end

    
    
end

