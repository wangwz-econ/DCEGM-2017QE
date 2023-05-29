function [region_logical, region_index] = test_fun_detect(vec)
% Given a column vector vec, this function detects the strictly decreasing part.
% This function facilitates the refinement of endogenous wealth grid in DC-EGM algorithm.

% Input:
%   vec: a N by 1 vector which has the same property as the endogenous wealth grid in DC-EGM
% Output:
%   region_logical: a logical value. true means that there is a strictly decreasing part needing refinement.
%   region_index: a scalar denoting the start index of the strictly decreasing region. It equals [] if no such region.
lag = [NaN; vec(1:end-1)];
diff = vec - lag;
n = sum(diff<0);

if n>0
  region_logical = true;
  region_index = find(diff<0, 1, 'first');
else
  region_logical = false;
  region_index = [];
end

end