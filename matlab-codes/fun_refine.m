function [res, index] = fun_refine(vector)

% This function deletes points in a column vector to produce a non-decreasing resulted vector.
% Important: The original vector should have a specific structure as in the DC-EGM QE paper.

% Inputs: 
%   vector, a column vector (generally non-decreasing, except for one discontinuous jump downward)
% Output: 
%   res, a column vector (its elements are non-decreasing!)
%   index, a row vector, storing the indices of remaining elements in the input vector 

lag = [NaN;vector(1:end-1)];
diff = vector - lag;
if all(diff>=0 | isnan(diff))
  res = vector;
  index = 1:length(vector);
else
  index_start = find(diff<0);
  temp = vector(index_start:end, :);
  res_1 = vector(1:index_start-1);
  res_2 = temp(temp>res_1(end));

  res = [res_1; res_2];
  deleted_index = index_start:1:(index_start + max(find(temp<=res_1(end)) - 1));
  index = setdiff(1:length(vector), deleted_index);
end

end