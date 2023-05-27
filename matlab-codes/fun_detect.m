function [indices] = fun_detect(vec)
% Given a column vector vec, this function detects the non-decreasing part.
% This function facilitates the refinement of endogenous wealth grid in DC-EGM algorithm.

% Input:
%   vec: a N by 1 vector which has the same property as the endogenous wealth grid in DC-EGM
% Output:
%   indices: a n by 2 matrix, each row indicates a region (by the first and last indices) 
%            where multiple FOCs are possible.
lag = [NaN; vec(1:end-1)];
diff = vec - lag;
n = sum(diff<0);

if n>0
  indices = NaN(n, 2);
  indices(:,1) = transpose(find(diff<0));
  for i =1:n
    res_1 = vec(1:(indices(i)-1));
    res_2 = vec(indices(i):end);
    index = min(find(res_2 > res_1(end)));
    indices(i,2) = indices(i) + index - 1;
  end

else
  indices = false;
end




end