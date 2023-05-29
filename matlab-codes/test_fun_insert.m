function [res_line] = test_fun_insert(line, index, x, y)
% This function inserts a point (x,y) at the index position to the original line and returns res_line.

x_ori = line(:, 1);
y_ori = line(:, 2); 

if length(x_ori) < index
  x_res = [x_ori; x]; 
  y_res = [y_ori; y]; 
  res_line = [x_res, y_res];
elseif index > 1
  x_line1 = x_ori(1:(index-1)); x_line2 = x_ori(index:end);
  y_line1 = y_ori(1:(index-1)); y_line2 = y_ori(index:end);
  x_res = [x_line1; x; x_line2];
  y_res = [y_line1; y; y_line2];
  res_line = [x_res, y_res];
elseif index==1
  x_res = [x; x_ori]; 
  y_res = [y; y_ori]; 
  res_line = [x_res, y_res];

end