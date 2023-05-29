function [res_fun, index] = fun_Choice2(func1, func2, grid_ends, varargin)
% Given two functions denoted by func1 and func2, 
%   this function returns the upper envelope of these two functions (denoted by res_fun).
% This is used to construct worker's choice by considering the upper envelope of the choice-
%   specific conditional value function.

% Inputs:
%   func1: N1 by 2 vectors, denoting a function f1
%   func2: N2 by 2 vectors, denoting another function f2
%   grid_ends: 1 by 2 vector, storing the minimum and maximum endpoints in the interpolation
% Other possible inputs (in order):
%   InterpMethod: character, denoting the interpolation method (default: 'linear')
%   N: the number of points used to denote the upper envelope function (default: 2000)
% Outputs:
%   res_fun: N by 2 vectors, denoting the upper envelope of functions f1 and f2
%   index: N by 1 vector consisting of 1 and 2, denoting whether the envelope is from f1 or f2 

grid_min = grid_ends(1);
grid_max = grid_ends(2);
n = nargin;
if n==3
  InterpMethod = 'linear';
  N = 2000;
elseif n==4
  InterpMethod = varargin{1};
  N = 2000;
elseif n==5
  InterpMethod = varargin{1};
  N = varargin{2};
else
  error 'Wrong Input Arguments in Function fun_Choice'
end

res_x = linspace(grid_min, grid_max, N)';
y1_re = interp1(func1(:,1), func1(:,2), res_x, InterpMethod, 'extrap');
y2_re = interp1(func2(:,1), func2(:,2), res_x, InterpMethod, 'extrap');
[res_y, index] = max([y1_re, y2_re], [], 2, 'includenan');

res_fun = [res_x, res_y];

end