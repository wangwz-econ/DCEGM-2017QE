function [res_x, res_y, index] = fun_UpperEnvelope(x1, y1, x2, y2, varargin)
% Given two functions denoted by x1, y1 and x2, y2 (four column vectors), this function returns the 
%   upper envelope of these two functions (denoted by res_x, res_y)

% Inputs:
%   x1,y1: N1 by 1 vectors, denoting a function f1
%   x2,y2: N2 by 1 vectors, denoting another function f2
% Other possible inputs (in order):
%   InterpMethod: character, denoting the interpolation method (default: 'linear')
%   N: the number of points used to denote the upper envelope function (default: 2000)
% Outputs:
%   res_x,res_y: N by 1 vectors, denoting the upper envelope of functions f1 and f2
%   index: N by 1 vector consisting of 1 and 2, denoting whether the envelope is from f1 or f2 

n = nargin;
if n==4
  InterpMethod = 'linear';
  N = 2000;
elseif n==5
  InterpMethod = varargin{1};
  N = 2000;
elseif n==6
  InterpMethod = varargin{1};
  N = varargin{2};
end


res_x = linspace(min([min(x1), min(x2)]), max([max(x1), max(x2)]), N)';
y1_re = interp1(x1, y1, res_x, InterpMethod, 'extrap');
y2_re = interp1(x2, y2, res_x, InterpMethod, 'extrap');
[res_y, index] = max([y1_re, y2_re], [], 2, 'includenan');

end