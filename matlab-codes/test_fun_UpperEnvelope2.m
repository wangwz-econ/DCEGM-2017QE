function [refined_func] = test_fun_UpperEnvelope2(func, varargin)
% This function conducts the refinement of the endogenous wealth grid!
% The specific process is described in the Algorithm 3: Upper Envelopment from the original QE paper.

% Inputs:
%   func: N by 2 matrix, 
%         first column denoting endogenous wealth grid, 
%         second column denoting choice-specific value function or consumption function
% Optional Inputs:
%   InterpMethod: character, the interpolation method used to extrapolate 
%                 the choice-specific value function (default is 'spline')
%   gridN: integer, number of points used to interpolate the overlapping region to compare 
%          choice-specific value function (defalue is 100)
%   insert: logical, 
%           if true, then insert a kink point where two regions of value function intersect
% Outputs:
%   refined_func: N1 by 2 matrix (N1 is not necessarily less than N because overlapping region 
%                                 in the endogenous wealth grid are interpolated many values)
%                 choice-specific value or consumption function after refining the wealth grid


n = nargin;
if n==1
  InterpMethod = 'linear';  gridN = 100; 
  insert = false;
elseif n==2
  InterpMethod = varargin{1};  gridN = 100;
  insert = false;
elseif n==3
  InterpMethod = varargin{1};  gridN = varargin{2};
  insert = false;
elseif n==4
  InterpMethod = varargin{1};  gridN = varargin{2};
  insert = true;
else
  error 'Wrong Input Arguments in Function fun_UpperEnvelope'
end

wealth = func(:, 1); value = func(:, 2);

[~, start] = test_fun_detect(wealth);
for step = 1:(length(wealth) - start - 1)
  increasing = (wealth(step+start) > wealth(step+start-1));
  if increasing, break;  end
end

wealth1 = wealth(1:start-1); value1 = value((1:start-1));
wealth3 = wealth((start+step):end); value3 = value((start+step):end);

line1 = [wealth(1:start-1), value(1:start-1)];
line2 = [wealth((start-1):(start+step-1)), value((start-1):(start+step-1))];
line3 = [wealth((start+step):end), value((start+step):end)];

overlap_min = wealth(start);
overlap_max = wealth(start-1);
overlap = linspace(overlap_min, overlap_max, gridN)';

extr_1 = interp1(line1(:,1), line1(:,2), overlap, InterpMethod, 'extrap');
extr_2 = interp1(line2(:,1), line2(:,2), overlap, InterpMethod, 'extrap');
extr_3 = interp1(line3(:,1), line3(:,2), overlap, InterpMethod, 'extrap');
[extr_max] = max([extr_1, extr_2, extr_3], [], 2, 'includenan');

wealth = [wealth1(wealth1<overlap_min); overlap; wealth3(wealth3>overlap_max)];
value  = [value1(wealth1<overlap_min); extr_max; value3(wealth3>overlap_max)];

refined_func = [wealth, value];


end