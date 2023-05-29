function [refined_value_cond, index_remaining, varargout] = test_fun_UpperEnvelope(value_cond, varargin)
% This function conducts the refinement of the endogenous wealth grid!
% The specific process is described in the Algorithm 3: Upper Envelopment from the original QE paper.

% Inputs:
%   value_cond: N by 2 matrix, 
%               first column denoting endogenous wealth grid, 
%               second column denoting choice-specific value function
% Optional Inputs:
%   InterpMethod: character, the interpolation method used to extrapolate 
%                 the choice-specific value function (default is 'spline')
%   gridN: integer, number of points used to interpolate the overlapping region to compare 
%          choice-specific value function (defalue is 2000)
%   insert: logical, 
%           if true, then insert a kink point where two regions of value function intersect
% Outputs:
%   refined_value_cond: N1 by 2 matrix (N1<N)
%                       choice-specific value function after refining the endogenous wealth grid
%   index_remaining: N1 by 1 vector
%                    indices (row indices) of remaining points in the original value_cond
% Optional Outputs:
%   varargout{1}: integer, the index in the wealth grid that a kink point should be inserted
%   varargout{2}: double, the wealth value of kink point 
%   varargout{3}: double, the choice-specific value function value of kink point

n = nargin;
if n==1
  InterpMethod = 'spline';  gridN = 100; 
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

wealth = value_cond(:, 1); value = value_cond(:, 2);

[~, start] = test_fun_detect(wealth);
for step = 1:(length(wealth) - start - 1)
  increasing = (wealth(step+start) > wealth(step+start-1));
  if increasing, break;  end
end

wealth1 = wealth(1:start-1); 
wealth3 = wealth((start+step):end); 

line1 = [wealth(1:start-1), value(1:start-1)];
line2 = [wealth((start-1):(start+step-1)), value((start-1):(start+step-1))];
line3 = [wealth((start+step):end), value((start+step):end)];

overlap_min = wealth(start); overlap_max = wealth(start-1);
overlap = linspace(overlap_min, overlap_max, gridN)';

extr_1 = interp1(line1(:,1), line1(:,2), overlap, InterpMethod, 'extrap');
extr_2 = interp1(line2(:,1), line2(:,2), overlap, InterpMethod, 'extrap');
extr_3 = interp1(line3(:,1), line3(:,2), overlap, InterpMethod, 'extrap');
[~, max_index] = max([extr_1, extr_2, extr_3], [], 2, 'includenan');

wealth_threshold_index = find(max_index==3, 1, 'first');
if ~isempty(wealth_threshold_index)
  wealth_threshold = overlap(wealth_threshold_index);
elseif isempty(wealth_threshold_index)
  wealth_threshold = overlap_max;
end

index_remaining = [find(wealth1<=wealth_threshold); find(wealth3>wealth_threshold)+start-step+1];
refined_value_cond = [wealth(index_remaining), value(index_remaining)];

if insert % no need to insert the kink point

  fun_kink = @(x) interp1(line1(:,1), line1(:,2), x, InterpMethod, 'extrap') - ...
                  interp1(line3(:,1), line3(:,2), x, InterpMethod, 'extrap');
  options = optimset('TolX', 1e-10);
  kink_wealth = fzero(fun_kink, wealth_threshold, options);
  kink_value = interp1(line1(:,1), line1(:,2), kink_wealth, InterpMethod, 'extrap');

  varargout{1} = find(wealth1<=wealth_threshold, 1, 'last') + 1;
  varargout{2} = kink_wealth;
  varargout{3} = kink_value;


end


end