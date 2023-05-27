% DC-EGM for the deterministic model
% The only difficulty is to construct the upper envelope of the endogenous grids and to 
%   interpolate a discountinuous and non-smooth function.

% The parameters are taken from the note of Figure 2 in the original paper.

clear all

% Economic Parameters
R = 1; beta = 0.98; y = 20; T = 20; delta = 1; 
u_D = @ (c) 1./c;     % marginal utility function
u_D_inv = @ (u) 1./u; % the inverse of marginal utility function

% Computational Parameters
gridN = 400; wealth_max = 400; wealth_min = 1e-5; option = optimset('TolX', 1e-5);
savings = linspace(wealth_min, wealth_max, gridN)'; % exogenous savings grid (row vector, 1 by gridN)
wealth_grid_w = R * savings + y;       % wealth grid if the agent chooses to work (1 by gridN)
wealth_grid_r = R * savings;           % wealth grid if the agent chooses to retire (1 by gridN)


% Variable policy_j, j = w, r is a T by 1 cell. 
% In each cell t, it stores a matrix denoting the worker's or retiree's value and policy function.
%   In each matrix of policy_r and policy_w, there are 3 columns:
%     the 1st column is the endogenous wealth grid,
%     the 2nd column is optimal consumption given the savings grid,
%     the 3rd column is unconditional value function.
policy_w = cell(T,1);
policy_r = cell(T,1);

% Variable choice_w is a T by 1 cell.
% In each cell t, it stores a matrix denoting the worker's period t choice-specific functions.
% To be specific, in each matrix of choice_w, there are 6 columns:
%   the 1st column is choice-specific optimal consumption given choice d_t=1 (keep working), 
%   the 2nd column is choice-specific optimal consumption given choice d_t=0 (retire), 
%   the 3rd column is the endogenous wealth grid given choice d_t=1 (keep working),
%   the 4th column is the endogenous wealth grid given choice d_t=0 (retire),
%   the 5th column is the choice-specific value function given choice d_t=1 (keep working),
%   the 6th column is the choice-specific value function given choice d_t=0 (retire).
choice_w = cell(T,1);

% The most important step is to obtain the upper envelope of the 5th and 6th column of choice_w for
%   each backward stepping t!
% Because we need to choose the maximum choice-specific value function across (keep working) and (retire) 
%   as our unconditional value funcion for the worker.

for t = T:-1:18
  if t == T
    policy_r{t}(:,1) = savings;
    policy_r{t}(:,2) = savings;
    policy_r{t}(:,3) = log(savings); % log(policy_r{t}(:,1))

    policy_w{t}(:,1) = savings;
    policy_w{t}(:,2) = savings;
    policy_w{t}(:,3) = log(savings); % log(policy_w{t}(:,1))
  else
    % Case 1: The retiree's problem.
    rhs_r = beta * R * u_D(interp1(policy_r{t+1}(:,1), policy_r{t+1}(:,2), wealth_grid_r, "linear", "extrap"));
    policy_r{t}(:,2) = u_D_inv(rhs_r);
    policy_r{t}(:,1) = savings + policy_r{t}(:,2);
    policy_r{t}(:,3) = log(policy_r{t}(:,2)) + ...
            beta * interp1(policy_r{t+1}(:,1), policy_r{t+1}(:,3), R*savings, "spline", "extrap");

    % Case 2: The worker's problem.
    % Step 1: Collect choice-specific optimal mappings in variable choice_w using Euler equations.
    rhs_w_w = beta * R * u_D(interp1(policy_w{t+1}(:,1), policy_w{t+1}(:,2), wealth_grid_w, "linear", "extrap"));
    rhs_w_r = beta * R * u_D(interp1(policy_r{t+1}(:,1), policy_r{t+1}(:,2), wealth_grid_r, "linear", "extrap"));

    choice_w{t}(:,1) = u_D_inv(rhs_w_w);
    choice_w{t}(:,2) = u_D_inv(rhs_w_r);
    choice_w{t}(:,3) = choice_w{t}(:,1) + savings;
    choice_w{t}(:,4) = choice_w{t}(:,2) + savings;
    choice_w{t}(:,5) = log(choice_w{t}(:,1)) - delta + ...
      beta * interp1(policy_w{t+1}(:,1), policy_w{t+1}(:,3), R*savings+y, "spline", "extrap");
    choice_w{t}(:,6) = log(choice_w{t}(:,2)) + ...
      beta * interp1(policy_r{t+1}(:,1), policy_r{t+1}(:,3), R*savings, "spline", "extrap");

    % Add credit constraint
    choice_w{t} = [NaN(1, size(choice_w{t}, 2)); choice_w{t}];
    choice_w{t}(1,1:4) = 0;

    % Step 2: Find the upper envelope of the choice-specific value functions

    % Step 2.1: First we need to detect regions where multiple FOCs are possible 
    %           based on the non-monotonicity of endogenous wealth grid.

    region = fun_detect(choice_w{t}(:,3));
    if ~region
      % Now the choice-specific endogenous wealth grid is non-decreasing, 
      %   we only need to contruct the upper envelope of two choice-specific value function!
      [temp_wealth, temp_value, index] = fun_UpperEnvelope(choice_w{t}(:,3), choice_w{t}(:,5), choice_w{t}(:,4), choice_w{t}(:,6), 'linear', gridN);
      temp_consum_w = interp1(choice_w{t}(:,3), choice_w{t}(:,1), temp_wealth, 'linear', 'extrap');
      temp_consum_r = interp1(choice_w{t}(:,4), choice_w{t}(:,2), temp_wealth, 'linear', 'extrap');
      temp_consum = [temp_consum_w(index==1); temp_consum_r(index==2)];

      policy_w{t}(:,1) = temp_wealth; 
      policy_w{t}(:,2) = temp_consum; 
      policy_w{t}(:,3) = temp_value; 
    else
      % Now the choice-specific endogenous wealth grid has non-monotonic part,
      %   we need first to delete the non-optimal region in the 

      temp_w_w = choice_w{t}(:,3); temp_c_w = choice_w{t}(:,1); temp_v_w = choice_w{t}(:,5);

      for i=1:size(region, 1)
        temp_w_w(region(i,1)-1:region(i,2)) = NaN;
        temp_c_w(region(i,1)-1:region(i,2)) = NaN;
        temp_v_w(region(i,1)-1:region(i,2)) = NaN;
        temp_w_w(isnan(temp_v_w)) = NaN;
        temp_c_w(isnan(temp_v_w)) = NaN;
      end
      temp_w_w(isnan(temp_w_w))=[]; temp_v_w(isnan(temp_v_w))=[]; temp_c_w(isnan(temp_c_w))=[]; 
      [temp_wealth, temp_value, index] = fun_UpperEnvelope(temp_w_w, temp_v_w, choice_w{t}(:,4), choice_w{t}(:,6), 'linear', gridN);
      temp_consum_w = interp1(temp_w_w, temp_c_w, temp_wealth, 'linear', 'extrap');
      temp_consum_r = interp1(choice_w{t}(:,4), choice_w{t}(:,2), temp_wealth, 'linear', 'extrap');
      temp_consum = [temp_consum_w(index==1); temp_consum_r(index==2)];

      policy_w{t}(:,1) = temp_wealth;
      policy_w{t}(:,2) = temp_consum;
      policy_w{t}(:,3) = temp_value;
      
    end

  end
end

plot(policy_w{18}(:,1), policy_w{18}(:,2))
plot(policy_w{1}(:,1), policy_w{1}(:,2))
plot(policy_w{10}(:,1), policy_w{10}(:,2))