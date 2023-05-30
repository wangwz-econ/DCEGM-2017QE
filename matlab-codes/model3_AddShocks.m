%% DESCRIPTION
% DC-EGM for the model with income and taste shocks
% The differences from model2_deterministic_retirement.m are:
%   The Euler equation is more complicated (integration over income shock and use next-period CCP).
%   There is no need to construct the upper envelope of choice-specific value functions; 
%     instead, compute the conditional choice probabilities.
%   Perform integration to calculate the expectation of the income process.
% Still, we need to refine the endogenous wealth grid to delete the suboptimal solutions.

% The parameters are taken from the note of Figure 4 in the original paper.

addpath("D:\Matlab_files\CompEcon2016\CEtools") 
% We need to use the toolbox provided with Applied Computational Economics and Finance.
% This toolbox is downloaded from http://compeconworkshop.org/.
% It is used to get the quadrature points and weights for a normal distribution.
% Specifically, qnwnorm function is what we need!

%% Step 1: Set up Parameters
clear all

% Economic Parameters
beta = 0.98; R = 1; y = 20; T = 20; delta = 1;
eta_var =  0.05;     % sigma_eta^2, variance of the income shock
eta_mu = -eta_var/2;  % mean of the income shock
epsilon = 0.15;
u_D = @ (c) 1./c;     % marginal utility function
u_D_inv = @ (u) 1./u; % the inverse of marginal utility function

% Computational Parameters
InterpMethod = 'spline';
ExtrapMethod = 'spline';
gridN = 2e+3;
quadN = 9;
wealth_min = 1e-5; wealth_max = 400;
savings = linspace(wealth_min, wealth_max, gridN)'; % exogenous savings grid (gridN by 1)
wealth_grid_r = R * savings;           % wealth grid if the agent chooses to retire (gridN by 1)

[log_eta, wgt_log_eta] = qnwnorm(quadN, eta_mu, eta_var); % both are column vectors
wealth_grid_w = repmat(R*savings,1,quadN) + repmat(transpose(y*exp(log_eta)), gridN, 1);
% wealth grid if the agent chooses to work (gridN by quadN)
% wealth_grid_w * wgt_log_eta can get a gridN by 1 vector which is the expectation of y*eta

policy_w = cell(T,1);
policy_r = cell(T,1);
% Variable policy_j, j = w, r is a T by 1 cell.
% In each cell t, it stores a matrix denoting the worker's or retiree's value and policy function.
%   In each matrix of a cell of policy_r and, there are 3 columns:
%     the 1st column is the endogenous wealth grid,
%     the 2nd column is optimal consumption given the savings grid,
%     the 3rd column is unconditional value function.
%   policy_w has a 4 columns:
%     the 1st column is the endogenous wealth grid,
%     the 2nd column is optimal consumption given the savings grid,
%     the 3th column is the ex ante value function,
%     the 4th column is the CCP that d_t=1.

choice_r = cell(T,1);
% Variable choice_r is a T by 1 cell.
% In each cell t, it stores a matrix denoting the worker's period t retirement-specific functions.
% To be specific, in each cell entry of choice_r, there are 3 columns:
%   the 1st column is choice-specific optimal consumption given choice d_t=0 (retire),
%   the 2nd column is the endogenous wealth grid given choice d_t=0 (retire),
%   the 3rd column is the choice-specific value function given choice d_t=0 (retire).

choice_w = cell(T,1);
% Variable choice_w is a T by 1 cell.
% In each cell t, it stores three matrices denoting the worker's period t working-specific functions.
% To be specific, in each cell entry of choice_w, there are 3 columns:
%   the 1st column is choice-specific optimal consumption given choice d_t=1 (keep working),
%   the 2nd column is the endogenous wealth grid given choice d_t=1 (keep working),
%   the 3rd column is the choice-specific value function given choice d_t=1 (keep working).
choice_w_refined = cell(T, 1);
% Store the refined choice-specific consumption and value functions.
%% Step 2: EGM Algorithm
% The most important step is to refine the endogenous wealth grids (2nd and 5th column)!
% Because a solution to the Euler equation is not necessary a solution that maximizes the Bellman equation!

for t = T:-1:1
  fprintf('t = %i', t)
  if t == T
    policy_r{t}(:,1) = savings;                  % the endogenous wealth grid
    policy_r{t}(:,2) = savings;                  % optimal consumption given the savings grid
    policy_r{t}(:,3) = log(policy_r{t}(:,1));    % unconditional value function

    policy_w{t}(:,1) = savings;                  % the endogenous wealth grid
    policy_w{t}(:,2) = NaN(gridN, 1);            % optimal consumption given the savings
    policy_w{t}(:,3) = epsilon*( log( exp(delta/epsilon) + 1 ) + (log(savings) - delta)./epsilon ); % ev
    policy_w{t}(:,4) = repmat(1/(1+exp(delta/epsilon)), gridN, 1);  % CCP when d_t=1
    policy_w{t}(:,5) = 1 - policy_w{t}(:,4);      % CCP when d_t=0

    choice_r{t}(:,1) = policy_r{t}(:,2);
    choice_r{t}(:,2) = policy_r{t}(:,1);
    choice_w_refined{t}(:,1) = savings;
    choice_w_refined{t}(:,2) = savings;

  else
    rhs_r = beta * R * u_D(interp1(policy_r{t+1}(:,1), policy_r{t+1}(:,2), wealth_grid_r, InterpMethod, "extrap"));

    % Case 1: The retiree's problem.
    policy_r{t}(:,2) = u_D_inv(rhs_r);
    policy_r{t}(:,1) = savings + policy_r{t}(:,2);
    policy_r{t}(:,3) = log(policy_r{t}(:,2)) + ...
      beta * interp1(policy_r{t+1}(:,1), policy_r{t+1}(:,3), wealth_grid_r, InterpMethod, "extrap");

    % Case 2: The worker's problem.
    % Step 1: Collect choice-specific optimal mappings in variable choice_w using Euler equations.

    % If the worker chooses to retire in period t:
    choice_r{t}(:,1) = policy_r{t}(:,2);
    choice_r{t}(:,2) = policy_r{t}(:,1);
    choice_r{t}(:,3) = policy_r{t}(:,3);

    % If the worker chooses to retire in period t:
    % to solve the smoothed euler equation, we need to collect t+1 period choice-specific consumption function and t+1 period ccp
    consum_w_tplus1 = griddedInterpolant(choice_w_refined{t+1}(:,2), choice_w_refined{t+1}(:,1), InterpMethod, ExtrapMethod);
    consum_r_tplus1 = griddedInterpolant(choice_r{t+1}(:,2), choice_r{t+1}(:,1), InterpMethod, ExtrapMethod);
    ccp_w_tplus1 = griddedInterpolant(policy_w{t+1}(:,1), policy_w{t+1}(:,4), InterpMethod, ExtrapMethod);
    ccp_r_tplus1 = griddedInterpolant(policy_w{t+1}(:,1), policy_w{t+1}(:,5), InterpMethod, ExtrapMethod);

    rhs_w = beta * R * ...
      (ccp_w_tplus1(wealth_grid_w)./consum_w_tplus1(wealth_grid_w) + ccp_r_tplus1(wealth_grid_w)./consum_r_tplus1(wealth_grid_w)) * wgt_log_eta; % gridN by quadN times quadN by 1
    choice_w{t}(:,1) = u_D_inv(rhs_w);
    choice_w{t}(:,2) = choice_w{t}(:,1) + savings;

    ev_tplus1 = griddedInterpolant(policy_w{t+1}(:,1), policy_w{t+1}(:,3), InterpMethod, ExtrapMethod);
    choice_w{t}(:,3) = log(choice_r{t}(:,1)) - delta + beta * ev_tplus1(wealth_grid_w) * wgt_log_eta;

    % Add credit constraint
    choice_w{t} = [NaN(1, size(choice_w{t}, 2)); choice_w{t}];
    choice_w{t}(1,[1 2]) = 0;
    choice_r{t} = [NaN(1, size(choice_r{t}, 2)); choice_r{t}];
    choice_r{t}(1,[1 2]) = 0;

    % Step 2: Refine the endogenous wealth grid for the worker's problem when choosing to keep working.
    value_w = cls_line(choice_w{t}(2:end,2), choice_w{t}(2:end,3)); % choice-specific value function when d_t=1
    consu_w = cls_line(choice_w{t}(:,2), choice_w{t}(:,1)); % choice-specific optimal consumption when d_t=1

    [value_w_refined, indx_del, new_dots] = cls_outer_refinement(value_w);
    consu_w_refined = cls_thinout(consu_w, indx_del+1); % delete the suboptimal solutions
    new_wealth = new_dots.x;
    new_consum = cls_interp(consu_w_refined, new_wealth);
    consu_w_refined = cls_grow(consu_w_refined, cls_line(new_wealth, new_consum), false);
    consu_w_refined = cls_sort(consu_w_refined);

    choice_w_refined{t}(:,2) = consu_w_refined.x;
    choice_w_refined{t}(:,1) = consu_w_refined.y;
    choice_w_refined{t}(:,3) = [NaN; value_w_refined.y];

    % Step 3: Collect the unconditional values
    policy_w{t}(:,1) = savings;

    % Compute the conditional choice probabilities
    csv_w_fun = griddedInterpolant(value_w_refined.x, value_w_refined.y, InterpMethod, ExtrapMethod);
    csv_r_fun = griddedInterpolant(choice_r{t}(2:end,2), choice_r{t}(2:end,3), InterpMethod, ExtrapMethod);

    csv_w = csv_w_fun(policy_w{t}(:,1)) ./ epsilon;
    csv_r = csv_r_fun(policy_w{t}(:,1)) ./ epsilon;

    normalization = max([csv_w, csv_r], [], 2);

    policy_w{t}(:,3) = epsilon * log( exp(csv_w - normalization) + exp(csv_r - normalization) ) + normalization;
    policy_w{t}(:,4) = exp(csv_w - normalization) ./ (exp(csv_w - normalization) + exp(csv_r - normalization));
    policy_w{t}(:,5) = 1 - policy_w{t}(:,4); 

  end % end the if statement to check if the current period less than T
  fprintf("  Calculation Successful!\n")
end % end the for loop that goes through t backward

%% Step 3: Draw the Optimal Consumption Functions in Different Periods

consum = choice_w_refined{15}(:,1); 
wealth = choice_w_refined{15}(:,2); 
consum = consum(wealth<120 & wealth>10); wealth = wealth(wealth<120 & wealth>10);
plot(wealth, consum)
