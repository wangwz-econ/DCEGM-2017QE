% DC-EGM for the deterministic model
% The only difficulty is to refine endogenous wealth grid by first detecting the non-motonotic region 
%   and second deleting certain points through the upper enevelope of choice-specific value function
%   in that region. 
% The function fun_UpperEnvelope is responsible for this step!

% The parameters are taken from the note of Figure 2 in the original paper.

%% Step 1: Set up Parameters
clear all

% Economic Parameters
R = 1; beta = 0.98; y = 20; T = 20; delta = 1; 
u_D = @ (c) 1./c;     % marginal utility function
u_D_inv = @ (u) 1./u; % the inverse of marginal utility function

% Computational Parameters
gridN = 2e+3; InterpMethod = 'pchip';
wealth_min = 1e-5; wealth_max = 400; 
savings = linspace(wealth_min, wealth_max, gridN)'; % exogenous savings grid (gridN by 1)
wealth_grid_w = R * savings + y;       % wealth grid if the agent chooses to work (gridN by 1)
wealth_grid_r = R * savings;           % wealth grid if the agent chooses to retire (gridN by 1)

% Variable policy_j, j = w, r is a T by 1 cell. 
% In each cell t, it stores a matrix denoting the worker's or retiree's value and policy function.
%   In each matrix of a cell of policy_r and policy_w, there are 3 columns:
%     the 1st column is the endogenous wealth grid,
%     the 2nd column is optimal consumption given the savings grid,
%     the 3rd column is unconditional value function.
policy_w = cell(T,1);
policy_r = cell(T,1);

% Variable choice_w is a T by 1 cell.
% In each cell t, it stores a matrix denoting the worker's period t choice-specific functions.
% To be specific, in each matrix of choice_w, there are 6 columns:
%   the 1st column is choice-specific optimal consumption given choice d_t=1 (keep working), 
%   the 2nd column is the endogenous wealth grid given choice d_t=1 (keep working),
%   the 3rd column is the choice-specific value function given choice d_t=1 (keep working),
%   the 4th column is choice-specific optimal consumption given choice d_t=0 (retire), 
%   the 5th column is the endogenous wealth grid given choice d_t=0 (retire),
%   the 6th column is the choice-specific value function given choice d_t=0 (retire).
choice_w = cell(T,1);

%% Step 2: EGM Algorithm
% The most important step is to refine the endogenous wealth grids (2nd and 5th column)!
% The underlying reason is that a solution to the Euler equation is not necessary the solution that 
%   maximizes the Bellman equation!

for t = T:-1:1
  fprintf('t = %i\n', t)
  if t == T
    policy_r{t}(:,1) = savings;                  % the endogenous wealth grid
    policy_r{t}(:,2) = savings;                  % optimal consumption given the savings grid
    policy_r{t}(:,3) = log(policy_r{t}(:,1));    % unconditional value function

    policy_w{t}(:,1) = savings;                  % the endogenous wealth grid
    policy_w{t}(:,2) = savings;                  % optimal consumption given the savings grid
    policy_w{t}(:,3) = log(policy_w{t}(:,1));    % unconditional value function
  else
    rhs_r = beta * R * u_D(interp1(policy_r{t+1}(:,1), policy_r{t+1}(:,2), wealth_grid_r, 'linear', "extrap"));
    rhs_w = beta * R * u_D(interp1(policy_w{t+1}(:,1), policy_w{t+1}(:,2), wealth_grid_w, 'linear', "extrap"));

    % Case 1: The retiree's problem.
    policy_r{t}(:,2) = u_D_inv(rhs_r);
    policy_r{t}(:,1) = savings + policy_r{t}(:,2);
    policy_r{t}(:,3) = log(policy_r{t}(:,2)) + ...
            beta * interp1(policy_r{t+1}(:,1), policy_r{t+1}(:,3), wealth_grid_r, 'linear', "extrap");

    % Case 2: The worker's problem.
    % Step 1: Collect choice-specific optimal mappings in variable choice_w using Euler equations.

    % If the worker chooses to keep working in period t:
    choice_w{t}(:,1) = u_D_inv(rhs_w);
    choice_w{t}(:,2) = choice_w{t}(:,1) + savings;
    choice_w{t}(:,3) = log(choice_w{t}(:,1)) - delta + ...
      beta * interp1(policy_w{t+1}(:,1), policy_w{t+1}(:,3), wealth_grid_w, 'linear', "extrap");

    % If the worker chooses to retire in period t:
    choice_w{t}(:,4) = policy_r{t}(:,2);
    choice_w{t}(:,5) = policy_r{t}(:,1);
    choice_w{t}(:,6) = policy_r{t}(:,3); 

%     % Add credit constraint
    choice_w{t} = [NaN(1, size(choice_w{t}, 2)); choice_w{t}];
    choice_w{t}(1,[1 2 4 5]) = 0;

    % Step 2: Refine the endogenous wealth grid for the worker's problem when choosing to keep working.
    temp_wealth = choice_w{t}(:,2); temp_value = choice_w{t}(:,3); temp_consum = choice_w{t}(:,1);
    [region_logical, region_index] = fun_detect(temp_wealth);
    ori_value_cond = [temp_wealth(~isnan(temp_value)), temp_value(~isnan(temp_value))];
    ori_consum_cond = [temp_wealth, temp_consum];
    it = 0;
    while true
      if ~region_logical % refinement completed or no need for refinement
        refined_value_cond = ori_value_cond;
        refined_consum_cond = ori_consum_cond;
        break

      elseif region_logical % refinement process
        ori_value_cond = fun_UpperEnvelope2(ori_value_cond, InterpMethod, 100);
        ori_consum_cond = fun_UpperEnvelope2(ori_consum_cond, InterpMethod, 100);
        temp_wealth = ori_value_cond(:,1);

%         ori_value_cond = [temp_wealth, temp_value];
%         [refined_value_cond, index_remaining] = fun_UpperEnvelope(ori_value_cond, InterpMethod, gridN);
%         temp_wealth = refined_value_cond(:,1);
%         temp_value = refined_value_cond(:,2);
%         temp_consum = temp_consum(index_remaining);

%         ori_value_cond = [temp_wealth, temp_value];
%         [refined_value_cond, index_remaining, kink_index, kink_wealth, kink_value] = fun_UpperEnvelope(ori_value_cond, InterpMethod, gridN, true);
%         kink_consum = interp1(refined_value_cond(1:region_index,1), temp_consum(1:region_index), kink_wealth, InterpMethod, 'extrap');
%         refined_value_cond_inserted = fun_insert(refined_value_cond, kink_index, kink_wealth, kink_value);
%         refined_consum_inserted =  fun_insert([refined_value_cond(:,1), temp_consum(index_remaining)], kink_index, kink_wealth, kink_consum);
%         temp_wealth = refined_value_cond_inserted(:,1);
%         temp_consum = refined_consum_inserted(:,2);
%         temp_value = refined_value_cond_inserted(:,2);

        [region_logical, region_index] = fun_detect(temp_wealth); 
        it = it+1;
      end % end the if statement to check whether we need to further refine the endogenous wealth grid
    end % end the refinement of the endogenous wealth grid
    fprintf("iteration to refinements: %i \n", it)
    % Step 3: Choose d_t to maximize the choice-specific value function.
    %         And then store the maximized choice-specific mappings as the unconditional mappings.
%     [temp_wealth_r, index_unique] = unique(choice_w{t}(:,5));
%     temp_valur_r = choice_w{t}(:,6); temp_valur_r = temp_valur_r(index_unique);
%     temp_consum_r = choice_w{t}(:,4); temp_consum_r = temp_consum_r(index_unique);
%     [wealth, value, index] = fun_Choice(temp_wealth, temp_value, temp_wealth_r, temp_valur_r, [wealth_min, wealth_max], 'linear', gridN);

    [value, index] = fun_Choice2(refined_value_cond, [choice_w{t}(:,5), choice_w{t}(:,6)] , [wealth_min, wealth_max], InterpMethod, gridN);
    consum_w = interp1(refined_consum_cond(:,1), refined_consum_cond(:,2), value(:,1), 'linear', 'extrap');
    consum_r = interp1(choice_w{t}(:,5), choice_w{t}(:,4), value(:,1), 'linear', 'extrap');
    consum = [consum_w(index==1); consum_r(index==2)];
    policy_w{t}(:,1) = value(:,1);
    policy_w{t}(:,2) = consum;
    policy_w{t}(:,3) = value(:,2);
  end % end the if statement to check if the current period less than T
end % end the for loop that goes through t backward

%%
hold off
figure1 = figure(1);
ax = gca;
hold on
a = policy_w{18}(:,1); b =  policy_w{18}(:,2);
a = a(b<=40); b = b(b<=40);
plot(a, b, 'Color', 'blue')
plot(policy_w{10}(:,1), policy_w{10}(:,2), 'Color', 'red')
plot(policy_w{1}(:,1), policy_w{1}(:,2), 'Color', 'black')
grid on
ax.YTick = 0:5:40;
