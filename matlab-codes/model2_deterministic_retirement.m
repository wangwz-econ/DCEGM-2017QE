%% DESCRIPTION
% DC-EGM for the deterministic model
% The only difficulty is to refine endogenous wealth grid by first detecting the non-motonotic region
%   and second deleting certain points through the upper enevelope of choice-specific value function
%   in that region.
% We use the cls_outer_refinement method of the cls_line class to implement this step!

% The parameters are taken from the note of Figure 2 in the original paper.

%% Step 1: Set up Parameters
clear all

% Economic Parameters
beta = 0.98; R = 1; y = 20; T = 20; delta = 1;
u_D = @ (c) 1./c;     % marginal utility function
u_D_inv = @ (u) 1./u; % the inverse of marginal utility function

% Computational Parameters
InterpMethod = 'spline';
gridN = 1e+3;
wealth_min = 1e-5; wealth_max = 400;
savings = linspace(wealth_min, wealth_max, gridN)'; % exogenous savings grid (gridN by 1)
wealth_grid_w = R * savings + y;       % wealth grid if the agent chooses to work (gridN by 1)
wealth_grid_r = R * savings;           % wealth grid if the agent chooses to retire (gridN by 1)

% Variable policy_j, j = w, r is a T by 1 cell.
% In each cell t, it stores a matrix denoting the worker's or retiree's value and policy function.
%   In each matrix of a cell of policy_r and policy_w, the first 3 columns are:
%     the 1st column is the endogenous wealth grid,
%     the 2nd column is optimal consumption given the savings grid,
%     the 3rd column is unconditional value function.
%   policy_w has a 4th column:
%     the 4th column is one's retirement decision (1 means to keep working, 0 means to retire).
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
% Because a solution to the Euler equation is not necessary a solution that maximizes the Bellman equation!

for t = T:-1:1
  fprintf('t = %i\n', t)
  if t == T
    policy_r{t}(:,1) = savings;                  % the endogenous wealth grid
    policy_r{t}(:,2) = savings;                  % optimal consumption given the savings grid
    policy_r{t}(:,3) = log(policy_r{t}(:,1));    % unconditional value function

    policy_w{t}(:,1) = savings;                  % the endogenous wealth grid
    policy_w{t}(:,2) = savings;                  % optimal consumption given the savings grid
    policy_w{t}(:,3) = log(policy_w{t}(:,1));    % unconditional value function
    policy_w{t}(:,4) = 0;                        % always choose to retire in the terminal period T
  else
    rhs_r = beta * R * u_D(interp1(policy_r{t+1}(:,1), policy_r{t+1}(:,2), wealth_grid_r, InterpMethod, "extrap"));
    rhs_w = beta * R * u_D(interp1(policy_w{t+1}(:,1), policy_w{t+1}(:,2), wealth_grid_w, InterpMethod, "extrap"));

    % Case 1: The retiree's problem.
    policy_r{t}(:,2) = u_D_inv(rhs_r);
    policy_r{t}(:,1) = savings + policy_r{t}(:,2);
    policy_r{t}(:,3) = log(policy_r{t}(:,2)) + ...
      beta * interp1(policy_r{t+1}(:,1), policy_r{t+1}(:,3), wealth_grid_r, InterpMethod, "extrap");

    % Case 2: The worker's problem.
    % Step 1: Collect choice-specific optimal mappings in variable choice_w using Euler equations.

    % If the worker chooses to keep working in period t:
    choice_w{t}(:,1) = u_D_inv(rhs_w);
    choice_w{t}(:,2) = choice_w{t}(:,1) + savings;
    choice_w{t}(:,3) = log(choice_w{t}(:,1)) - delta + ...
      beta * interp1(policy_w{t+1}(:,1), policy_w{t+1}(:,3), wealth_grid_w, InterpMethod, "extrap");

    % If the worker chooses to retire in period t:
    choice_w{t}(:,4) = policy_r{t}(:,2);
    choice_w{t}(:,5) = policy_r{t}(:,1);
    choice_w{t}(:,6) = policy_r{t}(:,3);

    % Add credit constraint
    choice_w{t} = [NaN(1, size(choice_w{t}, 2)); choice_w{t}];
    choice_w{t}(1,[1 2 4 5]) = 0;

    % Step 2: Refine the endogenous wealth grid for the worker's problem when choosing to keep working.
    value_w = cls_line(choice_w{t}(2:end,2), choice_w{t}(2:end,3)); % choice-specific value function when d_t=1
    consu_w = cls_line(choice_w{t}(:,2), choice_w{t}(:,1)); % choice-specific optimal consumption when d_t=1

    [value_w_refined, indx_del, new_dots] = cls_outer_refinement(value_w);
    consu_w_refined = cls_thinout(consu_w, indx_del+1); % delete the suboptimal solutions
    new_wealth = new_dots.x;
    new_consum = cls_interp(consu_w_refined, new_wealth);
    consu_w_refined = cls_grow(consu_w_refined, cls_line(new_wealth, new_consum), false);
    consu_w_refined = cls_sort(consu_w_refined);

    % Step 3: Compare the choice-specific value function across d_t=0, d_t=1
    value_w_refined_mat = [value_w_refined.x, value_w_refined.y];
    value_r_mat = [choice_w{t}(:,5), choice_w{t}(:,6)];
    [value, index] = fun_Choice2(value_w_refined_mat, value_r_mat, [wealth_min, wealth_max], InterpMethod, gridN);
    consum_w = transpose( cls_interp(consu_w_refined, value(:,1)) );
    consum_r = interp1(choice_w{t}(:,5), choice_w{t}(:,4), value(:,1), InterpMethod, 'extrap');
    consum = [consum_w(index==1); consum_r(index==2)];

    % Step 4: Store the uncondition value function and optimal consumption function using policy_w cell
    policy_w{t}(:,1) = value(:,1);
    policy_w{t}(:,2) = consum;
    policy_w{t}(:,3) = value(:,2);
    index(index==2) = 0;
    index = reshape(index, [], 1);
    policy_w{t}(:,4) = index;
  end % end the if statement to check if the current period less than T
end % end the for loop that goes through t backward

%% Step 3: Draw the Optimal Consumption Functions in Different Periods
hold off
figure1 = figure(1);
ax = gca;
hold on
grid on
ax.YTick = 0:5:40;
w18 = policy_w{18}(:,1); c18 =  policy_w{18}(:,2);
w18 = w18(c18<=40); c18 = c18(c18<=40);
w10 = policy_w{10}(:,1); c10 = policy_w{10}(:,2);
w1 = policy_w{1}(:,1); c1 = policy_w{1}(:,2);
plot(w18, c18, 'Color', 'blue')
text(w18(end), c18(end), 't=18', 'VerticalAlignment', 'bottom', 'Color', 'blue')
plot(w10, c10, 'Color', 'red')
text(w10(end), c10(end), 't=10', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'red')
plot(w1, c1, 'Color', 'black')
text(w1(end), c1(end), 't=1', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'black')
xlabel("Wealth")
ylabel("Optimal Consumption")
title("Figure 2 - Panel a")
% title("Figure 3 - Panel a")
subtitle("Optimal consumption rules in period t = 18, 10, 1")
saveas(figure1, "m2 _ policy function for workers.jpg")
saveas(figure1, "FIG2-PanelA _ policy function for workers.jpg")
% saveas(figure1, "FIG3-PanelA _ policy function for workers.jpg")
hold off

%% Step 4: Simulate a Lifecycle Path Given an Initial Wealth Level for R = 1/beta
% How to smartly store the retirement decision?
if R > 1 
  figure2 = figure(2);
  a0_grid = [1, 50, 100, 150];
  for i = 1: 4
    a0 = a0_grid(i);
    consum_path = NaN(T, 1);   % consumption
    wealth_path = NaN(T+1, 1);   % wealth level
    workin_path = NaN(T, 1);   % working status

    wealth_path(1, 1) = a0;
    for t = 1:1:T                     % loop through time periods for a particular individual
      consum_path(t, 1) = interp1(policy_w{t}(:,1), policy_w{t}(:,2), wealth_path(t, 1), InterpMethod, 'extrap');
      retire_index = find(policy_w{t}(:,4)==0, 1, 'first');
      wealth_t = policy_w{t}(:,1);
      workin_path(t, 1) = (wealth_path(t,1) < wealth_t(retire_index));
      if workin_path(t, 1)<=0
        wealth_path(t+1, 1) = R * (wealth_path(t, 1) - consum_path(t, 1));
      else
        wealth_path(t+1, 1) = R * (wealth_path(t, 1) - consum_path(t, 1)) + y;
      end
    end
    hold on
    plot(2:T, consum_path(2:end), "LineWidth", 2)
    text(wealth_path(T), consum_path(T), sprintf("A0=%d", a0), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right')
    retire_period = find(workin_path==0, 1, 'first');
    plot(retire_period, consum_path(retire_period, 1), "O", "MarkerEdgeColor","black", "MarkerFaceColor", 'black')
    grid on
    xticks(0:2:20)
    yticks(19:0.1:20)
    xlabel("time period t")
    ylabel("simulated consumption path")
  end
  title("Figure 3 - Panel b")
  subtitle("Simulated Consumption Profile")
  saveas(figure2, "FIG3-PanelB _ simulated consumption profile.jpg")
end


