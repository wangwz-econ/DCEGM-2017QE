% A Retirement Model with Consumption and Savings, Income Shocks and Credit Constraint
% The Originial Paper is "The Endogenous Grid Method for Discrete-Continuous Dynamic Choice Models 
% with (or without) Taste Shocks" (QE, 2017)

classdef model2_DCEGM_retirement < handle

properties (Access = public)
  label = "Consumption Model with Retirement";
  T = 25;           % time period
  gridN = 500;      % number of grid points
  wealth_max = 50;  % maximum wealth (end point of the gird)
  quadN = 5;        % # of quadrature points to calculate the integral
  simN = 10;        % number of trajectories about to simulate
  init = [10, 30];  % interval of initial wealth distribution
  r = 0.05;         % return to assets holding
  beta = 0.95;      % discount factor
  sigma = 0.25;     % sigma parameter in income shocks
  duw = 0.35;       % disutility of working
  theta = 1.95;     % CRRA of the utility function
  inc0 = 0.75;      % income equation: cons
  inc1 = 0.04;      % income equation: coef on age
  inc2 = 0.0002;    % income equation: coef on age^2
  c_min = 0.001;    % consumption floor (safety net in retirement)
  lambda = 0.02;    % scale of extreme value distributed taste shocks
end

properties (SetAccess=private, GetAccess=public)
  policy = polyline;  % polyline is another class used to conduct interpolation related methods
  value = polyline;
  sims = struct();    % store simulation trajectories
end

% working = 1 for working, 2 for retiring

methods (Access=public)
  % Definition of the model

  function model = model2_DCEGM_retirement(label,T,gridN,wealth_max,quadN,simN,init,r,beta,sigma,...
                                          duw,theta,inc0,inc1,inc2,c_min,lambda) 
    model.label = label;
    model.T = T;
    model.gridN = gridN;
    model.wealth_max = wealth_max;
    model.quadN = quadN;
    model.simN =simN;
    model.init = init;
    model.r = r;
    model.beta = beta;
    model.sigma = sigma;
    model.duw = duw;
    model.theta = theta;
    model.inc0 = inc0;
    model.inc1 = inc1;
    model.inc2 = inc2;
    model.c_min = c_min;
    model.lambda = lambda;
  end


  function u = utility(model, consumption, working) % utility function
    if model.theta==1
      u = log(consumption);
    else
      u = (consumption.^(1-model.theta)-1)/(1-model.theta);
    end
    u = u - model.duw * (working==1);
  end

  function mu = u_D(model, consumption)   % marginal utility function (u')
    mu = consumption .^ (-model.theta);
  end

  function cons = u_D_inv(model, mutil)   % inverse of the marginal utility function
    cons = mutil .^ (-1/model.theta);
  end

  function w = income(model, t, shock)    % income in period t with given normal shock
    % assume t=1 is age=20
    age = t + 19;
    w = exp(model.inc0 + model.inc1 * age - model.inc2 * (age.^2) + shock);
  end

  function w1 = budget(model, t, savings, shocks, working)   
    % wealth in period t+1 (under different realiztion of income shocks)
    % inputs: savings = 1x(gridN) row vector of savings
    %				  shocks = (quadN)x1 column vector of shocks
    % output: wealth in period t+1, (quadN by gridN) matrix of all possible next period wealth
    w1 = ones(size(shocks,1), 1) * savings(1,:) * (1 + model.r) + ...
         (working==1) * income(model, t+1, shocks(:,1) * ones(1, size(savings, 2)));
  end

  function mw1 = mbudget(model, t, savings, shocks, working) %% not so sure what this is?
    % derivatives of wealth in t+1 w.r.t. savings
    % inputs and outputs as above
    mw1 = repmat(1+model.r, size(shocks, 1), size(savings, 2));
  end

  function solve_dcegm(model)   % solve the model with DC-EGM algorithm
    if model.lambda < eps && model.sigma > eps
      error(sprintf(['Solution for the model with income shocks but without taste shocks is not supported!\n' ...
        'Can not use quadrature to calculate expectations of discontinuous functions.\n' ...
        'Set lambda > 0 or sigma = 0']))
    end
    model.policy = polyline;  % initialize previous optimal consumption function 
    model.value  = polyline;  % initialize previous value function

    % step 1: first solve the retiree problem

    % step 1.1: quadrature points and exogenous savings gridpoints
    % set up the quadrature points to calculte the integral (expectation over all possible
    % income shocks) and also set up the exogenous savings grid

    [quadp, quadw] = quadpoints(model.quadN, 0, 1);
    quadstnorm = norminv(quadp, 0, 1);
    savingsgrid = linspace(0, model.wealth_max, model.gridN);

    % step 1.2: main EGM loop

    for t = model.T:-1:1
      fprintf('t: %d \n', t)
      if t == model.T % in the terminal period T

        for id = 1:2
          model.policy(id, t) = polyline([0, model.wealth_max], [0, model.wealth_max], 'Policy function in period T');
          model.value(id, t)  = polyline([0, model.wealth_max], [0, NaN], 'Value function in period T');
          %vf(1)=0 is important, otherwise vf cannot be computed in the terminal period
        end

      else % in periods other than the terminal period
        for id=1:2 % for each decision, 1=remain working, 2=retiring
          wk1 = model.budget(t, savingsgrid, quadstnorm*model.sigma, id);
          wk1 = max(model.c_min, wk1); % to ensure minimum consumption and to avoid NaNs
          vl1 = model.value_function(1, t+1, wk(1:end)); % next period value of working, reshaped to row vector
          vl1(2,:) = model.value_function(2, t+1, wk(1:end));% next period value of retiring, row vector
          pr1 = model.chpr(vl1)*(id==1);%probability to remain working tomorrow (conditional on being a worker)
          cons1 = model.policy
        end
      end
    end


  end


end



end



