% A retirement model with consumption and savings, 
% income shocks and credit constraint

classdef model2_DCEGM_retirement < handle

properties (Access = public)
    label = "Consumption Model with Retirement";
    T = 25;
    gridN = 500;
    wealth_max = 50;
    quadN = 5;
    simN = 10;
    init = [10, 30];
    r = 0.05;
    beta = 0.95;
    sigma = 0.25;
    duw = 0.35;
    theta = 1.95;
    inc0 = 0.75;
    inc1 = 0.04;
    inc2 = 0.0002;
    c_min = 0.001;
    lambda = 0.02;
end

properties (SetAccess=private, GetAccess=public)
    policy = polyline;
    value = polyline;
    sims = struct()
end

% working = 1 for working, 2 for retiring

methods (Access=public)
    % Definition of the model
    % Vectorize!!
    
    function u = util(me, consumption, working)
        if me.theta==1
            u = log(consumption);
        else
            u = (consumption.^(1-me.theta)-1)/(1-me.theta);
        end 
        u = u - me.duw * (working==1);
    end

    function mu = mutil(me, consumption)
        if me.theta == 1
            mu = 1 ./ consumption;
        else
            mu = consumption .^ (-me.theta);
        end
    end

    function cons = imutil(me, mutil)
        if me.theta==1
            cons = 1 ./ mutil;
        else
            cons = mutil .^ (-1/me.theta);
        end
    end

    function w = income(me, it, shock)
        % income in period it with given normal shock
        % assume it=1 is age=20
        age = it + 19;
        w = exp(me.inc0 + me.inc1*age - me.inc2*age .*age + shock);
    end

    function w1 = budget(me, it, savings, shocks, working)
        %wealth in period t+1, where it=t
		    %inputs: savings = 1x(gridN) row vector of savings
		    %				 shocks = (quadN)x1 column vector of shocks
		    %output: 
        % w1 =quadN by gridN matrix of all possible next period wealth
        w1 = ones(size(shocks,1),1) * savings(1,:) * (1+me.r) + ...
            (working==1) *income(me, it+1, shocks(:,1)*ones(1, size(savings, 2)));
    end

    function mw1 = mbudget(me, it, savings, shocks, working)
        %derivatives of wealth in t+1 wrt savings
        %inputs and outputs as above
        mw1 = repmat(1+me.r, size(shocks, 1), size(savings, 2));
    end

    function solve_dcegm(me)
        %solve the model with DC-EGM algorithm
        if me.lambda < eps && me.sigma > eps
            error(sprintf(['Solution for the model with income shocks but with without taste shocks is not supported!\n' ...
                'Can not use quadrature to calculate expectations of discontinuous functions.\n' ...
                'Set lambda > 0 or sigma = 0']))
        end
        me.policy = polyline; % delete previous solution
        me.value = polyline; % delete previous solution

        % first solve the retiree problem

        [quadp, quadw] = quadpoints(me.quadN, 0, 1);
        quadstnorm = norminv(quadp, 0, 1);
        savingsgrid = linspace(0, me.wealth_max, me.gridN);

        % main EGM loop

        fprintf('t:')
        for it = me.T:-1:1
            fprintf(' %d', it)
            if it == me.T
                %terminal period
                for id = 1:2
                    me.policy(id,it)=polyline([0, me.wealth_max], [0, me.wealth_max], 'Policy function in period T');
                    me.value(id,it)=polyline([0, me.wealth_max], [0, NaN], 'Value function in period T');
                    %vf(1)=0 is important, otherwise vf cannot be computed
                    %in terminal period
                end
            else
                %not the terminal period
                for id=1:2 %for each decision, 1=remain working, 2=retiring
                    wk1 = me.budget(it, savingsgrid, quadstnorm*me.sigma, id);
                    wk1 = max(me.c_min wk1); %to insure minimum consumption and to avoid NaNs
                    vl1 = me.value_function(1, it+1, wk(1:end)); %next period value of working, reshapted to row vector
                    vl1(2,:) = me.value_function(2, it+1, wk(1:end));%next period value of retiring, row vector
                    pr1 = me.chpr(vl1)*(id==1);%probability to remain working tomorrow (conditional on being a worker)
                    cons1=me.policy
                end
            end
        end
  

    end


end



end



