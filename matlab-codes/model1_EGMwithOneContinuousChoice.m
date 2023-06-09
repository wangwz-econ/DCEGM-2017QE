% EGM with only one continuous state variable
% model is charaterized by constant savings rate and log normal
% distribution of wage income

clear all

% step 1: set the parameters

quadN = 10;       % number of quadrature points to calculate expectation
wealth_max = 10;  % maximum wealth
gridN = 100;      % number of grid points

T = 25;           % number of time periods
R = 1.05;         % gross return to savings
beta = 0.95;      % discount factor
y = 1;            % wage income
sigma = 0.25;     % sigma parameter in lognormal distribution

% quadrature notes and weights
[quadp, quadw]=quadpoints(quadN, 0, 1);  
% prepare quadrature points for calculation of expectations of Normal
quadstnorm = norminv(quadp, 0, 1);         

grid_savings = linspace(0, wealth_max, gridN); % savings gridpoints

% For each period t, get the optimal consumption function,
% characterized by two row vectors in the struct policy 

policy{T}.w = [0 wealth_max];        % terminal period wealth
policy{T}.c = [0 wealth_max];        % terminal period optimal consumption

for t = (T-1):-1:1
    w1 = y + exp(quadstnorm * sigma) * R * grid_savings;
    c1 = interp1(policy{t+1}.w, policy{t+1}.c, w1, 'linear', 'extrap');
    rhs = quadw' * (1 ./ c1);
    policy{t}.c = [0, 1./(beta*R*rhs)];
    policy{t}.w = [0, grid_savings + policy{t}.c(2:end)];
end


% Plot the optimal policy functions
figure1 = figure(1);
for t = (T):-1:1
	plot(policy{t}.w, policy{t}.c)
	hold all
end
set(gca, 'XLim',[0 wealth_max])
xlabel(gca, 'Wealth')
ylabel(gca,' Optimal consumption')
title(gca, 'Optimal consumption rules by age')
saveas(figure1, "m1 _ policy function by different ages.jpg")
hold on


% wealth0 = 5;
% path_consum = NaN(T, 1);
% path_asset  = NaN(T+1, 1);
% path_asset(1,1) = wealth0;
% 
% for t = 1:1:T-1
%   path_consum(t,1)  = interp1(policy{t}.w, policy{t}.c, path_asset(t,1), 'spline', 'extrap');
%   path_asset(t+1,1) = interp1(policy{t}.w, R*(y+policy{t}.c), path_asset(t,1), 'spline', 'extrap');
% end
% figure2 = figure(2);
% plot(path_consum)
% xlabel('Time')
% ylabel('Consumption Path')
% title('Consumption Path with Initial Wealth 5')
% saveas(figure2, "m1 _ consumption path.jpg")
