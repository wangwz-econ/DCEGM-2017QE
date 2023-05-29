clear all
load test_data_forFig1
%% GRAPH ILLUSTRATION
a = wealth; b = valuec;
InterpMethod = 'line'; option = optimset('TolX', 1e-5); gridN=2000;
f_forzero = @(x) interp1(a1, b1, x, InterpMethod, 'extrap') - interp1(a2, b2, x, InterpMethod, 'extrap');
kink = fzero(f_forzero, 10, option); % 30.5922
hold off

figure1 = figure(1);
hold on
plot(a1, b1, 'Color', 'red', 'Marker', "x", 'MarkerEdgeColor', 'black')
plot(a2, b2, 'Color', 'blue', 'Marker', "o", 'MarkerEdgeColor', 'black', 'MarkerSize', 4)
xline(kink)
xline([a(12) a(11)], "--", {num2str(a(12)), num2str(a(11))},'LabelOrientation', 'horizontal')
for i=1:11
  if i<8
    text(a1(i), b1(i), num2str(i), 'FontSize', 8, 'VerticalAlignment', 'bottom', 'Color', 'red')
  else
    text(a1(i), b1(i), num2str(i), 'FontSize', 8, 'VerticalAlignment', 'top', 'Color', 'red')
  end
end
for i=12:20
  text(a2(i-11), b2(i-11), num2str(i), 'FontSize', 8, 'VerticalAlignment', 'top', 'Color', 'blue')
end
title("Figure 1 - Panel d")
subtitle("Solutions to Euler equation are not necessarily solutions to the Bellman equation", "FontSize", 8)
xlabel("Wealth")
ylabel("Choice-Specific Value Function when Choosing to Work")
saveas(figure1, "FIG1-PanelD _ original.jpg")

%% FUNCTION TESTING 1 (no kink interpolation)

% value_cond = [a, b];
% [refined_value_cond, ~] = test_fun_UpperEnvelope(value_cond, 'linear', 2000);
% 
% figure2 = figure(2);
% hold on
% plot(a1, b1, 'Color', 'red', 'Marker', "x", 'MarkerEdgeColor', 'black')
% plot(a2, b2, 'Color', 'blue', 'Marker', "o", 'MarkerEdgeColor', 'black', 'MarkerSize', 4)
% plot(refined_value_cond(1:13,1), refined_value_cond(1:13,2), ...
%   'LineStyle', "-.", 'Color', "#7E2F8E", 'Marker', "*", 'MarkerEdgeColor', 'black', 'MarkerSize', 4)
% title("Figure 1 - Panel d")
% subtitle(["Construct the upper envelope of the choice-specific value function","Without kink insertion"], "FontSize", 8)
% xlabel("Wealth")
% ylabel("Choice-Specific Value Function when Choosing to Work")

%% FUNCTION TESTING 2 (with kink interpolation)

value_cond = [a, b];
[refined_value_cond, index, kink_index, kink_wealth, kink_value] = test_fun_UpperEnvelope(value_cond, 'linear', 2000, true);

test = test_fun_insert(refined_value_cond(1:13,:), kink_index, kink_wealth, kink_value);

figure3 = figure(3);
hold on
plot(a1, b1, 'Color', 'red', 'Marker', "x", 'MarkerEdgeColor', 'black')
plot(a2, b2, 'Color', 'blue', 'Marker', "o", 'MarkerEdgeColor', 'black', 'MarkerSize', 4)
plot(test(1:13,1), test(1:13,2), ...
  'LineStyle', "-.", 'Color', "#7E2F8E", 'Marker', "*", 'MarkerEdgeColor', 'black', 'MarkerSize', 4)
title("Figure 1 - Panel d")
subtitle(["Construct the upper envelope of the choice-specific value function","With kink insertion"], "FontSize", 8)
xlabel("Wealth")
ylabel("Choice-Specific Value Function when Choosing to Work")
saveas(figure3, "FIG1-PanelD _ after refinement with kink insertion.jpg")

%% FUNCTION TESTING
% value_cond = [a, b];
% refined_value_cond = test_fun_UpperEnvelope2(value_cond, 'linear', 100);
% 
% figure4 = figure(3);
% hold on
% plot(a1, b1, 'Color', 'red', 'Marker', "x", 'MarkerEdgeColor', 'black')
% plot(a2, b2, 'Color', 'blue', 'Marker', "o", 'MarkerEdgeColor', 'black', 'MarkerSize', 4)
% c = refined_value_cond(:,1); d = refined_value_cond(:,2);
% c = c(c<=40); d = d(c<=40);
% plot(c, d, ...
%    'Color', "#7E2F8E", 'Marker', "*", 'MarkerEdgeColor', 'black', 'MarkerSize', 4)
% title("Figure 1 - Panel d")
% subtitle(["Construct the upper envelope of the choice-specific value function","Without kink insertion"], "FontSize", 8)
% saveas(figure3, "FIG1-PanelD _ 4 _ test fun_UpperEnvelope2.jpg")


%% FUNCTION TESTING
% [refined_value_cond, index_remaining] = fun_UpperEnvelope(ori_value_cond, InterpMethod, 2000);
% 
% wealth = ori_value_cond(:, 1); value = ori_value_cond(:, 2);
% 
% [~, start] = fun_detect(wealth);
% for step = 1:(length(wealth) - start - 1)
%   increasing = (wealth(step+start) > wealth(step+start-1));
%   if increasing, break;  end
% end
% 
% wealth1 = wealth(1:start-1); 
% wealth3 = wealth((start+step):end); 
% 
% line1 = [wealth(1:start-1), value(1:start-1)];
% line2 = [wealth((start-1):(start+step-1)), value((start-1):(start+step-1))];
% line3 = [wealth((start+step):end), value((start+step):end)];
% 
% overlap_min = wealth(start); overlap_max = wealth(start-1);
% overlap = linspace(overlap_min, overlap_max, gridN)';
% 
% extr_1 = interp1(line1(:,1), line1(:,2), overlap, InterpMethod, 'extrap');
% extr_2 = interp1(line2(:,1), line2(:,2), overlap, InterpMethod, 'extrap');
% extr_3 = interp1(line3(:,1), line3(:,2), overlap, InterpMethod, 'extrap');
% [~, max_index] = max([extr_1, extr_2, extr_3], [], 2, 'includenan');
% 
% wealth_threshold_index = find(max_index==3, 1, 'first');

%% DEVELOPEMENT OF THE FUNCTION
% [~, index_start] = fun_detect(a);
% for step = 1:(length(a)-index_start - 1)
%   increasing = (a(step+index_start) > a(step+index_start-1));
%   if increasing, break;  end
% end
% 
% wealth1 = a(1:index_start-1); value1 = b(1:index_start-1);
% wealth3 = a((index_start+step):end); value3 =  b((index_start+step):end);
% 
% line1 = [a(1:index_start-1), b(1:index_start-1)];
% line2 = [a((index_start-1):(index_start+step-1)), b((index_start-1):(index_start+step-1))];
% line3 = [a((index_start+step):end), b((index_start+step):end)];
% 
% overlap_min = a(index_start); overlap_max = a(index_start-1);
% overlap = linspace(overlap_min, overlap_max, gridN)';
% 
% extr_1 = interp1(line1(:,1), line1(:,2), overlap, InterpMethod, 'extrap');
% extr_2 = interp1(line2(:,1), line2(:,2), overlap, InterpMethod, 'extrap');
% extr_3 = interp1(line3(:,1), line3(:,2), overlap, InterpMethod, 'extrap');
% [value, index] = max([extr_1, extr_2, extr_3], [], 2, 'includenan');
% wealth = [overlap(index==1); overlap(index==2); overlap(index==3)];
% value = [line1(:,2); value; line3(:,2)];
% [wealth, sort_index] = sort(wealth);
% value = value(sort_index);
% 
% [~, max_index] = max([extr_1, extr_2, extr_3], [], 2, 'includenan');
% wealth_threshold_index = find(max_index==3, 1, 'first');
% wealth_threshold = overlap(wealth_threshold_index);
% 
% index_remaining = [find(wealth1<=wealth_threshold); find(wealth3>wealth_threshold)+index_start-step+1];
% 
% res_line = [a(index_remaining), b(index_remaining)];
