% Clear environment
clear; clc; close all;
ERROR_MCPE_all = [];
S_min = [];
rng(0);
T_i = sort(randi([5, 13], 10, 1));

m = 2;
L = 5;
N = 10;

for sv = 1:10:101
    process = [0.5, sv];
    [ERROR_MCPE, u_MCPE_all] = function_MCPE(process);
    ERROR_MCPE_all = [ERROR_MCPE_all; ERROR_MCPE];
    for i = 1:length(u_MCPE_all)
        U_all{i} = reshape(u_MCPE_all{i}, T_i(i) * m, 1);
        H{i} = hankel_r(U_all{i}, L, T_i(i) - L + 1, m);
    end

    H_MCPE = [];

    for i = 1:N
        H_MCPE = [H_MCPE H{i}];
    end

    S = svd(H_MCPE);
    S_min = [S_min min(S)];
end

% Assume we have horizontal axis T, with length m
T = S_min;  % Horizontal axis example

% Assume we have data from three algorithms (e.g., several repeated experiments)
% Here, for demonstration, we randomly simulate: each algorithm does k experiments for each T

data_algo1 = ERROR_MCPE_all;  % Data from the first algorithm

ERROR_MCPE_all = [];
S_min = [];
for sv = 1:10:101
    process = [0.05, sv];
    [ERROR_MCPE, u_MCPE_all] = function_MCPE(process);
    ERROR_MCPE_all = [ERROR_MCPE_all; ERROR_MCPE];
    for i = 1:length(u_MCPE_all)
        U_all{i} = reshape(u_MCPE_all{i}, T_i(i) * m, 1);
        H{i} = hankel_r(U_all{i}, L, T_i(i) - L + 1, m);
    end

    H_MCPE = [];

    for i = 1:N
        H_MCPE = [H_MCPE H{i}];
    end

    S = svd(H_MCPE);
    S_min = [S_min min(S)];
end

data_algo2 = ERROR_MCPE_all;  % Data from the second algorithm

ERROR_MCPE_all = [];
S_min = [];
for sv = 1:10:101
    process = [0.005, sv];
    [ERROR_MCPE, u_MCPE_all] = function_MCPE(process);
    ERROR_MCPE_all = [ERROR_MCPE_all; ERROR_MCPE];
    for i = 1:length(u_MCPE_all)
        U_all{i} = reshape(u_MCPE_all{i}, T_i(i) * m, 1);
        H{i} = hankel_r(U_all{i}, L, T_i(i) - L + 1, m);
    end

    H_MCPE = [];

    for i = 1:N
        H_MCPE = [H_MCPE H{i}];
    end

    S = svd(H_MCPE);
    S_min = [S_min min(S)];
end

data_algo3 = ERROR_MCPE_all;  % Data from the third algorithm

% Calculate the mean and standard deviation for each algorithm at each T
mean_algo1 = mean(data_algo1, 2);
std_algo1 = std(data_algo1, 0, 2);

mean_algo2 = mean(data_algo2, 2);
std_algo2 = std(data_algo2, 0, 2);

mean_algo3 = mean(data_algo3, 2);
std_algo3 = std(data_algo3, 0, 2);

%===== Plotting section =====
figure; hold on; box on;

% (1) Shaded region can be drawn using fill:
%     Upper and lower bounds = mean ± standard deviation
upper1 = mean_algo1 + std_algo1;
lower1 = mean_algo1 - std_algo1;
x1_fill_MCPE = [T, fliplr(T)];  % For filling, concatenate
y1_fill_MCPE = [upper1', fliplr(lower1')];

fill(x1_fill_MCPE, y1_fill_MCPE, 'r', ...   % Shaded color: red
    'FaceAlpha', 0.2, ...       % Transparency of shaded area
    'EdgeColor', 'none');       % No edge line
plot(T, mean_algo1, '-o', ...   % Mean curve
    'Color', 'r', ...
    'LineWidth', 2, ...
    'MarkerFaceColor','r');
% (2) Algorithm 2
upper2 = mean_algo2 + std_algo2;
lower2 = mean_algo2 - std_algo2;
x2_fill_MCPE = [T, fliplr(T)];
y2_fill_MCPE = [upper2', fliplr(lower2')];

fill(x2_fill_MCPE, y2_fill_MCPE, 'b', ...
    'FaceAlpha', 0.2, ...
    'EdgeColor', 'none');
plot(T, mean_algo2, '-d', ...
    'Color', 'b', ...
    'LineWidth', 2, ...
    'MarkerFaceColor','b');

% (3) Algorithm 3
upper3 = mean_algo3 + std_algo3;
lower3 = mean_algo3 - std_algo3;
x3_fill_MCPE = [T, fliplr(T)];
y3_fill_MCPE = [upper3', fliplr(lower3')];

fill(x3_fill_MCPE, y3_fill_MCPE, 'g', ...
    'FaceAlpha', 0.2, ...
    'EdgeColor', 'none');
plot(T, mean_algo3, '-s', ...
    'Color', 'g', ...
    'LineWidth', 2, ...
    'MarkerFaceColor','g');

set(gca, 'YScale', 'log');          % Set the y-axis to logarithmic scale

ylim([1e-3, 1e-1]);                 % Set y-axis limits
yticks([1e-3, 1e-2, 1e-1]);         % Select desired ticks
yticklabels({'10^{-3}','10^{-2}','10^{-1}'});  % Set tick labels
% Add legend, labels, title, etc.
legend({'Algorithm1','Algorithm2','Algorithm3'}, 'Location', 'best');
xlabel('Length of each experiment T');
ylabel('Relative error');
title('Comparison of Different Methods');

% save('x1_fill_MCPE.mat','T');  % Save variable T to a file
% save('x2_fill_MCPE.mat','T');  % Save variable T to a file
% save('x3_fill_MCPE.mat','T');  % Save variable T to a file
% save('y1_fill_MCPE.mat','data_algo1');  % Save variable data_algo1 to a file
% save('y2_fill_MCPE.mat','data_algo2');  % Save variable data_algo2 to a file
% save('y3_fill_MCPE.mat','data_algo3');  % Save variable data_algo3 to a file
