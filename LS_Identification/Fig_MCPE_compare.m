% Clear environment
clear; clc; close all;
ERROR_MCPE_all = [];
S_min = [];
T_i = [5; 6; 7; 10; 11; 13; 14; 14; 14; 14];
m = 2;
L = 5;
N = 10;

for sv = 0.1:0.1:30
    process = [0.005, 40, sv, 1];
    [ERROR_MCPE, u_MCPE_all] = function_MCPE_compare(process);
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
T = 0.1:0.1:30;  % Horizontal axis example

% Assume we have data from three algorithms (e.g., several repeated experiments)
% Here, for demonstration, we randomly simulate: each algorithm does k experiments for each T

data_algo1 = ERROR_MCPE_all;  % Data from the first algorithm

% Calculate the mean and standard deviation for each algorithm at each T
mean_algo1 = mean(data_algo1, 2);
std_algo1 = std(data_algo1, 0, 2);

ERROR_MCPE_all = [];
S_min = [];
for sv = 0.1:0.1:30
    process = [0.005, 40, 1, sv];
    [ERROR_MCPE, u_MCPE_all] = function_MCPE_compare(process);
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

% Calculate the mean and standard deviation for each algorithm at each T
mean_algo2 = mean(data_algo2, 2);
std_algo2 = std(data_algo2, 0, 2);

ERROR_MCPE_all = [];
S_min = [];
for sv = 1
    process = [0.005, 40, 1, sv];
    [ERROR_MCPE, u_MCPE_all] = function_MCPE_compare(process);
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
plot(T, mean_algo1, '-', ...   % Mean curve
    'Color', 'r', ...
    'LineWidth', 2, ...
    'MarkerFaceColor','r');
plot(T(1:9:300), mean_algo1(1:9:300), 'o', ...   % Mean curve
    'Color', 'r', ...
    'LineWidth', 2, ...
    'MarkerFaceColor','r');

upper2 = mean_algo2 + std_algo2;
lower2 = mean_algo2 - std_algo2;
x2_fill_MCPE = [T, fliplr(T)];
y2_fill_MCPE = [upper2', fliplr(lower2')];

fill(x2_fill_MCPE, y2_fill_MCPE, 'b', ...
    'FaceAlpha', 0.2, ...
    'EdgeColor', 'none');
plot(T, mean_algo2, '-', ...
    'Color', 'b', ...
    'LineWidth', 2, ...
    'MarkerFaceColor','b');

plot(T(1:9:300), mean_algo2(1:9:300), 'o', ...   % Mean curve
    'Color', 'b', ...
    'LineWidth', 2, ...
    'MarkerFaceColor','b');

save('x1_fill_MCPE_compare.mat','T');  % Save variable T to a file
save('x2_fill_MCPE_compare.mat','T');  % Save variable T to a file
% save('x3_fill_MCPE.mat','T');  % Save variable T to a file
save('y1_fill_MCPE_compare.mat','data_algo1');  % Save variable data_algo1 to a file
save('y2_fill_MCPE_compare.mat','data_algo2');  % Save variable data_algo2 to a file
% save('y3_fill_MCPE.mat','data_algo3');  % Save variable data_algo3 to a file
