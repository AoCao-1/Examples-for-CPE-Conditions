% Clear environment
clear; clc; close all;
ERROR_HCPE_all = [];
S_min = [];
N = 10;
T_0 = 10;
T_i = [10; 10; 10; 7; 7; 8; 9; 10; 12; 12];
m = 2;
L = 5;
bar_p = 3;
alpha = ones(N, 1);

for sv = 1:10:101
    process = [0.5, sv];
    [ERROR_HCPE, u_HCPE_all] = function_HCPE(process);
    ERROR_HCPE_all = [ERROR_HCPE_all; ERROR_HCPE];

    T_max = T_0;  

    % Allocate and merge control input matrix, size m x T_max
    % The k-th column will represent the merged input at time step k
    u_sum = zeros(m, T_max);

    % Loop through each time step
    for k = 1:T_max
        % Perform weighted summation for all controllers
        total_u_k = zeros(m, 1);  % Merged quantity at time step k
        for i = 1:bar_p
            if k <= T_i(i)
                % If the i-th controller still has input at time k, accumulate
                total_u_k = total_u_k + alpha(i) * u_HCPE_all{i}(:, k);
            else
                % If k > T_i(i), this controller is considered to have no input
            end
        end
        % Store the result
        u_sum(:, k) = total_u_k;
    end

    U_all{1} = reshape(u_sum, T_0 * m, 1);
    H{1} = hankel_r(U_all{1}, L, T_0 - L + 1, m);
    for i = bar_p + 1:length(u_HCPE_all)
        U_all{i - bar_p + 1} = reshape(u_HCPE_all{i}, T_i(i) * m, 1);
        H{i - bar_p + 1} = hankel_r(U_all{i - bar_p + 1}, L, T_i(i) - L + 1, m);
    end

    H_HCPE = [];

    for i = 1:(N - bar_p + 1)
        H_HCPE = [H_HCPE H{i}];
    end

    S = svd(H_HCPE);
    S_min = [S_min min(S)];
end

% Assume we have a horizontal axis T, with length m
T = S_min;  % Horizontal axis example

% Assume data from three algorithms (e.g., results from repeated experiments)
% Here, we simulate randomly: each algorithm does k experiments for each T

data_algo1 = ERROR_HCPE_all;  % Data from the first algorithm

ERROR_HCPE_all = [];
S_min = [];
for sv = 1:10:101
    process = [0.05, sv];
    [ERROR_HCPE, u_HCPE_all] = function_HCPE(process);
    ERROR_HCPE_all = [ERROR_HCPE_all; ERROR_HCPE];
    T_max = T_0;  

    % Allocate and merge control input matrix, size m x T_max
    % The k-th column will represent the merged input at time step k
    u_sum = zeros(m, T_max);

    % Loop through each time step
    for k = 1:T_max
        % Perform weighted summation for all controllers
        total_u_k = zeros(m, 1);  % Merged quantity at time step k
        for i = 1:bar_p
            if k <= T_i(i)
                % If the i-th controller still has input at time k, accumulate
                total_u_k = total_u_k + alpha(i) * u_HCPE_all{i}(:, k);
            else
                % If k > T_i(i), this controller is considered to have no input
            end
        end
        % Store the result
        u_sum(:, k) = total_u_k;
    end

    U_all{1} = reshape(u_sum, T_0 * m, 1);
    H{1} = hankel_r(U_all{1}, L, T_0 - L + 1, m);
    for i = bar_p + 1:length(u_HCPE_all)
        U_all{i - bar_p + 1} = reshape(u_HCPE_all{i}, T_i(i) * m, 1);
        H{i - bar_p + 1} = hankel_r(U_all{i - bar_p + 1}, L, T_i(i) - L + 1, m);
    end

    H_HCPE = [];

    for i = 1:(N - bar_p + 1)
        H_HCPE = [H_HCPE H{i}];
    end

    S = svd(H_HCPE);
    S_min = [S_min min(S)];
end

data_algo2 = ERROR_HCPE_all;  % Data from the second algorithm

ERROR_HCPE_all = [];
S_min = [];
for sv = 1:10:101
    process = [0.005, sv];
    [ERROR_HCPE, u_HCPE_all] = function_HCPE(process);
    ERROR_HCPE_all = [ERROR_HCPE_all; ERROR_HCPE];
    T_max = T_0;  

    % Allocate and merge control input matrix, size m x T_max
    % The k-th column will represent the merged input at time step k
    u_sum = zeros(m, T_max);

    % Loop through each time step
    for k = 1:T_max
        % Perform weighted summation for all controllers
        total_u_k = zeros(m, 1);  % Merged quantity at time step k
        for i = 1:bar_p
            if k <= T_i(i)
                % If the i-th controller still has input at time k, accumulate
                total_u_k = total_u_k + alpha(i) * u_HCPE_all{i}(:, k);
            else
                % If k > T_i(i), this controller is considered to have no input
            end
        end
        % Store the result
        u_sum(:, k) = total_u_k;
    end

    U_all{1} = reshape(u_sum, T_0 * m, 1);
    H{1} = hankel_r(U_all{1}, L, T_0 - L + 1, m);
    for i = bar_p + 1:length(u_HCPE_all)
        U_all{i - bar_p + 1} = reshape(u_HCPE_all{i}, T_i(i) * m, 1);
        H{i - bar_p + 1} = hankel_r(U_all{i - bar_p + 1}, L, T_i(i) - L + 1, m);
    end

    H_HCPE = [];

    for i = 1:(N - bar_p + 1)
        H_HCPE = [H_HCPE H{i}];
    end

    S = svd(H_HCPE);
    S_min = [S_min min(S)];
end

data_algo3 = ERROR_HCPE_all;  

% Calculate the mean and standard deviation for each algorithm at each T
mean_algo1 = mean(data_algo1, 2);
std_algo1 = std(data_algo1, 0, 2);

mean_algo2 = mean(data_algo2, 2);
std_algo2 = std(data_algo2, 0, 2);

mean_algo3 = mean(data_algo3, 2);
std_algo3 = std(data_algo3, 0, 2);

%===== Plotting section =====
figure; hold on; box on;

% (1) Shaded area can be drawn using fill:
%     Upper and lower bounds = mean ± standard deviation
upper1 = mean_algo1 + std_algo1;
lower1 = mean_algo1 - std_algo1;
x1_fill_HCPE = [T, fliplr(T)];  % For fill, concatenate
y1_fill_HCPE = [upper1', fliplr(lower1')];

fill(x1_fill_HCPE, y1_fill_HCPE, 'r', ...   % Shaded color: red
    'FaceAlpha', 0.2, ...       % Transparency of shaded area
    'EdgeColor', 'none');       % No edge line
plot(T, mean_algo1, '-o', ...   % Mean curve
    'Color', 'r', ...
    'LineWidth', 2, ...
    'MarkerFaceColor','r');
% (2) Algorithm 2
upper2 = mean_algo2 + std_algo2;
lower2 = mean_algo2 - std_algo2;
x2_fill_HCPE = [T, fliplr(T)];
y2_fill_HCPE = [upper2', fliplr(lower2')];

fill(x2_fill_HCPE, y2_fill_HCPE, 'b', ...
    'FaceAlpha', 0.2, ...
    'EdgeColor', 'none');
plot(T, mean_algo2, '-d', ...
    'Color', 'b', ...
    'LineWidth', 2, ...
    'MarkerFaceColor','b');

% (3) Algorithm 3
upper3 = mean_algo3 + std_algo3;
lower3 = mean_algo3 - std_algo3;
x3_fill_HCPE = [T, fliplr(T)];
y3_fill_HCPE = [upper3', fliplr(lower3')];

fill(x3_fill_HCPE, y3_fill_HCPE, 'g', ...
    'FaceAlpha', 0.2, ...
    'EdgeColor', 'none');
plot(T, mean_algo3, '-s', ...
    'Color', 'g', ...
    'LineWidth', 2, ...
    'MarkerFaceColor','g');

set(gca, 'YScale', 'log');          % Set y-axis to logarithmic scale

ylim([1e-3, 1e-1]);                 % Set y-axis limits
yticks([1e-3, 1e-2, 1e-1]);         % Choose desired ticks
yticklabels({'10^{-3}','10^{-2}','10^{-1}'});  % Set tick labels
% Add legend, labels, and title
legend({'Algorithm1','Algorithm2','Algorithm3'}, 'Location', 'best');
xlabel('Length of each experiment T');
ylabel('Relative error');
title('Comparison of Different Methods');

% Save variables to .mat files
save('x1_fill_HCPE.mat', 'T');
save('x2_fill_HCPE.mat', 'T');
save('x3_fill_HCPE.mat', 'T');
save('y1_fill_HCPE.mat', 'data_algo1');
save('y2_fill_HCPE.mat', 'data_algo2');
save('y3_fill_HCPE.mat', 'data_algo3');
