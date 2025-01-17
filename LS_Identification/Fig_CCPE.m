% Clear environment
clear; clc; close all;
ERROR_CCPE_all = [];
S_min = [];
N = 10;
T_i = 25*ones(N,1);  % Set time steps for each controller
m = 2;
L = 5;
alpha = ones(N,1);  % Weighting factors for each controller

% Loop for different values of sv
for sv = 0.1:0.5:5.6
    process = [0.5, sv];
    [ERROR_CCPE, u_CCPE_all] = function_CCPE(process);  % Call CCPE function
    ERROR_CCPE_all = [ERROR_CCPE_all; ERROR_CCPE];  % Collect error data

    T_max = max(T_i);  % Maximum length of time steps

    % Allocate and merge control input matrix, size m x T_max
    % The k-th column represents the merged input at time step k
    u_sum = zeros(m, T_max);

    % Loop through each time step
    for k = 1:T_max
        total_u_k = zeros(m, 1);  % Merged input at time step k
        for i = 1:N
            if k <= T_i(i)
                % If the i-th controller still has input at time k, accumulate
                total_u_k = total_u_k + alpha(i) * u_CCPE_all{i}(:, k);
            else
                % If k > T_i(i), the controller's input is considered 0
            end
        end
        % Store the result
        u_sum(:, k) = total_u_k;
    end

    % Reshape and create Hankel matrix
    U_all = reshape(u_sum, T_max*m, 1);
    H_CCPE = hankel_r(U_all, L, T_max-L+1, m);

    % Singular Value Decomposition (SVD)
    S = svd(H_CCPE);
    S_min = [S_min, min(S)];  % Store the minimum singular values
end

% Assume we have a horizontal axis T
T = S_min;  % Horizontal axis example

% Assuming data from three algorithms (e.g., multiple experiments for each T)
% Here, we'll simulate random data for demonstration purposes

data_algo1 = ERROR_CCPE_all;  % Data from the first algorithm

% Repeat the process for Algorithm 2 and Algorithm 3
ERROR_CCPE_all = [];
S_min = [];
for sv = 0.1:0.5:5.6
    process = [0.05, sv];
    [ERROR_CCPE, u_CCPE_all] = function_CCPE(process);
    ERROR_CCPE_all = [ERROR_CCPE_all; ERROR_CCPE];
    T_max = max(T_i);

    u_sum = zeros(m, T_max);
    for k = 1:T_max
        total_u_k = zeros(m, 1);
        for i = 1:N
            if k <= T_i(i)
                total_u_k = total_u_k + alpha(i) * u_CCPE_all{i}(:, k);
            end
        end
        u_sum(:, k) = total_u_k;
    end

    U_all = reshape(u_sum, T_max*m, 1);
    H_CCPE = hankel_r(U_all, L, T_max-L+1, m);

    S = svd(H_CCPE);
    S_min = [S_min, min(S)];
end

data_algo2 = ERROR_CCPE_all;  % Data from the second algorithm

% Repeat for the third algorithm
ERROR_CCPE_all = [];
S_min = [];
for sv = 0.1:0.5:5.6
    process = [0.005, sv];
    [ERROR_CCPE, u_CCPE_all] = function_CCPE(process);
    ERROR_CCPE_all = [ERROR_CCPE_all; ERROR_CCPE];
    T_max = max(T_i);

    u_sum = zeros(m, T_max);
    for k = 1:T_max
        total_u_k = zeros(m, 1);
        for i = 1:N
            if k <= T_i(i)
                total_u_k = total_u_k + alpha(i) * u_CCPE_all{i}(:, k);
            end
        end
        u_sum(:, k) = total_u_k;
    end

    U_all = reshape(u_sum, T_max*m, 1);
    H_CCPE = hankel_r(U_all, L, T_max-L+1, m);

    S = svd(H_CCPE);
    S_min = [S_min, min(S)];
end

data_algo3 = ERROR_CCPE_all;  % Data from the third algorithm

% Calculate mean and standard deviation for each algorithm at each T
mean_algo1 = mean(data_algo1, 2);
std_algo1 = std(data_algo1, 0, 2);

mean_algo2 = mean(data_algo2, 2);
std_algo2 = std(data_algo2, 0, 2);

mean_algo3 = mean(data_algo3, 2);
std_algo3 = std(data_algo3, 0, 2);

% Plotting section
figure; hold on; box on;

% (1) Shaded area for Algorithm 1
upper1 = mean_algo1 + std_algo1;
lower1 = mean_algo1 - std_algo1;
x1_fill_CCPE = [T, fliplr(T)];
y1_fill_CCPE = [upper1', fliplr(lower1')];

fill(x1_fill_CCPE, y1_fill_CCPE, 'r', ...
    'FaceAlpha', 0.2, ...
    'EdgeColor', 'none');
plot(T, mean_algo1, '-o', ...
    'Color', 'r', ...
    'LineWidth', 2, ...
    'MarkerFaceColor','r');

% (2) Shaded area for Algorithm 2
upper2 = mean_algo2 + std_algo2;
lower2 = mean_algo2 - std_algo2;
x2_fill_CCPE = [T, fliplr(T)];
y2_fill_CCPE = [upper2', fliplr(lower2')];

fill(x2_fill_CCPE, y2_fill_CCPE, 'b', ...
    'FaceAlpha', 0.2, ...
    'EdgeColor', 'none');
plot(T, mean_algo2, '-d', ...
    'Color', 'b', ...
    'LineWidth', 2, ...
    'MarkerFaceColor','b');

% (3) Shaded area for Algorithm 3
upper3 = mean_algo3 + std_algo3;
lower3 = mean_algo3 - std_algo3;
x3_fill_CCPE = [T, fliplr(T)];
y3_fill_CCPE = [upper3', fliplr(lower3')];

fill(x3_fill_CCPE, y3_fill_CCPE, 'g', ...
    'FaceAlpha', 0.2, ...
    'EdgeColor', 'none');
plot(T, mean_algo3, '-s', ...
    'Color', 'g', ...
    'LineWidth', 2, ...
    'MarkerFaceColor','g');

set(gca, 'YScale', 'log');  % Set y-axis to logarithmic scale
ylim([1e-3, 1e-1]);  % Set y-axis range
yticks([1e-3, 1e-2, 1e-1]);  % Choose desired tick marks
yticklabels({'10^{-3}', '10^{-2}', '10^{-1}'});  % Set tick labels

% Add legend, labels, and title
legend({'Algorithm1', 'Algorithm2', 'Algorithm3'}, 'Location', 'best');
xlabel('Length of each experiment T');
ylabel('Relative error');
title('Comparison of Different Methods');

% Save variables to .mat files
save('x1_fill_CCPE.mat', 'T');
save('x2_fill_CCPE.mat', 'T');
save('x3_fill_CCPE.mat', 'T');
save('y1_fill_CCPE.mat', 'data_algo1');
save('y2_fill_CCPE.mat', 'data_algo2');
save('y3_fill_CCPE.mat', 'data_algo3');
