function [ERROR,u_all] = function_MCPE_compare(para)

% Set random seed to ensure reproducibility
rng(0);
process = para(1);
sv = para(2);
al1 = para(3);
al2 = para(4);

%% 1. Define the true state-space model

A_true = [1.178,  0.001,  0.511, -0.403;
         -0.051, 0.661, -0.011,  0.061;
          0.076, 0.335,  0.560,  0.382;
          0.000, 0.335,  0.089,  0.849];

B_true = [0.004, -0.087;
         0.467,  0.001;
         0.213, -0.235;
         0.213, -0.016];

%% 2. Define simulation parameters

m = size(B_true, 2);    % Input dimension (m=2)
n = size(A_true, 1);    % State dimension (n=4)
N = 10;                 % Number of trajectories
L = 5;                  % Parameter L

% Assume alpha_i = 1
alpha = al1 * ones(N, 1);
alpha(1) = al2;

% Range of trajectory lengths, determined by L and m
min_T = L;
max_T = (m + 1) * L - 1; % T_p < (m+1)*L, i.e., max_T = 149

% Generate N trajectory lengths, ensuring L <= T1 <= T2 <= ... <= TN < (m+1)*L
T_i = sort(randi([min_T, max_T], N, 1));
rng(0);

%% 3. Compute Δ_i

Delta = zeros(N, 1); % Initialize Delta with N elements
Delta(1) = L - 1;
T_i_s(1) = 0;
for i = 2:N
    Delta(i) = L - 1 - mod(sum(T_i(1:i-1)-L+1), L);
    T_i_s(i) = sum(T_i(1:i-1) - L + 1);
end

%% 4. Initialize cell arrays to store trajectory data

stateTrajectories = cell(N, 1);      % Store state for each trajectory, size n x T_i
inputSequences = cell(N, 1);         % Store input for each trajectory, size m x T_i
u_all = cell(N, 1);                  % Store control sequence for each trajectory
zBar = cell(N, 1);                   % Store trajectory's \bar{z}_i

%% 5. Generate multiple trajectories

u_all = gen_MCPE_signal(m, T_i, L, N, sv);
u_all{1} = 10 * u_all{1};
for i = 2:N
    u_all{i} = 1 * u_all{i};
end
for i = 1:N
    T_current = T_i(i);
    % Store control sequence in inputSequences
    inputSequences{i} = u_all{i};
    
    % Generate state trajectory
    x0_min = -1 * ones(n, 1);
    x0_max = 1 * ones(n, 1);
    x0 = x0_min + (x0_max - x0_min) .* rand(n, 1); % Randomly initialize initial state
    
    x = zeros(n, T_current);
    x(:, 1) = x0;
    
    % Generate process noise
    Q = process * eye(n); % Process noise covariance matrix
    w = mvnrnd(zeros(1, n), Q, T_current-1)'; % Each column is a time step of process noise
    
    % State update
    for k = 1:T_current-1
        x(:, k+1) = A_true * x(:, k) + B_true * u_all{i}(:, k) + w(:, k);
    end
    
    % Store state trajectory
    stateTrajectories{i} = x;
    
    % Store \bar{z}_i (i.e., zBar{i})
    zBar_i = [];
    for p_step = 1:m
        idx = p_step * L - 1;
        if idx <= T_current - 1
            zBar_i = [zBar_i, u_all{i}(:, idx + 1)];
        end
    end
    zBar{i} = zBar_i;
end

%% 8. Ordinary Least Squares (OLS) Identification based on multiple trajectories

% Prepare identification data
% Initialize empty regression matrix and target matrix
Phi_combined = [];      % Regression matrix, size (total data points) x (n + m)
X_next_combined = [];   % Target matrix, size (total data points) x n

for trajIdx = 1:N
    T_current = size(stateTrajectories{trajIdx}, 2);
    if T_current < 2
        continue;  % Skip trajectories that are too short
    end
    
    % Extract state and input data
    x = stateTrajectories{trajIdx};
    u = inputSequences{trajIdx};
    
    % Construct regression matrix [x1(k) x2(k) x3(k) x4(k) u1(k) u2(k)]
    Phi_ext = alpha(trajIdx) * [x(:, 1:T_current-1)', u(:, 1:T_current-1)'];  % (T_current-1) x (n + m)
    
    % Construct target matrix [x1(k+1) x2(k+1) x3(k+1) x4(k+1)]
    X_next = alpha(trajIdx) * x(:, 2:T_current)';  % (T_current-1) x n
    
    % Combine data
    Phi_combined = [Phi_combined; Phi_ext];
    X_next_combined = [X_next_combined; X_next];
end

% Check if there is enough data for identification
if isempty(Phi_combined) || isempty(X_next_combined)
    error('Not enough data for identification. Please check the trajectory generation section.');
end

% Perform Ordinary Least Squares estimation
AB_est = Phi_combined \ X_next_combined;  % (n + m) x n

% Extract A_est and B_est
A_est = AB_est(1:n, :)';  % First n rows correspond to A's coefficients, transposed to (n x n)
B_est = AB_est(n+1:end, :)';  % Last m rows correspond to B's coefficients, transposed to (n x m)

% Display estimation results
disp('True A matrix:');
disp(A_true);
disp('Estimated A matrix:');
disp(A_est);

disp('True B matrix:');
disp(B_true);
disp('Estimated B matrix:');
disp(B_est);

%% 9. Validate the identification results

% Select one trajectory for comparison
selectedTrajIdx = 1;  % Choose the 1st trajectory

% Get the input signal for the selected trajectory
u_selected = inputSequences{selectedTrajIdx};  % m x T_selected

% Get trajectory length
T_selected = size(stateTrajectories{selectedTrajIdx}, 2);

% Regenerate predicted trajectory
x_pred = zeros(n, T_selected);
x_pred(:, 1) = stateTrajectories{selectedTrajIdx}(:, 1);  % Use the true initial state for the trajectory

for k = 1:T_selected-1
    x_pred(:, k+1) = A_est * x_pred(:, k) + B_est * u_selected(:, k);
end

%% 10. Compute estimation error

% Compute error (MSE) for A and B matrices
ERROR_A = norm((A_true - A_est), 'fro');
ERROR_B = norm((B_true - B_est), 'fro');
ERROR = norm(([A_true B_true] - [A_est B_est]), "fro") / norm([A_true B_true], 'fro');
fprintf('MSE for A matrix: %.6f\n', ERROR_A);
fprintf('MSE for B matrix: %.6f\n', ERROR_B);
fprintf('MSE for [A B] matrix: %.6f\n', ERROR);

end
