function [ERROR, u_all] = function_CCPE(para)

%% ========== 0. Clear environment ==========

% Set random seed for reproducibility
rng(0);

process = para(1);
sv = para(2);

%% ========== 1. Define the true state-space model & basic parameters ==========
A_true = [1.178,  0.001,  0.511, -0.403;
         -0.051, 0.661, -0.011,  0.061;
          0.076, 0.335,  0.560,  0.382;
          0.000, 0.335,  0.089,  0.849];

B_true = [0.004, -0.087;
         0.467,  0.001;
         0.213, -0.235;
         0.213, -0.016];

n = size(A_true,1);     % State dimension
m = size(B_true,2);     % Input dimension
p = 10;                 % Number of controllers (trajectories)
L = 5;                  % Parameter L

% Assume alpha_i = 1
alpha = ones(p,1);

%% ========== 2. Define trajectory lengths T_i ==========
% Require mL <= T_i < (m+1)*L
min_T = m*L;          % Minimum trajectory length
max_T = (m+1)*L - 1;  % Maximum trajectory length

T_i = 25*ones(p,1); 
% T_i(i) represents the length of the i-th controller/trajectory

%% ========== 3. Initialize u_all{i} & state trajectory data structures ==========
u_all = cell(p,1);                % Control input sequences for each controller
stateTrajectories = cell(p,1);    % State trajectories for each controller
zBar = cell(p,1);                 % Store \bar{z}_i

%% ========== 4. Generate controller inputs (core part) ==========

u_all = gen_CCPE_signal(m, T_i, L, p, alpha, sv);

%% ========== 5. Generate state trajectories (multi-controller merge) ==========
% If generating a state sequence independently for each trajectory i, 
% do so for i = 1:p separately.
for i = 1:p
    T_current = T_i(i);
    x0 = -1 + 2*rand(n,1);   % Initial state
    x = zeros(n, T_current);
    x(:,1) = x0;

    % Process noise
    Q = process*eye(n);
    w = mvnrnd(zeros(n,1), Q, T_current-1)';  % (n x (T_current-1))

    % State update
    for k = 1:T_current-1
        % Use the control input of the i-th trajectory only:
        % x(k+1) = A*x(k) + B*u_all{i}(:,k) + w(k)
        x(:,k+1) = A_true*x(:,k) + B_true*u_all{i}(:,k) + w(:,k);
    end
    stateTrajectories{i} = x;
end

%% ========== 7. Ordinary Least Squares (OLS) Identification ==========

% Combine all trajectories to construct regression matrix & target matrix
Phi_combined = zeros(n + m, T_current - 1)';   % (total_samples x (n+m))
X_next_combined = zeros(n, T_current - 1)';     % (total_samples x n)

for i = 1:p
    T_current = T_i(i);
    x_i = stateTrajectories{i}; 
    u_i = u_all{i};

    if T_current < 2, continue; end
    % Construct the regression matrix: [x(k); u(k)]^T
    Phi_i = [x_i(:,1:T_current-1); u_i(:,1:T_current-1)]';  
    % => (T_current-1) x (n+m)
    
    % Target matrix: x(k+1)
    X_next_i = x_i(:,2:T_current)';  % (T_current-1) x n

    Phi_combined = Phi_combined + alpha(i) * Phi_i;
    X_next_combined = X_next_combined + alpha(i) * X_next_i;
end

% If insufficient data for regression, throw an error
if isempty(Phi_combined) || isempty(X_next_combined)
    error('Insufficient data for identification after merging.');
end

% Perform OLS
AB_est = Phi_combined \ X_next_combined;  % (n+m) x n

A_est = AB_est(1:n,:)';    % (n x n)
B_est = AB_est(n+1:end,:)';% (n x m)

disp('True A ='); disp(A_true);
disp('Estimated A ='); disp(A_est);
disp('True B ='); disp(B_true);
disp('Estimated B ='); disp(B_est);

%% ========== 8. Verification ==========

% Optionally, select the first trajectory for testing
selectedTraj = 1;
x_true = stateTrajectories{selectedTraj};
u_this = u_all{selectedTraj};
T_current = T_i(selectedTraj);

x_pred = zeros(n, T_current);
x_pred(:,1) = x_true(:,1);

for k = 1:T_current-1
    x_pred(:,k+1) = A_est * x_pred(:,k) + B_est * u_this(:,k);
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
