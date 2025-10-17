%% ========== 0. Clean Environment ==========
clear; clc; close all;

% Set random seed for reproducibility
rng(0);

% Set options
tol_opt = 1e-4;

%% ========== 1. Define the True State-Space Model & Basic Parameters ==========
A_true = [1.178,  0.001,  0.511, -0.403;
         -0.051, 0.661, -0.011,  0.061;
          0.076, 0.335,  0.560,  0.382;
          0.000, 0.335,  0.089,  0.849];

B_true = [0.004, -0.087;
         0.467,  0.001;
         0.213, -0.235;
         0.213, -0.016];

C_true = eye(4);
D_true = zeros(4, 2);

n = size(A_true, 1);     % State dimension
m = size(B_true, 2);     % Input dimension
p = 10;                  % Number of controllers (trajectories)
L = 7;                   % CPE order L
n_p = 4;                 % Output dimension
T = 100;                 % Closed-loop horizon (simulation length)

% Assume alpha_i = 0.1
alpha = ones(p, 1);
M = 1;                   % Number of consecutive applications of optimal input (multi-step)
noise_max = 0.001;      % Measurement noise on y

%% ========== 2. Define Trajectory Lengths T_i ==========
min_T = L;
max_T = (m + 1) * L - 1;  % T_p < (m+1)*L, max_T = 149

% Generate N trajectory lengths, ensuring L <= T1 <= T2 <= ... <= TN < (m+1)*L
T_i = 30 * ones(p, 1); 
N = sum(T_i) - (p - 1) * (L - 1);
rng(0);

%% ========== 3. Compute Δ_i ==========
Delta = zeros(p, 1);  % Initialize Delta with N+1 elements
Delta(1) = L - 1; 
for i = 2:p
    Delta(i) = L - 1 - mod(sum(T_i(1:i-1) - L + 1), L);
    T_i_s(i) = sum(T_i(1:i-1) - L + 1);
end

%% ========== 4. Initialize Cell Arrays for Trajectory Data ==========
stateTrajectories = cell(p, 1);      % Store state trajectories, size n x T_i
inputSequences = cell(p, 1);         % Store input sequences, size m x T_i
u_all = cell(p, 1);                  % Store control sequences
zBar = cell(p, 1);                   % Store \bar{z}_i

%% ========== 5. Generate Multiple Trajectories ==========
u_all = gen_MCPE_signal(m, T_i, L, p, 0.2);

for i = 1:p
    T_current = T_i(i);
    inputSequences{i} = u_all{i};

    % Initialize initial state randomly
    x0_min = -ones(n, 1);
    x0_max = ones(n, 1);
    x0 = x0_min + (x0_max - x0_min) .* rand(n, 1);
    
    % Initialize state trajectory
    x = zeros(n, T_current);
    x(:, 1) = x0;
    
    % Generate process noise
    Q = 0.0 * eye(n);  % Process noise covariance matrix
    w = mvnrnd(zeros(1, n), Q, T_current - 1)';  % Process noise at each time step
    
    % State update loop
    for k = 1:T_current-1
        x(:, k+1) = A_true * x(:, k) + B_true * u_all{i}(:, k) + w(:, k);
    end
    
    % Store state trajectory
    stateTrajectories{i} = x;
end

% Construct Hankel matrices for each trajectory
for i = 1:length(u_all)
    U_all{i} = reshape(u_all{i}, T_i(i) * m, 1);
    X_all{i} = reshape(stateTrajectories{i}, T_i(i) * n, 1);
    H1{i} = hankel_r(U_all{i}, L, T_i(i) - L + 1, m);
    H2{i} = hankel_r(X_all{i}, L, T_i(i) - L + 1, n);
end

% Combine Hankel matrices
H_u = [];
H_y = [];
for i = 1:length(u_all)
    H_u = [H_u alpha(i) * H1{i}];
    H_y = [H_y alpha(i) * H2{i}];
end

% Terminal constraints
u_term = zeros(m, 1);
y_term = C_true * inv(eye(n) - A_true) * B_true * u_term;

%% ========== Set up MPC ==========
% Cost matrices
R = 1e-2 * eye(m);  % Input weighting
Q = 3 * eye(n_p);    % Output weighting
S = zeros(m, n_p);
Pi = [kron(eye(L), R) kron(eye(L), S);
      kron(eye(L), S') kron(eye(L), Q)];

% Cost for QP
cost_alpha = 0;
H = 2 * [cost_alpha * eye(N - L + 1), zeros(N - L + 1, (m + n_p) * L);
         zeros((m + n_p) * L, N - L + 1), Pi];
f = [zeros(N - L + 1, 1); -2 * kron(eye(L), R) * repmat(u_term, L, 1); 
     -2 * kron(eye(L), Q) * repmat(y_term, L, 1)];

% Inequality constraints
u_max = inf * ones(m, 1);
u_min = -inf * ones(m, 1);
y_max = inf * ones(n_p, 1);
y_min = -inf * ones(n_p, 1);
ub = [inf * ones(N - L + 1, 1); repmat(u_max, L, 1); repmat(y_max, L, 1)];
lb = [-inf * ones(N - L + 1, 1); repmat(u_min, L, 1); repmat(y_min, L, 1)];

% Equality constraints
Hu = H_u;
Hy = H_y;
B = [Hu, -eye(m * L), zeros(m * L, n_p * L);
     Hy, zeros(n_p * L, m * L), -eye(n_p * L);
     zeros(m * n, N - L + 1), [eye(m * n), zeros(m * n, m * (L - n))], zeros(m * n, n_p * L);
     zeros(n_p * n, N - L + 1), zeros(n_p * n, m * L), [eye(n_p * n), zeros(n_p * n, n_p * (L - n))]];

%% ========== Initial I/O Trajectories ==========
u_init = 0.8 * ones(m, n);  % Initial input
x0 = [0.4; 0.4; 0.5; 0.5]; % Initial state
x_init = zeros(n, n);
x_init(:, 1) = x0;

for i = 1:n-1
    x_init(:, i+1) = A_true * x_init(:, i) + B_true * u_init(:, i);
end

y_init = C_true * x_init;

% Closed-loop storage variables
u_cl = zeros(m, T);
u_cl(:, 1:n) = u_init;
y_cl = zeros(n_p, T);
y_cl(:, 1:n) = y_init;
y_cl_noise = y_cl;
x_cl = zeros(n, T);
x_cl(:, 1) = x0;

% Simulate first n steps
for j = 1:n
    x_cl(:, j+1) = A_true * x_cl(:, j) + B_true * u_cl(:, j);
end

% Open-loop storage variables
u_ol = zeros(m * L, T);
y_ol = zeros(n_p * L, T);
sigma_ol = zeros(n_p * L, T);
alpha_ol = zeros(N - L + 1, T);
u_init_store = zeros(m * n, T);
y_init_store = zeros(n_p * n, T);

% Candidate solution storage variables
u_cand = u_ol;
y_cand = y_ol;
alpha_cand = alpha_ol;
fval_cand = zeros(1, T);

sol_store = zeros((m + n_p) * L + N - L + 1, T);
sol_cand = sol_store;

tic;
%% ========== MPC Iterations ==========
for j = n+1:M:T
    % Update equality constraints
    c = [zeros((m + n_p) * L, 1); u_init(:); y_init(:)];
    
    %% Solve with quadprog
    [sol, fval(j), EXITFLAG] = quadprog(H, f, [], [], B, c, lb, ub, []);
    
    if EXITFLAG < 0
        error('Optimization problem not solved exactly...');
    end
    
    fval(j) = fval(j) + repmat(y_term, L, 1)' * kron(eye(L), Q) * repmat(y_term, L, 1) + ...
        repmat(u_term, L, 1)' * kron(eye(L), R) * repmat(u_term, L, 1);
    
    sol_store(:, j) = sol;
    alpha_ol(:, j) = sol(1:N - L + 1);
    u_ol(:, j) = sol(N - L + 2:N - L + 1 + m * L);
    y_ol(:, j) = sol(N - L + 2 + m * L:N - L + 1 + (m + n_p) * L);

    u_init_store(:, j - n) = u_init(:);
    y_init_store(:, j - n) = y_init(:);

    % Simulate closed loop
    for k = j:j + M - 1
        u_cl(:, k) = u_ol(n * m + 1 + (k - j) * m:n * m + m + (k - j) * m, j);
        x_cl(:, k + 1) = A_true * x_cl(:, k) + B_true * u_cl(:, k);
        y_cl(:, k) = C_true * x_cl(:, k) + D_true * u_cl(:, k);
        y_cl_noise(:, k) = y_cl(:, k) .* (ones(n_p, 1) + noise_max * (-1 + 2 * rand));
        
        % Set new initial conditions
        u_init = [u_init(:, 2:end), u_cl(:, k)];
        y_init = [y_init(:, 2:end), y_cl_noise(:, k)];
    end
end

elapsedTime = toc

%% Plot
figure
subplot(2, 2, 1)
plot(1:T, u_cl(1, 1:end))
hold on
plot(1:T, u_term(1) * ones(1, T))
legend('u_1', 'u_{1,eq}')

subplot(2, 2, 2)
plot(1:T, u_cl(2, 1:end))
hold on
plot(1:T, u_term(2) * ones(1, T))
legend('u_2', 'u_{2,eq}')

subplot(2, 2, 3)
plot(1:T, y_cl(1, 1:end))
hold on
plot(1:T, y_term(1) * ones(1, T))
legend('y_1', 'y_{1,eq}')

subplot(2, 2, 4)
plot(1:T, y_cl(2, 1:end))
hold on
plot(1:T, y_term(2) * ones(1, T))
legend('y_2', 'y_{2,eq}')
