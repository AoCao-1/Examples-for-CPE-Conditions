%% ========== 0. Clear Environment ==========
clear; clc; close all;

% Set the random seed for reproducibility
rng(0);

%% ========== 1. Define True State Space Model & Basic Parameters ==========
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
alpha = 0.01*ones(p,1);

%% ========== 2. Define Trajectory Length T_i ==========
% Requirement: mL <= T_i < (m+1)*L
min_T = m*L;          % 20
max_T = (m+1)*L - 1;  % 29

T_i = 25*ones(p,1); 
% T_i(i) represents the length of the i-th controller/trajectory

%% ========== 3. Initialize u_all{i} & Trajectory Data Structures for Each Controller ==========
u_all = cell(p,1);                % Input sequences for each controller
stateTrajectories = cell(p,1);    % State trajectories for each controller
zBar = cell(p,1);                 % Store \bar{z}_i

%% ========== 4. Generate Controller Inputs (Core Part) ==========
u_all = gen_CCPE_signal(m, T_i, L, p, alpha, 1);

%% ========== 5. Generate State Trajectories (Merge Multiple Controllers) ==========
% If generating a state sequence for each "trajectory i" independently,
% loop over i = 1:p
for i = 1:p
    T_current = T_i(i);
    inputSequences{i} = u_all{i};
    x0 = -1 + 2*rand(n,1);   % Initial state
    x = zeros(n, T_current);
    x(:,1) = x0;

    % Process noise
    Q = 0*eye(n);
    w = mvnrnd(zeros(n,1), Q, T_current-1)';  % (n x (T_current-1))

    % State update
    for k = 1:T_current-1
        % Only use the control input corresponding to trajectory i => x(k+1) = A*x(k) + B*u_all{i}(:,k) + w(k)
        % To combine p controllers into a single trajectory, sum them
        x(:,k+1) = A_true*x(:,k) + B_true*u_all{i}(:,k) + w(:,k);
    end
    stateTrajectories{i} = x;
end

marn_x = 0;
marn_x1 = 0;
marn_u = 0;
N = p;
for trajIdx = 1:N
    T_current = size(stateTrajectories{trajIdx},2);
    if T_current < 2
        continue;  % Skip trajectories with insufficient length
    end
    
    % Extract state and input data
    x = stateTrajectories{trajIdx};
    u = inputSequences{trajIdx};
    
    % Build the regression matrix [x1(k) x2(k) x3(k) x4(k) u1(k) u2(k)]
    x_cur = alpha(trajIdx)*[x(:,1:T_current-1)];  % (T_current-1) x (n + m)
    u_cur = alpha(trajIdx)*u(:,1:T_current-1);
    % Build the target matrix [x1(k+1) x2(k+1) x3(k+1) x4(k+1)]
    x_next = alpha(trajIdx)*x(:,2:T_current);  % (T_current-1) x n
    
    % Merge the data
    marn_x = marn_x + x_cur;
    marn_x1 = marn_x1 + x_next;
    marn_u = marn_u + u_cur;
end

% The purpose and usage of specific functions will be explained later
%% Initialize LMI System
setlmis([]) 
%% Define Decision Variables
X = lmivar(2, [24, 4]); 
 
%% Add LMI Terms (only the upper right or lower left corner needs to be specified)
lmiterm([-1 1 1 X], marn_x, 1); 
lmiterm([-1 1 2 X], marn_x1, 1); 
lmiterm([-1 2 2 X], marn_x, 1);
lmisys = getlmis;
[aa, xopt] = gevp(lmisys, 1);
Q = dec2mat(lmisys, xopt, 1); % Obtain matrix X
K = marn_u * Q * inv(marn_x * Q);
abs(eig(A_true + B_true * K))
