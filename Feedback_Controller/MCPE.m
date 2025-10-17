% Clear Environment
clear; clc; close all;

% Set the random seed to ensure reproducibility
rng(0);

%% 1. Define the True State Space Model

A_true = [1.178,  0.001,  0.511, -0.403;
         -0.051, 0.661, -0.011,  0.061;
          0.076, 0.335,  0.560,  0.382;
          0.000, 0.335,  0.089,  0.849];

B_true = [0.004, -0.087;
         0.467,  0.001;
         0.213, -0.235;
         0.213, -0.016];

%% 2. Define Simulation Parameters

m = size(B_true, 2);    % Input dimension (m=2)
n = size(A_true, 1);    % State dimension (n=4)
N = 10;                 % Number of trajectories
L = 5;                  % Parameter L

% Assume alpha_i = 1
alpha = 0.1*ones(N,1);
alpha(1) = 1;

% Trajectory length range, determined by L and m
min_T = L;
max_T = (m + 1) * L - 1;  % T_p < (m+1)*L, i.e., max_T = 149

% Generate N trajectory lengths, ensuring L <= T1 <= T2 <= ... <= TN < (m+1)*L
T_i = sort(randi([min_T, max_T], N, 1));

rng(0);

%% 3. Calculate Δ_i

Delta = zeros(N, 1);  % Initialize Delta with N+1 elements
Delta(1) = L - 1; 
for i = 2:N
    Delta(i) = L - 1 - mod(sum(T_i(1:i-1) - L + 1), L);
    T_i_s(i) = sum(T_i(1:i-1) - L + 1);
end

%% 4. Initialize Cell Arrays to Store Trajectory Data

stateTrajectories = cell(N, 1);      % Store state for each trajectory, size is n x T_i
inputSequences = cell(N, 1);         % Store input for each trajectory, size is m x T_i
u_all = cell(N, 1);                  % Store control sequence for each trajectory
zBar = cell(N, 1);                   % Store \bar{z}_i for each trajectory

%% 5. Generate Multiple Trajectories

u_all = gen_MCPE_signal(m, T_i, L, N, 4);
u_all{1} = 10*u_all{1};
for i = 2:N
    u_all{i} = 10*u_all{i};
end

for i = 1:N
    T_current = T_i(i);
    % if (i == 1)
    %   u_all{i} = 0.01*u_all{i};
    % else
    %   u_all{i} = 0.01*u_all{i};  
    % end
    % Store control sequence into inputSequences
    inputSequences{i} = u_all{i};
 
    % Generate state trajectory
    x0_min = -1*ones(n, 1);
    x0_max = 1*ones(n, 1);
    x0 = x0_min + (x0_max - x0_min) .* rand(n, 1); % Randomly initialize initial state
    
    x = zeros(n, T_current);
    x(:, 1) = x0;
    
    % Generate process noise
    Q = 0.0 * eye(n); % Process noise covariance matrix
    w = mvnrnd(zeros(1, n), Q, T_current-1)'; % Each column is process noise for one time step
    
    % State update
    for k = 1:T_current-1
        x(:, k+1) = A_true * x(:, k) + B_true * u_all{i}(:, k) + w(:, k);
    end
    
    % Store state trajectory
    stateTrajectories{i} = x;
end

marn_x = [];
marn_x1 = [];
marn_u = [];
for trajIdx = 1:N
    T_current = size(stateTrajectories{trajIdx}, 2);
    if T_current < 2
        continue;  % Skip trajectories with insufficient length
    end
    
    % Extract state and input data
    x = stateTrajectories{trajIdx};
    u = inputSequences{trajIdx};
    
    % Build regression matrix [x1(k) x2(k) x3(k) x4(k) u1(k) u2(k)]
    x_cur = alpha(trajIdx) * [x(:, 1:T_current-1)];  % (T_current-1) x (n + m)
    u_cur = alpha(trajIdx) * u(:, 1:T_current-1);
    % Build target matrix [x1(k+1) x2(k+1) x3(k+1) x4(k+1)]
    x_next = alpha(trajIdx) * x(:, 2:T_current);  % (T_current-1) x n
    
    % Merge data
    marn_x = [marn_x x_cur];
    marn_x1 = [marn_x1 x_next];
    marn_u = [marn_u u_cur];
end

% Specific functions' meanings and usages will be explained later
%% Initialize LMI System
setlmis([]) 
%% Define Decision Variables
X = lmivar(2, [98, 4]); 
 
%% Add LMI Terms (only the upper right or lower left corner needs to be specified)
lmiterm([-1 1 1 X], marn_x, 1); 
lmiterm([-1 1 2 X], marn_x1, 1); 
lmiterm([-1 2 2 X], marn_x, 1);
lmisys = getlmis;
[aa, xopt] = gevp(lmisys, 1);
Q = dec2mat(lmisys, xopt, 1); % Obtain matrix X
K = marn_u * Q * inv(marn_x * Q);
abs(eig(A_true + B_true * K))
