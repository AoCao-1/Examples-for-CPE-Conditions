%% ========== 0. Clear Environment ==========
clear; clc; close all;

rng(0);  % Fix random seed

%% ========== 1. Define True System Model ==========
A_true = [1.178,  0.001,  0.511, -0.403;
             -0.051, 0.661, -0.011,  0.061;
              0.076, 0.335,  0.560,  0.382;
              0.000, 0.335,  0.089,  0.849];

B_true = [0.004, -0.087;
              0.467,  0.001;
              0.213, -0.235;
              0.213, -0.016];

n = size(A_true, 1);  % State dimension
m = size(B_true, 2);  % Input dimension

% ---- Set p, bar{p}, L ----
p = 10;            % Total number of controllers
bar_p = 3;         % First bar_p controllers
L = 5;
alpha_CCPE = ones(bar_p,1);  % CCPE alpha values
% ---- Define alpha_i ----
alpha = ones(p,1);

% ---- User specified T_0, and T_i for each trajectory ----
%   Assume T_0 is the trajectory length for the first bar_p controllers (same length); 
%   and the rest controllers have their own T_i
T_0 = 10;  % User specified: if T_0 >= (m+1)*L => CCPE; else => rank=1 sum=0
%   If T_0=25, (m+1)*L=3*5=15 => T_0>15 => CCPE
%   If T_0=12, it's <15 => rank=1+sum=0

% The remaining p - bar_p controllers have different trajectory lengths, here for example:
T_rest = [7; 7; 8; 9; 10; 12; 12];  % Size (p - bar_p, 1) = 4 x 1
% Combine T_i
% Signals for i=1..bar_p => T_0
%        for i=bar_p+1..p => T_rest(i - bar_p)
T_i = zeros(p,1);
for i=1:bar_p
    T_i(i) = T_0;
end
for i=bar_p+1:p
    T_i(i) = T_rest(i-bar_p);
end

fprintf('Trajectory lengths T_i = '); disp(T_i(:)');

%% ========== 2. Initialize Storage Structures ==========
u_all = cell(p,1);   % Control inputs
x_all = cell(p,1);   % State trajectories

%% ========== 3. Generate Inputs for First bar_p Controllers (Case1: CCPE / Case2: rank=1 & sum=0) ==========
if T_0 >= (m+1)*L
    % ----------- Use CCPE Method -----------
    fprintf('*** First %d controllers: Use CCPE method because T_0=%d >= (m+1)*L=%d\n', ...
             bar_p, T_0, (m+1)*L);

    u_CCPE = gen_CCPE_signal(m, T_0*ones(bar_p,1), L, bar_p, alpha_CCPE, 1);

    for i=1:bar_p
        % Call CCPE generation function (or embedded logic)
        % Here simplified: ensure non-zero at certain steps, others 0
        u_all{i} = u_CCPE{i};
    end

else
    % ----------- Use "rank=1 + sum=0" Design -----------
    fprintf('*** First %d controllers: Use rank=1 & sum=0 method because T_0=%d < (m+1)*L=%d\n', ...
             bar_p, T_0, (m+1)*L);

    u_CCPE = gen_rank1_sum0_signal(m, T_0*ones(bar_p,1), L, bar_p, alpha_CCPE, 1);

    for i=1:bar_p
        u_all{i} = u_CCPE{i};
    end
end

%% ========== 4. Use Combined Signals for First bar_p Controllers as Initial Signal ==========
% The sum_{i=1}^{bar_p} alpha_i z_i is used as the "initial signal"
% Here we combine the inputs for the first bar_p controllers
% User may want: sum_{i=1:bar_p} alpha_i u_all{i} => named init_signal
%   Size m x T_0 (assuming all bar_p have the same length)
init_signal = zeros(m, T_0);
for k = 1:T_0
    sum_k = zeros(m,1);
    for i=1:bar_p
        sum_k = sum_k + alpha(i)*u_all{i}(:,k);
    end
    init_signal(:,k) = sum_k;
end

%% ========== 5. Use MCPE Method for the Remaining p - bar_p Controllers ==========
fprintf('*** Remaining %d controllers: Use MCPE method.\n', p - bar_p);
    
u_MCPE = gen_MCPE_signal_2(m, [T_0; T_i(bar_p+1:p)], L, p - bar_p+1, init_signal, 1);

for i=bar_p+1 : p
    u_all{i} = u_MCPE{i-bar_p+1};
end

% Generate state trajectories
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
        % Use the i-th trajectory's corresponding control input => x(k+1)=A*x(k)+B*u_all{i}(:,k)+w(k)
        % If all p controllers are applied to the same trajectory, sum them
        x(:,k+1) = A_true*x(:,k) + B_true*u_all{i}(:,k) + w(:,k);
    end
    stateTrajectories{i} = x;
end

marn_x = 0;
marn_x1 = 0;
marn_u = 0;
N = p;

% Build regression matrices
for trajIdx = 1:N
    T_current = size(stateTrajectories{trajIdx}, 2);
    if T_current < 2
        continue;  % Skip trajectories that are too short
    end
    
    % Extract state and input data
    x = stateTrajectories{trajIdx};
    u = inputSequences{trajIdx};
    
    % Build regression matrix [x1(k) x2(k) x3(k) x4(k) u1(k) u2(k)]
    x_cur = alpha(trajIdx) * [x(:,1:T_current-1)];  % (T_current-1) x (n + m)
    u_cur = alpha(trajIdx) * u(:,1:T_current-1);
    % Build target matrix [x1(k+1) x2(k+1) x3(k+1) x4(k+1)]
    x_next = alpha(trajIdx) * x(:,2:T_current);  % (T_current-1) x n
    
    if trajIdx < bar_p+1
        % Combine data for the first bar_p controllers
        marn_x = marn_x + x_cur;
        marn_x1 = marn_x1 + x_next;
        marn_u = marn_u + u_cur;
    else
        marn_x = [marn_x x_cur];
        marn_x1 = [marn_x1 x_next];
        marn_u = [marn_u u_cur];    
    end
end

% Specific functions' meanings and usages will be explained later
%% Initialize LMI System
setlmis([]) 
%% Define Decision Variables
X = lmivar(2, [67, 4]); 

%% Add LMI Terms (only the upper right or lower left corner needs to be specified)
lmiterm([-1 1 1 X], marn_x, 1); 
lmiterm([-1 1 2 X], marn_x1, 1); 
lmiterm([-1 2 2 X], marn_x, 1);
lmisys = getlmis;
[aa, xopt] = gevp(lmisys, 1);
Q = dec2mat(lmisys, xopt, 1); % Obtain matrix X
K = marn_u * Q * inv(marn_x * Q);
abs(eig(A_true + B_true * K))
