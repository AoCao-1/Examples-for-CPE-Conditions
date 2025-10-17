%% ========== 0. Clear Environment ==========
clear; clc; close all;

% Set the random seed for reproducibility
rng(0);

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
D_true = zeros(4,2);

n = size(A_true,1);     % State dimension
m = size(B_true,2);     % Input dimension
p = 10;                  % Number of controllers (trajectories)
L = 7;                   % CPE order L
n_p = 4;                 % Output dimension
T = 100;                 % "Closed-loop horizon" (simulation length)

% Assume alpha_i = 0.1
alpha = ones(p,1);
M = 1;                   % Number of consecutive applications of optimal input (multi-step)
noise_max = 0.001;      % Measurement noise on y

%% ========== 2. Define Trajectory Length T_i ==========
% The requirement is mL <= T_i < (m+1)*L
% Trajectory length range, determined by L and m
min_T = L;
max_T = (m + 1) * L - 1; % T_p < (m+1)*L, i.e., max_T = 149
bar_p = 3;              % First bar_p controllers
alpha_CCPE = ones(bar_p,1);
T_0 = 30;               % User-defined: If T_0 >= (m+1)*L => CCPE; else => rank=1 sum=0
    % The following are the lengths for the p - bar_p remaining controllers, demonstrating:
    T_rest = 30*ones(p-bar_p);  % size = (p - bar_p, 1) = 4 x 1
    % Combine T_i
    % Signals i=1..bar_p => T_0
    %        i=bar_p+1..p => T_rest(i - bar_p)
    T_i = zeros(p,1);
    for i=1:bar_p
        T_i(i) = T_0;
    end
    for i=bar_p+1:p
        T_i(i) = T_rest(i-bar_p);
    end

N = sum(T_i(bar_p+1:p))-(p-bar_p)*(L-1)+T_0;

%% ========== 2. Initialize Storage Structures ==========
    u_all = cell(p,1);   % Control input
    x_all = cell(p,1);   % State trajectory

%% ========== 3. Generate Inputs for First bar_p Controllers (Case1: CCPE / Case2: rank=1 & sum=0) ==========
    if T_0 >= (m+1)*L
        % ----------- Use CCPE Method -----------
        fprintf('*** First %d controllers: Using CCPE method because T_0=%d >= (m+1)*L=%d\n', ...
                 bar_p, T_0, (m+1)*L);

        u_CCPE = gen_CCPE_signal(m, T_0*ones(bar_p,1), L, bar_p, alpha_CCPE,1);

        for i=1:bar_p
            % Directly call your CCPE generation function (or embedded logic)
            % This is simplified: only a few key steps ensure non-zero, others are zero
            u_all{i} = u_CCPE{i};
        end

    else
        % ----------- Use "rank=1 + sum=0" Design -----------
        fprintf('*** First %d controllers: Using rank=1 & sum=0 design because T_0=%d < (m+1)*L=%d\n', ...
                 bar_p, T_0, (m+1)*L);

        u_CCPE = gen_rank1_sum0_signal(m, T_0*ones(bar_p,1), L, bar_p, alpha_CCPE,1);

        for i=1:bar_p
            u_all{i} = u_CCPE{i};
        end
    end

%% ========== 4. Combine the First bar_p Controllers' Signals as Initial Signal ==========
% According to requirements, \sum_{i=1}^{\bar{p}} alpha_i z_i is used to create the "initial signal"
% Here, we can "superimpose" the inputs for p=bar_p
% The user might want: sum_{i=1:bar_p} alpha_i u_all{i} => called init_signal
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
fprintf('*** Last %d controllers: Using MCPE method.\n', p - bar_p);

u_MCPE = gen_MCPE_signal_2(m, [T_0;T_i(bar_p+1:p)], L, p - bar_p+1, init_signal,1);

for i=bar_p+1 : p
    u_all{i} = u_MCPE{i-bar_p+1};
end

for i=1:p
    T_cur = T_i(i);
    x_i = zeros(n, T_cur);
    x_i(:,1) = randn(n,1)*0.5; % Initialize state

    Q = 0.0*eye(n);
    w = mvnrnd(zeros(n,1), Q, T_cur-1)';  % Process noise

    for k=1:T_cur-1
        x_i(:, k+1) = A_true*x_i(:,k) + B_true*u_all{i}(:,k) + w(:,k);
    end
    x_all{i} = x_i;
end

T_max = T_0;

% Allocate matrix for combined control inputs, size m x T_max
% Later, the k-th column represents the combined input at time step k
u_sum = zeros(m, T_max);
x_sum = zeros(n, T_max);
% Loop through each time step
for k = 1:T_max
    % Sum weighted inputs for all controllers
    total_u_k = zeros(m,1);  % Combined value at time k
    total_x_k = zeros(n,1);  % Combined value at time k
    for i = 1:bar_p
        if k <= T_i(i)
            % If controller i has input at time k, accumulate it
            total_u_k = total_u_k + alpha(i) * u_all{i}(:, k);
            total_x_k = total_x_k + alpha(i) * x_all{i}(:, k);
        else
            % If k > T_i(i), that controller has insufficient length => treated as 0
        end
    end
    % Store the results
    u_sum(:, k) = total_u_k;
    x_sum(:, k) = total_x_k;
end

U_all{1} = reshape(u_sum,T_0*m,1);
X_all{1} = reshape(x_sum,T_0*n,1);
H1{1} = hankel_r(U_all{1},L,T_0-L+1,m);
H2{1} = hankel_r(X_all{1},L,T_0-L+1,n);

for i = bar_p+1:length(u_all)
    U_all{i-bar_p+1} = reshape(u_all{i},T_i(i)*m,1);
    X_all{i-bar_p+1} = reshape(x_all{i},T_i(i)*n,1);
    H1{i-bar_p+1} = hankel_r(U_all{i-bar_p+1},L,T_i(i)-L+1,m);
    H2{i-bar_p+1} = hankel_r(X_all{i-bar_p+1},L,T_i(i)-L+1,n);
end

H_u = [];
H_y = [];
for i = 1:(p-bar_p+1)
    H_u = [H_u H1{i}];
    H_y = [H_y H2{i}];
end

% Terminal constraints
u_term = 0*ones(m,1);
y_term = C_true*inv(eye(n)-A_true)*B_true*u_term;

%% Set up MPC
% Cost matrices
R = 1e-2*eye(m); % Input weighting
Q = 3*eye(n_p); % Output weighting
S = zeros(m,n_p);
Pi = [kron(eye(L),R) kron(eye(L),S);
      kron(eye(L),S') kron(eye(L),Q)];

% Cost for QP
cost_alpha=0;
H = 2*[cost_alpha*eye(N-L+1) zeros(N-L+1,(m+n_p)*L);
         zeros((m+n_p)*L,N-L+1) Pi]; % Cost
f = [zeros(N-L+1,1);-2*kron(eye(L),R)*repmat(u_term,L,1);-2*kron(eye(L),Q)*repmat(y_term,L,1)];

% Inequality constraints
u_max = inf*ones(m,1); u_min = -inf*ones(m,1);
y_max = inf*ones(n_p,1); y_min = -inf*ones(n_p,1);

ub = [inf(N-L+1,1);repmat(u_max,L,1);repmat(y_max,L,1)];
lb = [-inf(N-L+1,1);repmat(u_min,L,1);repmat(y_min,L,1)];

% Equality constraints
Hu = H_u;
Hy = H_y;

B = [Hu -eye(m*L) zeros(m*L,n_p*L); % Dynamics
     Hy zeros(n_p*L,m*L) -eye(n_p*L);
     zeros(m*n,N-L+1) [eye(m*n) zeros(m*n,m*(L-n))] zeros(m*n,n_p*L); % Initial input
     zeros(n_p*n,N-L+1) zeros(n_p*n,m*L) [eye(n_p*n) zeros(n_p*n,n_p*(L-n))]]; % Initial output

% Initial I/O trajectories
u_init = 0.8*ones(m,n); % Initial input
x0 = [0.4;0.4;0.5;0.5];
x_init = zeros(n,n); x_init(:,1) = x0;
for i = 1:n-1
   x_init(:,i+1) = A_true*x_init(:,i)+B_true*u_init(:,i);
end
y_init = C_true*x_init;

% Closed-loop storage variables
u_cl = zeros(m,T); u_cl(:,1:n) = u_init;
y_cl = zeros(n_p,T); y_cl(:,1:n) = y_init;
y_cl_noise = y_cl;
x_cl = zeros(n,T); x_cl(:,1) = x0;

% Simulate first nu steps
for j = 1:n
    x_cl(:,j+1) = A_true*x_cl(:,j) + B_true*u_cl(:,j);
end

% Open-loop storage variables
u_ol = zeros(m*L,T);
y_ol = zeros(n_p*L,T);
sigma_ol = zeros(n_p*L,T);
alpha_ol = zeros(N-L+1,T);
u_init_store = zeros(m*n,T);
y_init_store = zeros(n_p*n,T);
% Candidate solution storage variables
u_cand = u_ol;
y_cand = y_ol;
alpha_cand = alpha_ol;
fval_cand = zeros(1,T);

sol_store = zeros((m+n_p)*L+N-L+1,T);    
sol_cand = sol_store;
tic;
% MPC iterations
for j = n+1:M:T

    % Update equality constraints
    c = [zeros((m+n_p)*L,1);u_init(:);y_init(:)];

    %% Solve with quadprog
    [sol, fval(j), EXITFLAG] = quadprog(H,f,[],[],B,c,lb,ub,[]);
    if EXITFLAG<0%=1
        error('Optimization problem not solved exactly..') 
    end
    fval(j) = fval(j)+repmat(y_term,L,1)'*kron(eye(L),Q)*repmat(y_term,L,1)+repmat(u_term,L,1)'*kron(eye(L),R)*repmat(u_term,L,1);
    sol_store(:,j) = sol;
    alpha_ol(:,j) = sol(1:N-L+1);
    u_ol(:,j) = sol(N-L+2:N-L+1+m*L);
    y_ol(:,j) = sol(N-L+2+m*L:N-L+1+(m+n_p)*L);

    u_init_store(:,j-n) = u_init(:);
    y_init_store(:,j-n) = y_init(:);

    % Simulate closed-loop
    for k = j:j+M-1
        u_cl(:,k) = u_ol(n*m+1+(k-j)*m:n*m+m+(k-j)*m,j);
        x_cl(:,k+1) = A_true*x_cl(:,k) + B_true*u_cl(:,k);
        y_cl(:,k) = C_true*x_cl(:,k)+D_true*u_cl(:,k);
        y_cl_noise(:,k) = y_cl(:,k).*(ones(n_p,1)+noise_max*(-1+2*rand));
    % Set new initial conditions
        u_init = [u_init(:,2:end) u_cl(:,k)];
        y_init = [y_init(:,2:end) y_cl_noise(:,k)];
    end
end
elapsedTime = toc;

%% Plot
figure
subplot(2,2,1)
plot(1:T,u_cl(1,1:end))
hold on
plot(1:T,u_term(1)*ones(1,T))
legend('u_1','u_{1,eq}')

subplot(2,2,2)
plot(1:T,u_cl(2,1:end))
hold on
plot(1:T,u_term(2)*ones(1,T))
legend('u_2','u_{2,eq}')

subplot(2,2,3)
plot(1:T,y_cl(1,1:end))
hold on
plot(1:T,y_term(1)*ones(1,T))
legend('y_1','y_{1,eq}')

subplot(2,2,4)
plot(1:T,y_cl(2,1:end))
hold on
plot(1:T,y_term(2)*ones(1,T))
legend('y_2','y_{2,eq}')
