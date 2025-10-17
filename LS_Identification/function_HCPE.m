function [ERROR, u_all] = function_HCPE(para)

%% ========== 0. Clear environment ==========
process = para(1);
sv = para(2);

rng(0);  % Set random seed

%% ========== 1. Define the true system model ==========
A_true = [1.178,  0.001,  0.511, -0.403;
         -0.051, 0.661, -0.011,  0.061;
          0.076, 0.335,  0.560,  0.382;
          0.000, 0.335,  0.089,  0.849];

B_true = [0.004, -0.087;
          0.467,  0.001;
          0.213, -0.235;
          0.213, -0.016];

n = size(A_true,1);  % State dimension
m = size(B_true,2);  % Input dimension

% ---- Set p, bar_p, L ----
p = 10;            % Total number of controllers
bar_p = 3;         % First bar_p controllers
L = 5;             % Length of the control signals
alpha_CCPE = ones(bar_p, 1);
% ---- Define alpha_i ----
alpha = ones(p, 1);

% ---- User-specified T_0 and T_i for each trajectory ----
% Assume T_0 is the same for the first bar_p controllers; 
% the rest have individual T_i
T_0 = 10;  % User-specified: if T_0 >= (m+1)*L => CCPE; else => rank=1 sum=0
% If T_0=25, (m+1)*L=3*5=15 => T_0>15 => CCPE
% If T_0=12, T_0 < 15 => rank=1+sum=0

% Lengths of the rest of the trajectories for controllers p - bar_p
T_rest = [7; 7; 8; 9; 10; 12; 12];  % size = (p - bar_p, 1) = 4 x 1
% Combine T_i
T_i = zeros(p, 1);
for i = 1:bar_p
    T_i(i) = T_0;
end
for i = bar_p + 1:p
    T_i(i) = T_rest(i - bar_p);
end

fprintf('Trajectory lengths T_i = '); disp(T_i(:)');

%% ========== 2. Initialize storage structures ==========
u_all = cell(p, 1);   % Control input
x_all = cell(p, 1);   % State trajectories

%% ========== 3. Generate input for the first bar_p controllers (Case1: CCPE / Case2: rank=1 & sum=0) ==========

if T_0 >= (m + 1) * L
    % ----------- Use CCPE method -----------
    fprintf('*** First %d controllers: Using CCPE method, because T_0=%d >= (m+1)*L=%d\n', ...
             bar_p, T_0, (m + 1) * L);

    u_CCPE = gen_CCPE_signal(m, T_0 * ones(bar_p, 1), L, bar_p, alpha_CCPE, sv);

    for i = 1:bar_p
        u_all{i} = u_CCPE{i};
    end
else
    % ----------- Use "rank=1 + sum=0" method -----------
    fprintf('*** First %d controllers: Using rank=1 & sum=0 method, because T_0=%d < (m+1)*L=%d\n', ...
             bar_p, T_0, (m + 1) * L);

    u_CCPE = gen_rank1_sum0_signal(m, T_0 * ones(bar_p, 1), L, bar_p, alpha_CCPE, sv);

    for i = 1:bar_p
        u_all{i} = u_CCPE{i};
    end
end

%% ========== 4. Merge signals from the first bar_p controllers as the initial signal ==========
% Here, we merge the input signals for the first bar_p controllers
% to create the initial signal: sum_{i=1}^{bar_p} alpha_i u_all{i} => init_signal
init_signal = zeros(m, T_0);
for k = 1:T_0
    sum_k = zeros(m, 1);
    for i = 1:bar_p
        sum_k = sum_k + alpha(i) * u_all{i}(:, k);
    end
    init_signal(:, k) = sum_k;
end

%% ========== 5. Use MCPE method for the remaining p - bar_p controllers ==========
fprintf('*** For the remaining %d controllers: Using MCPE method.\n', p - bar_p);

u_MCPE = gen_MCPE_signal_2(m, [T_0; T_i(bar_p + 1:p)], L, p - bar_p + 1, init_signal, sv);

for i = bar_p + 1:p
    u_all{i} = u_MCPE{i - bar_p + 1};
end

%% ========== 6. Generate state trajectories ==========
for i = 1:p
    T_cur = T_i(i);
    x_i = zeros(n, T_cur);
    x_i(:, 1) = randn(n, 1) * 0.5; % Initial state

    Q = process * eye(n);
    w = mvnrnd(zeros(n, 1), Q, T_cur - 1)';  % Process noise

    for k = 1:T_cur - 1
        x_i(:, k + 1) = A_true * x_i(:, k) + B_true * u_all{i}(:, k) + w(:, k);
    end
    x_all{i} = x_i;
end

%% ========== 7. Ordinary Least Squares (OLS) Identification (Merge all trajectories) ==========
Phi = [];  
X_next = [];
for i = 1:p
    T_cur = T_i(i);
    if T_cur < 2, continue; end
    X_i = x_all{i};
    U_i = u_all{i};

    Pi = [ X_i(:, 1:T_cur - 1); U_i(:, 1:T_cur - 1) ]';
    Xi_next = X_i(:, 2:T_cur)';

    Phi = [Phi; Pi];
    X_next = [X_next; Xi_next];
end

if isempty(Phi)
    error('Not enough data for identification!');
end

AB_est = Phi \ X_next;    % (n+m) x n
A_est = AB_est(1:n, :)';
B_est = AB_est(n+1:end, :)';

disp('True A ='); disp(A_true);
disp('Estimated A ='); disp(A_est);
disp('True B ='); disp(B_true);
disp('Estimated B ='); disp(B_est);

%% 10. Compute estimation error

% Compute error (MSE) for A and B matrices
ERROR_A = norm((A_true - A_est), 'fro');
ERROR_B = norm((B_true - B_est), 'fro');
ERROR = norm(([A_true B_true] - [A_est B_est]), "fro") / norm([A_true B_true], 'fro');
fprintf('MSE for A matrix: %.6f\n', ERROR_A);
fprintf('MSE for B matrix: %.6f\n', ERROR_B);
fprintf('MSE for [A B] matrix: %.6f\n', ERROR);

end
