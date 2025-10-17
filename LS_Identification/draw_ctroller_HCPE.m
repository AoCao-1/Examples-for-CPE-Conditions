clear; clc; close all;
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

    u_CCPE = gen_rank1_sum0_signal(m, T_0 * ones(bar_p, 1), L, bar_p, alpha_CCPE, 0.5);

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

u_MCPE = gen_MCPE_signal_2(m, [T_0; T_i(bar_p + 1:p)], L, p - bar_p + 1, init_signal, 0.5);

for i = bar_p + 1:p
    u_all{i} = u_MCPE{i - bar_p + 1};
end

u_HCPE_all = u_all;
N = 10;
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

save('u_HCPE_all');
  %%%%%%%%%%%%%%%%%%%%%%%%
%% ============ 美观且语义清晰的绘图（控制值 + 长度检查） ============
% 假定此时已有： m, L, N, T_i (Nx1), u_MCPE_all (cell array, 每项为 m x T_i(i) )
max_length = (m + 1) * L - 1;   % 长度阈值, 例如 14

% 时间轴总长度（把段串联在一条时间线上）
total_len = sum(T_i);

% 生成渐变颜色
cmap = parula(N);

% 创建 figure, 上下两行 + 一个小的长度条图放在下面（用 tiledlayout）
figure('Units','normalized','Position',[0.05 0.05 0.9 0.75]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
% 我们使用前两个大格作为两个输入的时间序列，第3格显示长度条形图
% ---- u1 ----
nexttile([1 1]); hold on;
start_idx = 0;
for i = 1:N
    tseg = start_idx + (1:T_i(i));        % 用 1-based 时间步更直观
    plot(tseg, u_HCPE_all{i}(1,:), '-', 'LineWidth', 1.6, 'Color', cmap(i,:));
    % 段中间标注长度
    text_mean = start_idx + round(T_i(i)/2);
    text(text_mean, max(ylim)-0.05*range(ylim), sprintf('T_%d=%d', i, T_i(i)),...
        'HorizontalAlignment','center','FontSize',9,'Color',cmap(i,:));
    % 竖线分隔（不在最后一段留最后界线）
    xline(start_idx + T_i(i), ':', 'Color',[0.7 0.7 0.7], 'LineWidth', 0.8);
    start_idx = start_idx + T_i(i);
end
xlabel('Time step (concatenated segments)');
ylabel('u_1');
title('Control input: dimension 1 (u_1)');
grid on; box on;
% 使图的 y 范围对称/留白更好看
yl = ylim; ylim([yl(1)-0.1*range(yl), yl(2)+0.12*range(yl)]);

% ---- u2 ----
nexttile([1 1]); hold on;
start_idx = 0;
for i = 1:N
    tseg = start_idx + (1:T_i(i));
    plot(tseg, u_HCPE_all{i}(2,:), '--', 'LineWidth', 1.6, 'Color', cmap(i,:));
    % 标注段序号或长度（可选）
    text_mean = start_idx + round(T_i(i)/2);
    text(text_mean, min(ylim)+0.05*range(ylim), sprintf('T_%d=%d', i, T_i(i)),...
        'HorizontalAlignment','center','FontSize',9,'Color',cmap(i,:));
    xline(start_idx + T_i(i), ':', 'Color',[0.7 0.7 0.7], 'LineWidth', 0.8);
    start_idx = start_idx + T_i(i);
end
xlabel('Time step (concatenated segments)');
ylabel('u_2');
title('Control input: dimension 2 (u_2)');
grid on; box on;
yl = ylim; ylim([yl(1)-0.1*range(yl), yl(2)+0.12*range(yl)]);

% ---- 底部：每段长度条形图 + 阈值线 ----
nexttile([1 1]); hold on;
barh_idx = 1:N;          % 段索引作横坐标
hbar = bar(barh_idx, T_i, 'FaceColor','flat');
for i = 1:N
    hbar.CData(i,:) = cmap(i,:);   % 用同一色图标记
end
% 画阈值线（水平虚线）
yline_val = max_length;
xlim([0.5 N+0.5]);
ylabel('Length (T_i)');
xlabel('Segment index');
title('Segment lengths vs threshold');
% 把阈值画成红色虚线并标注
plot([0.5 N+0.5], [yline_val yline_val], 'r--', 'LineWidth', 1.6);
text(N+0.6, yline_val, sprintf(' (m+1)L-1 = %d', yline_val), 'Color','r', 'FontWeight','bold');

% 在每个条上写长度数字
for i = 1:N
    text(i, T_i(i)+0.3, num2str(T_i(i)), 'HorizontalAlignment','center', 'FontSize',9);
end
grid on; box on;

% 检查是否有超出阈值的段（理论上不应有）
violations = find(T_i >= max_length);
if ~isempty(violations)
    % 高亮显示超出阈值的条目为红色，并给出提示
    for k = 1:length(violations)
        idx = violations(k);
        hbar.CData(idx,:) = [1 0.2 0.2]; % 红色
        text(idx, T_i(idx)+0.8, 'EXCEED', 'HorizontalAlignment','center', 'Color','r', 'FontWeight','bold');
    end
    warning('存在段长度 >= (m+1)L-1（请检查）: segments %s', mat2str(violations'));
else
    % 在图上写确认信息
    text(N+0.2, max(T_i)/2, 'All segments satisfy T_i < (m+1)L-1', 'Color',[0 0.5 0], 'FontWeight','bold');
end

% 美化整体标题与 colorbar (用作图例索引)
sgtitle('MCPE control sequences and segment lengths','FontSize',13,'FontWeight','bold');

% 可选：如果你想要整体时间轴坐标一致（方便对齐），可以把两个上方子图的 xlim 设为相同：
xlim_total = [1 total_len];
% 将两个子图的 xlim 统一（假如你愿意）
ax = findall(gcf,'Type','axes');
% ax 顺序可能不是我们期望的，按标题匹配设置 xlim（安全做法）
for k = 1:length(ax)
    if contains(get(get(ax(k),'Title'),'String'),'dimension 1') || contains(get(get(ax(k),'Title'),'String'),'dimension 2')
        xlim(ax(k), xlim_total);
    end
end

hold off;