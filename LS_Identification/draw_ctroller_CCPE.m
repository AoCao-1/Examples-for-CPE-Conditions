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
N = 10;                 % Number of trajectories
L = 5;                  % Parameter L

% Assume alpha_i = 1
alpha = ones(N,1);  % Weighting factors for each controller


rng(0);
T_i = 25*ones(N,1);  % Set time steps for each controller


%% 4. Initialize cell arrays to store trajectory data

stateTrajectories = cell(N, 1);      % Store state for each trajectory, size n x T_i
inputSequences = cell(N, 1);         % Store input for each trajectory, size m x T_i
u_CCPE_all = cell(N, 1);                  % Store control sequence for each trajectory
zBar = cell(N, 1);                   % Store trajectory's \bar{z}_i

%% 5. Generate multiple trajectories
u_CCPE_all = gen_CCPE_signal(m, T_i, L, N, alpha, 0.5);
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

save('u_CCPE_all');

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
    plot(tseg, u_CCPE_all{i}(1,:), '-', 'LineWidth', 1.6, 'Color', cmap(i,:));
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
    plot(tseg, u_CCPE_all{i}(2,:), '--', 'LineWidth', 1.6, 'Color', cmap(i,:));
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