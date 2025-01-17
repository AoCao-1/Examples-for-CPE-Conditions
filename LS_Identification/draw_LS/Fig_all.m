% 清空环境
clear; clc; close all;
load('y1_fill_MCPE.mat')
load('y2_fill_MCPE.mat')
load('y3_fill_MCPE.mat')

load('x1_fill_MCPE.mat')
load('x2_fill_MCPE.mat')
load('x3_fill_MCPE.mat')


%===== 绘图部分 =====

subplot(2, 2, 1);
plot(T, data_algo1, '-o', ...   % 平均曲线
    'Color', 'r', ...
    'LineWidth', 2, ...
    'MarkerFaceColor','r');
hold on
plot(T, data_algo2, '-d', ...
    'Color', 'b', ...
    'LineWidth', 2, ...
    'MarkerFaceColor','b');
hold on
plot(T, data_algo3, '-s', ...
    'Color', 'g', ...
    'LineWidth', 2, ...
    'MarkerFaceColor','g');

set(gca, 'YScale', 'log');          % 将y轴设置为对数坐标



ylim([1e-3, 1e-1]);                 % 设置y轴范围
yticks([1e-3, 1e-2, 1e-1]);         % 选择想要显示的刻度
yticklabels({'10^{-3}','10^{-2}','10^{-1}'});  % 设置刻度标签
% 添加图例、标签、标题等
legend({'$\sigma_{w} =0.5$','$\sigma_{w} =0.05$','$\sigma_{w} =0.005$'}, 'Location', 'best','Interpreter', 'latex');
xlabel('Smallest singular value of $H^{mos}_{n+1}\left( \{ U_{i}\}^{N}_{i=1} \right) $','Interpreter', 'latex');
ylabel('$||{\hat G}-G||/||G||$', 'Interpreter', 'latex');




load('y1_fill_CCPE.mat')
load('y2_fill_CCPE.mat')
load('y3_fill_CCPE.mat')

load('x1_fill_CCPE.mat')
load('x2_fill_CCPE.mat')
load('x3_fill_CCPE.mat')

%===== 绘图部分 =====

subplot(2, 2, 2);
plot(T, data_algo1, '-o', ...   % 平均曲线
    'Color', 'r', ...
    'LineWidth', 2, ...
    'MarkerFaceColor','r');
hold on
plot(T, data_algo2, '-d', ...
    'Color', 'b', ...
    'LineWidth', 2, ...
    'MarkerFaceColor','b');
hold on
plot(T, data_algo3, '-s', ...
    'Color', 'g', ...
    'LineWidth', 2, ...
    'MarkerFaceColor','g');

set(gca, 'YScale', 'log');          % 将y轴设置为对数坐标

ylim([1e-3, 1e-1]);                 % 设置y轴范围
yticks([1e-3, 1e-2, 1e-1]);         % 选择想要显示的刻度
yticklabels({'10^{-3}','10^{-2}','10^{-1}'});  % 设置刻度标签
% 添加图例、标签、标题等
legend({'$\sigma_{w} =0.5$','$\sigma_{w} =0.05$','$\sigma_{w} =0.005$'}, 'Location', 'best','Interpreter', 'latex');
xlabel('Smallest singular value of $H^{cum}_{n+1}\left( \{ U_{i}\}^{N}_{i=1} \right) $','Interpreter', 'latex');
ylabel('$||{\hat G}-G||/||G||$', 'Interpreter', 'latex');






load('y1_fill_HCPE.mat')
load('y2_fill_HCPE.mat')
load('y3_fill_HCPE.mat')
load('x1_fill_HCPE.mat')
load('x2_fill_HCPE.mat')
load('x3_fill_HCPE.mat')
%===== 绘图部分 =====

subplot(2, 2, 3);
plot(T, data_algo1, '-o', ...   % 平均曲线
    'Color', 'r', ...
    'LineWidth', 2, ...
    'MarkerFaceColor','r');
hold on
plot(T, data_algo2, '-d', ...
    'Color', 'b', ...
    'LineWidth', 2, ...
    'MarkerFaceColor','b');
hold on
plot(T, data_algo3, '-s', ...
    'Color', 'g', ...
    'LineWidth', 2, ...
    'MarkerFaceColor','g');

set(gca, 'YScale', 'log');          % 将y轴设置为对数坐标

ylim([1e-3, 1e-1]);                 % 设置y轴范围
yticks([1e-3, 1e-2, 1e-1]);         % 选择想要显示的刻度
yticklabels({'10^{-3}','10^{-2}','10^{-1}'});  % 设置刻度标签
% 添加图例、标签、标题等
legend({'$\sigma_{w} =0.5$','$\sigma_{w} =0.05$','$\sigma_{w} =0.005$'}, 'Location', 'best','Interpreter', 'latex');
xlabel('Smallest singular value of $H^{hyb}_{n+1}\left( \{ U_{i}\}^{N}_{i=1} \right) $','Interpreter', 'latex');
ylabel('$||{\hat G}-G||/||G||$', 'Interpreter', 'latex');


load('y1_fill_MCPE_compare.mat')
load('y2_fill_MCPE_compare.mat')

load('x1_fill_MCPE_compare.mat')
load('x2_fill_MCPE_compare.mat')

%===== 绘图部分 =====

subplot(2, 2, 4);
plot(T, data_algo1, '-', ...   % 平均曲线
    'Color', 'r', ...
    'LineWidth', 2, ...
    'MarkerFaceColor','r');
hold on
plot(T, data_algo2, '-', ...
    'Color', 'b', ...
    'LineWidth', 2, ...
    'MarkerFaceColor','b');
hold on
plot(T(1:18:300), data_algo1(1:18:300), 'o', ...   % 平均曲线
    'Color', 'r', ...
    'LineWidth', 2, ...
    'MarkerFaceColor','r');

plot(T(1:18:300), data_algo2(1:18:300), 'o', ...   % 平均曲线
    'Color', 'b', ...
    'LineWidth', 2, ...
    'MarkerFaceColor','b');

% 添加图例、标签、标题等
set(gca, 'YScale', 'log');          % 将y轴设置为对数坐标

ylim([1e-3, 1e-1]);                 % 设置y轴范围
yticks([1e-3, 1e-2, 1e-1]);         % 选择想要显示的刻度
yticklabels({'10^{-3}','10^{-2}','10^{-1}'});  % 设置刻度标签

legend({'$\alpha_{1} =1$','$\alpha_{i} =1,i=\{2,...,10\}$'}, 'Location', 'best','Interpreter', 'latex');
xlabel('Changes in weighting factors','Interpreter', 'latex');
ylabel('$||{\hat G}-G||/||G||$', 'Interpreter', 'latex');

