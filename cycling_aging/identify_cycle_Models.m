%% 自动识别锂离子电池循环老化程序
% 2024/08/05 liuyang
% Ref from Paul.Gasper@nrel.gov
% Article: Paul Gasper et al 2022 J. Electrochem. Soc. 169 080518

clear; clc; close all;
addpath('Functions')
addpath('slanCM') % 绘图工具https://github.com/slandarer/slanColor

% 结果对比
% load fitted_models

% 在引导重采样过程中,有时会出现病态的雅可比矩阵
% 通过查看90%的置信区间,这些结果被忽略了
warning('Off')

%% Loading and formatting data
% 原始循环老化数据:
load Naumann data

%% 在data表格加入U_a:石墨平衡电位
% 石墨平衡电位来自Safari and Delacourt,  2011, JES, https://doi.org/10.1149/1.3567007. 
% 石墨锂离子化学计量来自Schimpe et. al. using graphite/lithium coin cell meausrements.
% x_a: stoichiometric fraction of lithium intercalated into graphite as a
% function of state-of-charge (SOC).

x_a = @(soc) 8.5e-3 + soc.*(7.8e-1 - 8.5e-3);
% U_a: graphite to lithium open cirucit potential as a function of x_a
U_a = @(x_a) 0.6379 + 0.5416.*exp(-305.5309.*x_a) + 0.044.*tanh(-1.*(x_a-0.1958)./0.1088) - 0.1978.*tanh((x_a-1.0571)./0.0854) - 0.6875.*tanh((x_a+0.0117)./0.0529) - 0.0175.*tanh((x_a-0.5692)./0.0875);
% Add a column to the table:
data.U_a = U_a(x_a(data.soc)); % 原始数据table最后一列加入负极平衡电势

%%
% 创建训练集和验证集:
cellNums = unique(data.cellNum); % 根据电池数量生成cellNum×1的数组

% Consistent rng seed for reproducible test/train sets:
rng(6686)
data_partition = cvpartition(length(cellNums),'HoldOut',0.3); % 交叉验证数据集划分
trainset = cellNums(training(data_partition));
data_train = data(any(data.cellNum==trainset',2),:);
validationset = cellNums(test(data_partition));
data_validation = data(any(data.cellNum==validationset',2),:);

% 创建外推数据,EFCs为12000圈:
sim_cellNums = [21; 22; 23]; % 新创建的电池编号
temps = 40; % 新创建的工况温度
soc = 0.5; % 新工况SOC
dod = 0.5;
Cavg = [0.5, 1, 1.5];
Ua = U_a(x_a(soc));

t = [0:100:12000]'; % 12000圈
data_sim = table(); % 创建空的做外推的table
for i = 1:length(sim_cellNums)
    temp = table(); % 暂存table
    temp.cellNum = sim_cellNums(i).*ones(length(t),1);
    temp.t = t; % 时间单位为EFCs
    temp.t_days = t./(12/Cavg(i)); % 单位为day
    temp.TdegC = temps.*ones(length(t),1); % °C温度
    temp.TdegK = (temps+273.15).*ones(length(t),1);
    temp.soc = soc.*ones(length(t),1);
    temp.dod = dod.*ones(length(t),1);
    temp.Cavg = Cavg(i).*ones(length(t),1);
    temp.U_a = Ua.*ones(length(t),1);
    data_sim = [data_sim; temp];
end

% 删除除data data_train data_validation data_sim外的其它变量
clearvars -except data data_train data_validation data_sim 

%% 绘制原始数据容量衰减图
cellNums = data.cellNum;
x = data.t;
x_days = data.tdays;
y = data.qdis;
s = 20;

cycle_aging_raw_plot(x, x_days, y, data, cellNums, s);
   
clearvars -except data data_train data_validation data_sim

%% Naumann et. al., 2020 J. Power Sources 451 227666模型
disp("开始复现文献模型 ref. Naumann et. al., 2020 JPS.") % 复现base模型, q=1-(a*Crate+b)*(c*(DOD-0.6)^3+d)*EFC^0.5
% Some constants:
TdegK_ref = 298.15; % 25 deg C
Ua_ref = 0.123; % V, Ua @ 50% SOC
Rug = 8.314; % Universal gas constant
F = 96485; % Faraday's constant

% base模型 q=1-(a*Crate+b)*(c*(DOD-0.6)^3+d)*EFC^0.5
% Input vars: EFC, Crate, DOD 输入变量
% Input params: a,b,c,d 输入参数
k_cyc = @(p,x) (p(1).*x(:,1) + p(2)).*(p(3).*(x(:,2) - 0.6).^3 + p(4));
p_k_cyc = [0.063,0.0971,4.0253,1.0923]; % 子模型系数a,b,c,d(β1)

% Capacity loss model:
Naumann_model = @(p,x) 1 - k_cyc(p,x(:,2:3)).*x(:,1).^(0.5)/100; % x1=>EFCs, x2=>Crate, x3=>DOD
p_Naumann_model = p_k_cyc;

% Evaluate predictions for the entire data set or only DOD cases
% 所有数据拟合结果
x = [data.t, data.Cavg, data.dod]; % 影响因素
qdis_pred = Naumann_model(p_Naumann_model,x);
% Residuals
R = data.qdis - qdis_pred; % 残差,真实值-预测值
% Mean absolute error:
MAE = sum(abs(R))/length(R);
fprintf("全部数据拟合MAE: %0.3g%%\n", MAE*100);

%%
% 绘制Naumann_model拟合结果:
cellNums = unique(data.cellNum)';
% colors = lines(length(cellNums));
colors = colormap(slanCM('Spectral')); %'RdBu'
colors = colors(1:12:end,:);

ax1 = subplot(3,2,1:4); hold on; box on; grid on; % capacity vs time
ax2 = subplot(3,2,5); hold on; box on; grid on; % residuals vs time
ax3 = subplot(3,2,6); box on; grid on; % residuals histogram
for cellNum = cellNums
    mask = data.cellNum == cellNum;
    p1 = plot(ax1, data.t(mask), data.qdis(mask), 'ok', 'MarkerFaceColor', colors(cellNum,:), 'MarkerSize', 6);
    p2 = plot(ax1, data.t(mask), qdis_pred(mask), '-', 'Color', colors(cellNum,:), 'LineWidth', 2);
    plot(ax2, data.t(mask), R(mask), '-', 'Color', colors(cellNum,:), 'LineWidth', 2);
end
histogram(ax3, R, 'BinWidth', 0.02, 'Orientation', 'horizontal', 'FaceColor', 'k')
% Format residuals plots:
yline(ax2, 0, '--k', 'LineWidth', 2);
RLim = max(abs(ax2.YLim)); ax2.YLim = [-RLim, RLim];
yline(ax3, 0, '--k', 'LineWidth', 2); ax3.YLim = [-RLim, RLim];
% Plot decorations:
xlabel(ax1, 'EFCs'); ylabel(ax1, 'Relative capacity');
legend(ax1, [p1,p2], {'Data','EFCs^{0.5}'}, 'Location', 'southwest')
title(ax1,sprintf("Baseline model, q_{dis} = 1 - k(C_{rate},DOD)*EFCs^{0.5}, MAE=%0.3g%%", MAE*100))
xlabel(ax2, 'EFCs'); ylabel(ax2, 'Residual error');
xlabel(ax3, 'Counts'); ylabel(ax3, 'Residual error'); grid on;

set(gcf, 'units','inches','PaperPosition', [0 0 12 10]);
print(gcf, 'Naumann_model','-r600','-dpng') % 注意更改保存位置

%%
% 仅评估DOD组,样本7-13
data_DOD = table();
cellNums = data.cellNum;
for cellNum = unique(cellNums,'stable')'
    if (cellNum >= 7) && (cellNum <= 13)
        mask = cellNums == cellNum;
        data_temp = table();
        data_temp.cellNum = data.cellNum(mask);
        data_temp.t = data.t(mask);
        data_temp.Cavg = data.Cavg(mask);
        data_temp.dod = data.dod(mask);
        data_temp.qdis = data.qdis(mask);
        data_DOD = [data_DOD; data_temp];
    end
end
x = [data_DOD.t, data_DOD.Cavg, data_DOD.dod]; % 影响因素  
qdis_pred = Naumann_model(p_Naumann_model,x);
% Residuals
R = data_DOD.qdis - qdis_pred; % 残差,真实值-预测值
% Mean absolute error:
MAE = sum(abs(R))/length(R);
fprintf("全部数据拟合MAE: %0.3g%%\n", MAE*100);

%%
cellNums = unique(data_DOD.cellNum)';
colors = colormap(slanCM('RdBu')); %'RdBu'
colors = colors(1:20:end,:);

ax1 = subplot(3,2,1:4); hold on; box on; grid on; % capacity vs time
ax2 = subplot(3,2,5); hold on; box on; grid on; % residuals vs time
ax3 = subplot(3,2,6); box on; grid on; % residuals histogram
for cellNum = cellNums
    mask = data_DOD.cellNum == cellNum;
    p1 = plot(ax1, data_DOD.t(mask), data_DOD.qdis(mask), 'ok', 'MarkerFaceColor', colors(cellNum,:), 'MarkerSize', 6);
    p2 = plot(ax1, data_DOD.t(mask), qdis_pred(mask), '-', 'Color', colors(cellNum,:), 'LineWidth', 2);
    plot(ax2, data_DOD.t(mask), R(mask), '-', 'Color', colors(cellNum,:), 'LineWidth', 2);
end
histogram(ax3, R, 'BinWidth', 0.02, 'Orientation', 'horizontal', 'FaceColor', 'k')
% Format residuals plots:
yline(ax2, 0, '--k', 'LineWidth', 2);
RLim = max(abs(ax2.YLim)); ax2.YLim = [-RLim, RLim];
yline(ax3, 0, '--k', 'LineWidth', 2); ax3.YLim = [-RLim, RLim];
% Plot decorations:
xlabel(ax1, 'EFCs'); ylabel(ax1, 'Relative capacity');
legend(ax1, [p1,p2], {'Data','EFCs^{0.5}'}, 'Location', 'southwest')
title(ax1,sprintf("Baseline model, q_{dis} = 1 - k(C_{rate},DOD)*EFCs^{0.5}, MAE=%0.3g%%", MAE*100))
xlabel(ax2, 'EFCs'); ylabel(ax2, 'Residual error');
xlabel(ax3, 'Counts'); ylabel(ax3, 'Residual error'); grid on;

set(gcf, 'units','inches','PaperPosition', [0 0 12 10]);
print(gcf, 'Naumann_model_DOD','-r600','-dpng') % 注意更改保存位置


%% 利用交叉验证和自举检验重采样优化模型
% Optimize with cross-validation and bootstrap resampling:
% 训练集数据
x = [data_train.t, data_train.Cavg, data_train.dod];
y = data_train.qdis;
cellNums = data_train.cellNum; % 进行训练的样本电池编号

fitOpt.CV = 'LeaveOut'; fitOpt.bootstrap = 'On'; fitOpt.bootstrapIterations = 1000; % 留一法
Naumann_train = optimize(x, y, cellNums, Naumann_model, p_Naumann_model, fitOpt); % 优化算法,输出的结果
disp("通过交叉验证和自举重采样重新优化模型...")
fprintf("优化后的模型在训练集MAE: %0.3g%%\n", Naumann_train.MAE*100); % 在训练集上为2.99%

clearvars -except data data_train data_validation data_sim...
    TdegK_ref Ua_ref Rug F Naumann_model Naumann_train...
    data_DOD

%% Reoptimization of baseline model
% 在训练数据集上对基准模型q=1-β1*t^β2进行重新优化,以便与自动识别的模型进行比较
disp(" ")
disp("优化后的降解EFCs^0.5模型")

% Step.1 简化的q=1-β1*t^0.5拟合全局参数β1
% Define a local equation for determining the time-invariant capacity fade rates, beta_1.
% Input vars: t
% Input params: beta_1
sqrt_model = @(p,x) 1 - p(1).*x(:,1).^(0.5); % 简化的t^0.5模型
p0_sqrt_model = 0.001; % 初始值,针对每个case优化得到一个β1

% Assemble data:
x = data_train.t;
y = data_train.qdis;
cellNums = data_train.cellNum;
% Optimize:
sqrt_local_fit = optimize_local(x, y, cellNums, sqrt_model, p0_sqrt_model, []);
fprintf("t^0.5局部模型对训练集MAE: %0.3g%%\n", sqrt_local_fit.MAE*100); % 2.17%

% t^0.5模型拟合效果
% Plot the local fit result:
plotOpt.layout = 'single axis';
plotOpt.labels = {'EFCs^{0.5} fit'};
plotOpt.colors = {'r'};
plotOpt.xlabel = 'EFCs';
plotOpt.ylabel = 'Relative capacity';
plotOpt.title = sprintf("Local fits of q_{dis} = 1 - \\beta_1*EFCs^{0.5}, MAE=%0.3g%%", sqrt_local_fit.MAE*100);
plot_capacity_fits(x(:,1), sqrt_local_fit, data_train, plotOpt)
set(gcf, 'units','inches','PaperPosition', [0 0 12 8]);
print(gcf, 'sqrt_local_train','-r600','-dpng') % 注意更改保存位置

% Step.2 根据Step.1的β1再训练子模型β1向量
% 带有ArrTfl_mod子模型的beta_1模型:
% 可以通过删除常数来简化Naumann的模型
% 这使得参数单位的可解释性稍差,但使模型的结构易于与自动识别的模型进行比较
% Input vars: TdegK, U_a, Crate, DOD
% Input params: γ1,γ2,γ3,γ4,γ5,γ6,γ7,γ8
ArrTfl_mod_model = @(p,x) p(1).*exp(p(2).*(1./x(:,1))).*(exp(p(3).*x(:,2))+p(4))...
    .*(p(5).*x(:,3) + p(6)).*(p(7).*(x(:,4) - 0.6).^3 + p(8));

p0_ArrTfl_mod = [0.001, 0, 0, 0, 0.1382, 0.0599, 0.4968, 1.2796];

% 组装数据:
% 每个数据系列只有一个beta_1值,获取每个数据系列的不变数据变量的值进行训练:
submodel_data = assemble_invariant_data(data_train);
x = [submodel_data.TdegK, submodel_data.U_a, submodel_data.Cavg, submodel_data.dod];
y = sqrt_local_fit.p(:,1); % 来自t^0.5拟合的β1值
cellNums = submodel_data.cellNum;

% Optimize with cross-validation:
fitOpt.CV = 'LeaveOut';
ArrTfl_mod_train = optimize(x, y, cellNums, ArrTfl_mod_model, p0_ArrTfl_mod, fitOpt);
% Plot sub-model fit result:
figure; hold on; box on; grid on;
plot(submodel_data.TdegC, ArrTfl_mod_train.y, 'ok', 'LineWidth', 1.5)
plot(submodel_data.TdegC, ArrTfl_mod_train.y_fit, 'xr', 'LineWidth', 1.5)
xlabel('Temperature (\circC)'); ylabel('\beta_1 (EFCs^{-0.5})');
legend('Locally fit values', 'Sub-model prediction', 'Location', 'northwest')
title(sprintf("Modified Arrhenius Tafel sub-model, MAPE=%0.3g%%", ArrTfl_mod_train.MAPE*100))


% Step.3 用训练好的子模型β1去拟合子模型的局部参数γi
% Use the beta_1 sub-model to construct a global model
% Input vars: t, TdegK, U_a
% Input params: beta_1(γi)
sqrt_ArrTflmod_model = @(p,x) 1 - ArrTfl_mod_model(p, x(:,2:5)).*x(:,1).^(0.5);
p0_sqrt_ArrTflmod_model = ArrTfl_mod_train.p;

% Assemble data:
x = [data_train.t, data_train.TdegK, data_train.U_a, data_train.Cavg, data_train.dod];
y = data_train.qdis;
cellNums = data_train.cellNum;
% Optimize with cross-validation and bootstrap resampling:
fitOpt.CV = 'LeaveOut'; fitOpt.bootstrap = 'On'; fitOpt.bootstrapIterations = 1000;
disp("通过交叉验证和自举重采样重新优化t^0.5模型...")
sqrt_ArrTflmod_train = optimize(x, y, cellNums, sqrt_ArrTflmod_model, p0_sqrt_ArrTflmod_model, fitOpt);
fprintf("t^0.5全局优化模型在训练集MAE: %0.3g%%\n", sqrt_ArrTflmod_train.MAE*100);

% Plot global fit results:绘制全局优化结果
plotOpt.layout = 'individual axes';
plotOpt.labels = {'Fit'};
plotOpt.colors = {'r'};
plotOpt.xlabel = 'EFCs';
plotOpt.ylabel = 'Relative capacity';
plotOpt.title = sprintf("q_{dis} = 1 - ArrTfl_{mod}*EFCs^{0.5}, MAE=%0.3g%%", sqrt_ArrTflmod_train.MAE*100);
plotOpt.confidenceInterval = [5 95];
plot_capacity_fits(x(:,1), sqrt_ArrTflmod_train, data_train, plotOpt)
set(gcf, 'units','inches','PaperPosition', [0 0 25 30]);
print(gcf, 'sqrt_ArrTflmod_train','-r600','-dpng');

% 局部模型的系统误差会传播到全局模型
% 由于beta_1子模型的不准确,全局模型的总体MAE为2.6%,而局部模型的MAE为2.17%
% Bootstrapping gives useful information such as the distributions of parameter values.
figure; plotmatrix(sqrt_ArrTflmod_train.p_boot);

% 由于某些引导迭代上的条件不好的雅可比矩阵
% 一些参数具有极值,90%置信区间应该看起来更正常
p_boot_90CI = [];
for i = 1:size(sqrt_ArrTflmod_train.p_boot,2)
    p_boot = sqrt_ArrTflmod_train.p_boot(:,i);
    CI = prctile(p_boot, [5 95]);
    mask = p_boot >= CI(1) & p_boot <= CI(2);
    p_boot_90CI = [p_boot_90CI, p_boot(mask)];
end
figure; plotmatrix(p_boot_90CI);
% 参数值看起来不相关
% 这意味着每个参数都对数据的不同特征进行建模

% Step.4 验证集评估
% Evaluate the data on the validation set:
x = [data_validation.t, data_validation.TdegK, data_validation.U_a, data_validation.Cavg, data_validation.dod];
y = data_validation.qdis;
sqrt_ArrTflmod_validation = evaluate(x, y, sqrt_ArrTflmod_model, sqrt_ArrTflmod_train);
fprintf("t^0.5全局优化模型在验证集: %0.3g%%\n", sqrt_ArrTflmod_validation.MAE*100); % 2.21%
plotOpt.title = sprintf("q_{dis} = 1 - ArrTfl_{mod}*EFCs^{0.5}, Validation MAE=%0.3g%%", sqrt_ArrTflmod_validation.MAE*100);
plot_capacity_fits(x(:,1), sqrt_ArrTflmod_validation, data_validation, plotOpt)
set(gcf, 'units','inches','PaperPosition', [0 0 25 18]);
print(gcf, 'sqrt_ArrTflmod_validation','-r600','-dpng');


% Step.5 外推数据,新工况预测12000圈
% Simulate 12000 cycs aging:
x = [data_sim.t, data_sim.TdegK, data_sim.U_a, data_sim.Cavg, data_sim.dod];
sqrt_ArrTflmod_sim = simulate(x, sqrt_ArrTflmod_model, sqrt_ArrTflmod_train);
% Plot the simulation result:
x = data_sim.t;
plotOpt.xlabel = 'EFCs';
plotOpt.title = '1,2000 cycs simulation, q_{dis} = 1 - ArrTfl_{mod}*EFCs^{0.5}';
plot_capacity_sim(x, sqrt_ArrTflmod_sim, data_sim, plotOpt);
set(gcf, 'units','inches','PaperPosition', [0 0 12 6]);
print(gcf, 'sqrt_ArrTflmod_sim','-r600','-dpng');

clearvars -except data data_train data_validation data_sim...
    TdegK_ref Ua_ref Rug F Naumann_model Naumann_train...
    sqrt_local_fit sqrt_ArrTflmod_model sqrt_ArrTflmod_train sqrt_ArrTflmod_validation sqrt_ArrTflmod_sim...
    data_DOD

%% Check validity of model simplification
% Plot a comparison of the Naumann model and simplified model predictions:
plotOpt.layout = 'individual axes';
plotOpt.labels = {'ArrTfl_{mod} Fit','Naumann fit'};
plotOpt.colors = {'k','g'};
plotOpt.xlabel = 'EFCs';
plotOpt.ylabel = 'Relative discharge capacity';
plotOpt.title = sprintf("EFCs^{0.5} (ArrTfl_{mod}) MAE=%0.3g%%, Naumann EFCs^{0.5} model MAE=%0.3g%%", sqrt_ArrTflmod_train.MAE*100, Naumann_train.MAE*100);
plotOpt.confidenceInterval = [5 95];
plot_capacity_fits(data_train.t, [sqrt_ArrTflmod_train, Naumann_train], data_train, plotOpt)
set(gcf, 'units','inches','PaperPosition', [0 0 25 30]);
print(gcf, 'sqrt_ArrTflmod_train vs Naumann','-r600','-dpng');

clearvars -except data data_train data_validation data_sim ...
    sqrt_local_fit sqrt_ArrTflmod_model sqrt_ArrTflmod_train sqrt_ArrTflmod_validation sqrt_ArrTflmod_sim...
    Naumann_model Naumann_train data_DOD

%% 双层优化-符号回归幂律模型: q=α1-β1*t^α2
% 全局参数:α1,α2 局部参数β1
% Step.1 双层优化-全局参数优化,优化得到αi和βi
disp(" ")
disp("识别优化幂律模型.")
% Define an equation for bi-level optimization:
% Input vars: t
% Input params: alpha_1, alpha_2, beta_1
power_model_bilevel = @(p_gbl,p_lcl,x) p_gbl(1) - p_lcl(1).*x(:,1).^p_gbl(2);
p0_gbl = [1,0.5]; % α1,α2全局参数
p0_lcl = 0.001; % β1子模型参数
% Assemble the data:
x = data_train.t;
y = data_train.qdis;
cellNums = data_train.cellNum;
% Optimize:
power_model_bilevel_training = optimize_bilevel(x, y, cellNums, power_model_bilevel, p0_gbl, p0_lcl);
fprintf("幂律模型在训练集误差: %0.3g%%\n", power_model_bilevel_training.MAE*100);
% Plot the local fit results:
plotOpt.layout = 'single axis';
plotOpt.labels = {'Fit'};
plotOpt.colors = {'b'};
plotOpt.xlabel = 'Time (days)';
plotOpt.ylabel = 'Relative discharge capacity';
plotOpt.title = sprintf("q_{dis} = \\alpha_1 - \\beta_1*EFCs^{(\\alpha_2)}, Local MAE=%0.3g%%", power_model_bilevel_training.MAE*100);
plotOpt.confidenceInterval = [];
plot_capacity_fits(x(:,1), power_model_bilevel_training, data_train, plotOpt)

% Step.2 Lasso自动识别子模型局部参数β1{γ1,γ2,γ3,γ4}
% Automatically identify beta_1 sub-model:
disp("利用Lasso L1正则化识别子模型局部参数...")
% Assemble data:
% Only one beta_1 value per data series, grab the value of invariant data
% variables for each data series to train with:
submodel_data = assemble_invariant_data(data_train);
% x = [submodel_data.TdegK, submodel_data.soc, submodel_data.U_a, submodel_data.Cavg, submodel_data.dod];
% y = power_model_bilevel_training.p_lcl(:,1);
x = [submodel_data.Cavg, submodel_data.dod];
y = power_model_bilevel_training.p_lcl(:,1);
cellNums = submodel_data.cellNum;

% Step.2.1 生成可能得描述符,特征为T, SOC及U_a, Cavg, DOD
% Generate possible descriptors for linear and multiplicative models:
% x0 = {submodel_data.TdegK, [submodel_data.soc, submodel_data.U_a], submodel_data.Cavg, submodel_data.dod};
% x0_vars = {{'TdegK'}, {'soc', 'U_a'}, {'Cavg'}, {'dod'}}; % XA,XB,XC,XD

x0 = {submodel_data.Cavg, submodel_data.dod};
x0_vars = {{'Cavg'}, {'dod'}}; % XA,XB

[xLin, xLin_vars] = generate_features_linear(x0, x0_vars); % 线性符号库
[xMult, xMult_vars] = generate_features_multiplicative(x0, x0_vars); % 乘法符号库

%% Step.2.2 识别子模型β1参数
% 设置采用SISSO或者Lasso算法来筛选最优的符号算子
% 当符号算子较多100+,而样本数据较少时,应考虑使用SISSO
select_alg = 'SISSO'; % Lasso

if strcmp(select_alg,'SISSO')

    disp("Sub-models for beta_1")
    disp('下面进行SISSO算法优化.')
    % SISSO算法
    % nNonzeroCoefs: 非零系数的最大数量
    % nFeaturesPerSisIter:每次迭代选择的特征数量(每个非零系数一次迭代)
    % nNonzeroCoefs*nFeaturesPerSisIter要小于等于特征数量
    nNonzeroCoefs = 2; nFeaturesPerSisIter = 6;
    fprintf("Searching for models up to %d dimemsions, considering %d new features per iteration.\n",nNonzeroCoefs, nFeaturesPerSisIter)
    siLnFitInfo = SissoRegressor(nNonzeroCoefs, nFeaturesPerSisIter); % 线性
    siLnFitInfo = fitSisso(siLnFitInfo, xLin, y);
    printModels(siLnFitInfo, xLin_vars);

    siMultFitInfo = SissoRegressor(nNonzeroCoefs, nFeaturesPerSisIter); % 乘法
    siMultFitInfo = fitSisso(siMultFitInfo, xMult, log(y)); % 乘法通过取对数转为线性
    printModels(siMultFitInfo, xMult_vars);

    % 根据SiSSO结果构造子模型函数及系数值
    % 线性
    siLn_eq = identify_descriptors_SISSO(xLin_vars, siLnFitInfo, 'lin'); % 传入变量符号及SISSO优化对象
    siLn_model = construct_func_handle(siLn_eq.eq_1SE, [x0_vars{:}]); % 生成模型

    % 乘法
    siMul_eq = identify_descriptors_SISSO(xMult_vars, siMultFitInfo, 'mult'); 
    siMul_model = construct_func_handle(siMul_eq.eq_1SE, [x0_vars{:}]);

    % Reoptimize the models:
    fitOpt.CV = 'LeaveOut'; fitOpt.bootstrap = 'Off';
    siLn_model_train = optimize(x, y, cellNums, siLn_model, siLn_eq.p_1SE, fitOpt);
    siMul_model_train = optimize(x, y, cellNums, siMul_model, siMul_eq.p_1SE, fitOpt);
    fprintf("Linear model R2adj: %0.3g, Multiplicative model R2adj: %0.3g\n", siLn_model_train.R2adj, siMul_model_train.R2adj);

else
    disp("Sub-models for beta_1")
    disp('下面进行Lasso算法优化.')
    CV = 4; % LASSO优化4折交叉验证
    plotOpt.CVplot = 'On'; plotOpt.lambdaplot = 'Off'; % The lambdaplot (see lasso help) is often uninterpretable.
    linearFitInfo = identify_descriptors_linear(xLin, xLin_vars, y, CV, plotOpt);
    linear_model = construct_func_handle(linearFitInfo.eq_1SE, [x0_vars{:}]);
    
    multFitInfo = identify_descriptors_multiplicative(xMult, xMult_vars, y, CV, plotOpt);
    mult_model = construct_func_handle(multFitInfo.eq_1SE, [x0_vars{:}]);

    % Reoptimize the models:
    fitOpt.CV = 'LeaveOut'; fitOpt.bootstrap = 'Off';
    linear_model_train = optimize(x, y, cellNums, linear_model, linearFitInfo.p_1SE, fitOpt);
    mult_model_train = optimize(x, y, cellNums, mult_model, multFitInfo.p_1SE, fitOpt);
    fprintf("Linear model R2adj: %0.3g, Multiplicative model R2adj: %0.3g\n", linear_model_train.R2adj, mult_model_train.R2adj);
end

%%
% Step.2.3 绘图比较三种模型
% Plot a comparison of the results:
figure; hold on; box on; grid on;
if strcmp(select_alg,'SISSO')
    plot(submodel_data.dod.*100, y, 'ok', 'MarkerSize', 6, 'LineWidth', 1.5) % y为β1
    plot(submodel_data.dod.*100, siLn_model_train.y_fit, 'xr', 'MarkerSize', 6, 'LineWidth', 1.5);
    plot(submodel_data.dod.*100, siMul_model_train.y_fit, '+b', 'MarkerSize', 6, 'LineWidth', 1.5);
    xlabel('DOD (%)'); ylabel(sprintf("\\beta_1 (EFCs^{-%0.2g})", power_model_bilevel_training.p_gbl(2)));
    legend('Locally fit values', 'Linear sub-model', 'Multiplicative sub-model', 'Location', 'northwest')
    title(sprintf("Linear sub-model: R^2_{adj}=%0.2g, MAPE=%0.3g%% \nMultiplicative sub-model: R^2_{adj}=%0.2g, MAPE=%0.3g%%",...
        siLn_model_train.R2adj, siLn_model_train.MAPE*100, siMul_model_train.R2adj, siMul_model_train.MAPE*100))

else
    plot(submodel_data.TdegC, y, 'ok', 'MarkerSize', 6, 'LineWidth', 1.5) % y为β1
    plot(submodel_data.TdegC, linear_model_train.y_fit, 'xr', 'MarkerSize', 6, 'LineWidth', 1.5);
    plot(submodel_data.TdegC, mult_model_train.y_fit, '+b', 'MarkerSize', 6, 'LineWidth', 1.5);
    xlabel('Temperature (\circC)'); ylabel(sprintf("\\beta_1 (days^{-%0.2g})", power_model_bilevel_training.p_gbl(2)));
    legend('Locally fit values', 'Linear sub-model', 'Multiplicative sub-model', 'Location', 'northwest')
    title(sprintf("Linear sub-model: R^2_{adj}=%0.2g, MAPE=%0.3g%% \nMultiplicative sub-model: R^2_{adj}=%0.2g, MAPE=%0.3g%%",...
        linear_model_train.R2adj, linear_model_train.MAPE*100, mult_model_train.R2adj, mult_model_train.MAPE*100))
end
%%
% Step.3 使用乘法模型构建全局模型
% Build a global model using the multiplicative model:
% Input vars: t, TdegK, soc, U_a
% Input params: alpha_0, alpha_2, beta_1(gamma_0, gamma_1, gamma_2)

% power_model = @(p,x) p(1) - siMul_model(p(3:5),x(:,[2:3])).*(x(:,1).^p(2));
% p0 = [power_model_bilevel_training.p_gbl, siMul_eq.p_1SE];

% 线性
power_model_ln = @(p,x) p(1) - siLn_model(p(3:5),x(:,[2:3])).*(x(:,1).^p(2));
p0 = [power_model_bilevel_training.p_gbl, siLn_eq.p_1SE];
% Assemble the data:
x = [data_train.t, data_train.Cavg, data_train.dod];
y = data_train.qdis;
cellNums = data_train.cellNum;

% Optimize with cross-validation and bootstrap resampling:
fitOpt.CV = 'LeaveOut'; fitOpt.bootstrap = 'On'; fitOpt.bootstrapIterations = 1000;
disp("自举采样和交叉验证优化幂律模型...")
power_model_train = optimize(x, y, cellNums, power_model, p0, fitOpt);
fprintf("优化后的幂律模型在训练集MAE: %0.3g%%\n", power_model_train.MAE*100);
% Plot global fit results:
plotOpt.layout = 'individual axes';
plotOpt.labels = {'Fit'};
plotOpt.colors = {'b'};
plotOpt.xlabel = 'EFCs';
plotOpt.ylabel = 'Relative capacity';
plotOpt.title = sprintf("q_{dis} = \\alpha_0 - \\beta_1(C_{avg},DOD)*t^{(\\alpha_2)}, Global MAE=%0.3g%%", power_model_train.MAE*100);
plotOpt.confidenceInterval = [5 95];
plot_capacity_fits(x(:,1), power_model_train, data_train, plotOpt)
set(gcf, 'units','inches','PaperPosition', [0 0 25 30]);
print(gcf, 'siMult_train','-r600','-dpng');

%%
% Test model on validation data:
x = [data_validation.t, data_validation.Cavg, data_validation.dod];
y = data_validation.qdis;
power_model_validation = evaluate(x, y, power_model, power_model_train);
fprintf("幂律模型全局模型对验证集MAE: %0.3g%%\n", power_model_validation.MAE*100);
plotOpt.title = sprintf("q_{dis} = \\alpha_0 - \\beta_1(C_{avg},DOD)*t^{(\\alpha_2)}, Global Validation MAE=%0.3g%%", power_model_validation.MAE*100);
plot_capacity_fits(x(:,1), power_model_validation, data_validation, plotOpt)
set(gcf, 'units','inches','PaperPosition', [0 0 25 18]);
print(gcf, 'siMult_validation','-r600','-dpng');

% Simulate 12000 cycs aging:
x = [data_sim.t, data_sim.Cavg, data_sim.dod];
power_model_sim = simulate(x, power_model, power_model_train);
% Plot the simulation result:
x = data_sim.t;
plotOpt.xlabel = 'EFCs';
plotOpt.title =  sprintf("1,2000 cycs simulation, q_{dis} = \\alpha_0 - \\beta_1(C_{avg},DOD)*t^{(\\alpha_2)}");
plot_capacity_sim(x, power_model_sim, data_sim, plotOpt);
set(gcf, 'units','inches','PaperPosition', [0 0 12 6]);
print(gcf, 'siMult_sim','-r600','-dpng');

%%
clearvars -except data data_train data_validation data_sim ...
    sqrt_local_fit sqrt_ArrTflmod_model sqrt_ArrTflmod_train sqrt_ArrTflmod_validation sqrt_ArrTflmod_sim ...
    power_model_bilevel_training power_model power_model_train power_model_validation power_model_sim...
    siLnFitInfo siLn_model siLn_eq siMultFitInfo siMul_model siMul_eq...
    siLn_model_train siMul_model_train

%%
















%% Automatically identify sigmoidal model
disp(" ")
disp("识别优化S型模型.")
% Define an equation for bi-level optimization:
% Input vars: t
% Input params: alpha_0, beta_1, alpha_2, beta_3
sigmoidal_model_bilevel = @(p_gbl,p_lcl,x) p_gbl(1) - 2.*p_lcl(1).*(0.5-(1./(1+exp((p_gbl(2).*x(:,1)).^p_lcl(2)))));
p0_gbl = [1, 1e-2];
p0_lcl = [0.1, 0.5];
% Assemble the data:
x = data_train.t;
y = data_train.qdis;
cellNums = data_train.cellNum;
% Optimize:
sigmoidal_model_bilevel_training = optimize_bilevel(x, y, cellNums, sigmoidal_model_bilevel, p0_gbl, p0_lcl);
fprintf("S型局部模型对训练集MAE: %0.3g%%\n", sigmoidal_model_bilevel_training.MAE*100);
% Plot the local fit results:
plotOpt.layout = 'single axis';
plotOpt.labels = {'Fit'};
plotOpt.colors = {'r'};
plotOpt.xlabel = 'Time (days)';
plotOpt.ylabel = 'Relative discharge capacity';
plotOpt.title = sprintf("q_{dis} = \\alpha_0 - 2*\\beta_1*(0.5-(1/(1+exp((\\alpha_2*t)^{(\\beta_3)})))), Local MAE=%0.3g%%", sigmoidal_model_bilevel_training.MAE*100);
plotOpt.confidenceInterval = [];
plot_capacity_fits(x(:,1), sigmoidal_model_bilevel_training, data_train, plotOpt)


% Automatically identify beta_1 sub-model:
disp("利用Lasso L1正则化识别子模型局部参数...")
disp("Sub-model for beta_1")
% Assemble data:
% Only one beta_1 value per data series, grab the value of invariant data
% variables for each data series to train with:
submodel_data = assemble_invariant_data(data_train);
x = [submodel_data.TdegK, submodel_data.soc, submodel_data.U_a];
y = sigmoidal_model_bilevel_training.p_lcl(:,1);
cellNums = submodel_data.cellNum;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 采用符号描述回归,加法线性和乘法
% Generate possible descriptors for linear and multiplicative models:
x0 = {submodel_data.TdegK, [submodel_data.soc, submodel_data.U_a]};
x0_vars = {{'TdegK'}, {'soc', 'U_a'}};
[xLin, xLin_vars] = generate_features_linear(x0, x0_vars);
[xMult, xMult_vars] = generate_features_multiplicative(x0, x0_vars);
% Identify linear sub-model descriptors
CV = 4; % 4 fold cross-validation for LASSO optimization.
plotOpt.CVplot = 'On'; plotOpt.lambdaplot = 'Off'; % The lambdaplot (see lasso help) is often uninterpretable.
b1_linearFitInfo = identify_descriptors_linear(xLin, xLin_vars, y, CV, plotOpt);
b1_linear_model = construct_func_handle(b1_linearFitInfo.eq_1SE, [x0_vars{:}]);
b1_multFitInfo = identify_descriptors_multiplicative(xMult, xMult_vars, y, CV, plotOpt);
b1_mult_model = construct_func_handle(b1_multFitInfo.eq_1SE, [x0_vars{:}]);

% Reoptimize the models:
fitOpt.CV = 'LeaveOut'; fitOpt.bootstrap = 'Off';
b1_linear_model_train = optimize(x, y, cellNums, b1_linear_model, b1_linearFitInfo.p_1SE, fitOpt);
b1_mult_model_train = optimize(x, y, cellNums, b1_mult_model, b1_multFitInfo.p_1SE, fitOpt);
fprintf("Linear model R2adj: %0.3g, Multiplicative model R2adj: %0.3g\n", b1_linear_model_train.R2adj, b1_mult_model_train.R2adj);
% Plot a comparison of the results:
figure; hold on; box on; grid on;
plot(submodel_data.TdegC, y, 'ok', 'MarkerSize', 6, 'LineWidth', 1.5)
plot(submodel_data.TdegC, b1_linear_model_train.y_fit, 'xr', 'MarkerSize', 6, 'LineWidth', 1.5);
plot(submodel_data.TdegC, b1_mult_model_train.y_fit, '+b', 'MarkerSize', 6, 'LineWidth', 1.5);
xlabel('Temperature (\circC)'); ylabel("\beta_1");
legend('Locally fit values', 'Linear sub-model', 'Multiplicative sub-model', 'Location', 'northwest')
title(sprintf("Linear sub-model: R^2_{adj}=%0.2g, MAPE=%0.3g%% \nMultiplicative sub-model: R^2_{adj}=%0.2g, MAPE=%0.3g%%",...
    b1_linear_model_train.R2adj, b1_linear_model_train.MAPE*100, b1_mult_model_train.R2adj, b1_mult_model_train.MAPE*100))
% Multiplicative model is best.

% Automatically identify beta_3 sub-model:
disp("Sub-model for beta_2")
% Assemble data:
% Only one beta_2 value per data series, grab the value of invariant data
% variables for each data series to train with:
x = [submodel_data.TdegK, submodel_data.soc, submodel_data.U_a];
y = sigmoidal_model_bilevel_training.p_lcl(:,2);
cellNums = submodel_data.cellNum;
% Generate possible descriptors for linear and multiplicative models:
x0 = {submodel_data.TdegK, [submodel_data.soc, submodel_data.U_a]};
x0_vars = {{'TdegK'}, {'soc', 'U_a'}};
[xLin, xLin_vars] = generate_features_linear(x0, x0_vars);
[xMult, xMult_vars] = generate_features_multiplicative(x0, x0_vars);

% Identify linear sub-model descriptors
CV = 4; % 4 fold cross-validation for LASSO optimization.
plotOpt.CVplot = 'On'; plotOpt.lambdaplot = 'Off'; % The lambdaplot (see lasso help) is often uninterpretable.
b3_linearFitInfo = identify_descriptors_linear(xLin, xLin_vars, y, CV, plotOpt);
b3_linear_model = construct_func_handle(b3_linearFitInfo.eq_1SE, [x0_vars{:}]);
b3_multFitInfo = identify_descriptors_multiplicative(xMult, xMult_vars, y, CV, plotOpt);
b3_mult_model = construct_func_handle(b3_multFitInfo.eq_1SE, [x0_vars{:}]);
% Reoptimize the models:
fitOpt.CV = 'LeaveOut'; fitOpt.bootstrap = 'Off';
b3_linear_model_train = optimize(x, y, cellNums, b3_linear_model, b3_linearFitInfo.p_1SE, fitOpt);
b3_mult_model_train = optimize(x, y, cellNums, b3_mult_model, b3_multFitInfo.p_1SE, fitOpt);
fprintf("Linear model R2adj: %0.3g, Multiplicative model R2adj: %0.3g\n", b3_linear_model_train.R2adj, b3_mult_model_train.R2adj);
% Plot a comparison of the results:
figure; hold on; box on; grid on;
plot(submodel_data.TdegC, y, 'ok', 'MarkerSize', 6, 'LineWidth', 1.5)
plot(submodel_data.TdegC, b3_linear_model_train.y_fit, 'xr', 'MarkerSize', 6, 'LineWidth', 1.5);
plot(submodel_data.TdegC, b3_mult_model_train.y_fit, '+b', 'MarkerSize', 6, 'LineWidth', 1.5);
xlabel('Temperature (\circC)'); ylabel("\beta_3)");
legend('Locally fit values', 'Linear sub-model', 'Multiplicative sub-model', 'Location', 'northwest')
title(sprintf("Linear sub-model: R^2_{adj}=%0.2g, MAPE=%0.3g%% \nMultiplicative sub-model: R^2_{adj}=%0.2g, MAPE=%0.3g%%",...
    b3_linear_model_train.R2adj, b3_linear_model_train.MAPE*100, b3_mult_model_train.R2adj, b3_mult_model_train.MAPE*100))


% Build a global model using the multiplicative model:
% Input vars: t, TdegK, soc, U_a
% Input params: alpha_0, alpha_2, beta_1(gamma_0, gamma_1, ...), beta_3(...)
sigmoidal_model = @(p,x) p(1) - 2.*b1_mult_model(p(3:8),x(:,2:4)).*(0.5-(1./(1+exp((p(2).*x(:,1)).^b3_mult_model(p(9:end),x(:,2:4))))));
p0 = [sigmoidal_model_bilevel_training.p_gbl, b1_multFitInfo.p_1SE, b3_multFitInfo.p_1SE];
% Assemble the data:
x = [data_train.t, data_train.TdegK, data_train.soc, data_train.U_a];
y = data_train.qdis;
cellNums = data_train.cellNum;
% Optimize with cross-validation and bootstrap resampling:
fitOpt.CV = 'LeaveOut'; fitOpt.bootstrap = 'On'; fitOpt.bootstrapIterations = 1000;
disp("自举采样和交叉验证优化S型模型...")
disp("注意: 自举采样和交叉验证优化S型模型可能需要几分钟.")
sigmoidal_model_train = optimize(x, y, cellNums, sigmoidal_model, p0, fitOpt);
fprintf("S型模型在训练集上的MAE: %0.3g%%\n", sigmoidal_model_train.MAE*100);

% Plot global fit results:
plotOpt.layout = 'individual axes';
plotOpt.labels = {'Fit'};
plotOpt.colors = {'r'};
plotOpt.xlabel = 'Time (days)';
plotOpt.ylabel = 'Relative discharge capacity';
plotOpt.title = sprintf("q_{dis} = \\alpha_0 - 2*\\beta_1*(0.5-(1/(1+exp((\\alpha_2*t)^{(\\beta_3)})))), Global MAE=%0.3g%%", sigmoidal_model_train.MAE*100);
plotOpt.confidenceInterval = [5 95];
plot_capacity_fits(x(:,1), sigmoidal_model_train, data_train, plotOpt)

% Test model on validation data:
x = [data_validation.t, data_validation.TdegK, data_validation.soc, data_validation.U_a];
y = data_validation.qdis;
sigmoidal_model_validation = evaluate(x, y, sigmoidal_model, sigmoidal_model_train);
fprintf("S型模型在验证集上的MAE: %0.3g%%\n", sigmoidal_model_validation.MAE*100);
plotOpt.title = sprintf("q_{dis} = \\alpha_0 - 2*\\beta_1*(0.5-(1/(1+exp((\\alpha_2*t)^{(\\beta_3)})))), Global Validation MAE=%0.3g%%", sigmoidal_model_validation.MAE*100);
plot_capacity_fits(x(:,1), sigmoidal_model_validation, data_validation, plotOpt)

% Simulate 20 years aging:
x = [data_sim.t, data_sim.TdegK, data_sim.soc, data_sim.U_a];
sigmoidal_model_sim = simulate(x, sigmoidal_model, sigmoidal_model_train);
% Plot the simulation result:
x = data_sim.t_years;
plotOpt.xlabel = 'Time (years)';
plotOpt.title = '20 year simulation, q_{dis} = \alpha_0 - 2*\beta_1*(0.5-(1/(1+exp((\alpha_2*t)^{(\beta_3)}))))';
plot_capacity_sim(x, sigmoidal_model_sim, data_sim, plotOpt);

clearvars -except data data_train data_validation data_sim ...
    sqrt_local_fit sqrt_ArrTflmod_model sqrt_ArrTflmod_train sqrt_ArrTflmod_validation sqrt_ArrTflmod_sim ...
    power_model_bilevel_training power_model power_model_train power_model_validation power_model_sim ...
    sigmoidal_model_bilevel_training sigmoidal_model sigmoidal_model_train sigmoidal_model_validation sigmoidal_model_sim

%% Compare various models
% Compare global model capacity predictions:
plotOpt.labels = {'t^{0.5} (ArrTfl_{mod})', 'Power Law', 'Sigmoidal'};
plotOpt.colors = {'k', 'b', 'r'};
plotOpt.xlabel = 'Time (days)';
plotOpt.ylabel = 'Relative discharge capacity';

% 获取当前文件夹路径
currentFolder = pwd;
% 设置保存图形的次级文件夹名称
subFolderName = 'Output_results';
% 创建次级文件夹（如果尚不存在）
if ~exist(subFolderName, 'dir')
    mkdir(subFolderName);
end
% 构建完整的保存路径
fullSubFolderPath = fullfile(currentFolder, subFolderName);

% Training data:
plotOpt.title = "Comparison of global models - Training data";
fitResults = [sqrt_ArrTflmod_train, power_model_train, sigmoidal_model_train];
% Individual plots with confidence intervals:
plotOpt.layout = 'individual axes'; plotOpt.confidenceInterval = [5 95];
plot_capacity_fits(data_train.t, fitResults, data_train, plotOpt)

% % 保存图片
% set(gcf, 'units','inches','PaperPosition', [0 0 40 50]);
% print(gcf, [fullSubFolderPath,'\Train_data1'],'-r300','-dpng') % 注意更改保存位置
% Single plot w/o confidence intervals (nice for see global residuals):
plotOpt.layout = 'single axis'; plotOpt.confidenceInterval = [];
plot_capacity_fits(data_train.t, fitResults, data_train, plotOpt)

% Validation data:
plotOpt.title = "Comparison of global models - Validation data";
fitResults = [sqrt_ArrTflmod_validation, power_model_validation, sigmoidal_model_validation];
% Individual plots with confidence intervals:
plotOpt.layout = 'individual axes'; plotOpt.confidenceInterval = [5 95];
plot_capacity_fits(data_validation.t, fitResults, data_validation, plotOpt)
% Single plot w/o confidence intervals (nice for see global residuals):
plotOpt.layout = 'single axis'; plotOpt.confidenceInterval = [];
plot_capacity_fits(data_validation.t, fitResults, data_validation, plotOpt)

% 20 year simulation:
plotOpt.title = "Comparison of global models - 20 year simulation";
fitResults = [sqrt_ArrTflmod_sim, power_model_sim, sigmoidal_model_sim];
% Individual plots with confidence intervals:
plotOpt.layout = 'individual axes'; plotOpt.confidenceInterval = [5 95];
plot_capacity_sim(data_sim.t, fitResults, data_sim, plotOpt)
% Single plot w/o confidence intervals (nice for see global residuals):
plotOpt.layout = 'single axis'; plotOpt.confidenceInterval = [];
plot_capacity_sim(data_sim.t, fitResults, data_sim, plotOpt)

% Comparison of fit metrics from training:
fitResults = [sqrt_ArrTflmod_train, power_model_train, sigmoidal_model_train];
labels = categorical(plotOpt.labels);
MSE = [fitResults.MSE];
MSE_CV = [fitResults.MSE_CV];
R2adj = [fitResults.R2adj];
% MSE plot:
figure; hold on; box on; grid on;
plot(labels, MSE, 'ok', 'MarkerSize', 6)
plot(labels, MSE_CV, 'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 6)
ylabel('Mean squared error')
legend('MSE', 'MSE_{CV}', 'Location', 'best')
% R2adj plot:
figure; box on; grid on;
plot(labels, R2adj, 'dr', 'MarkerFaceColor', 'r', 'MarkerSize', 6)
ylabel('Adj. coeff. of determination')

% Clean up:
clearvars -except data data_train data_validation data_sim ...
    sqrt_local_fit sqrt_ArrTflmod_model sqrt_ArrTflmod_train sqrt_ArrTflmod_validation sqrt_ArrTflmod_sim ...
    power_model_bilevel_training power_model power_model_train power_model_validation power_model_sim ...
    sigmoidal_model_bilevel_training sigmoidal_model sigmoidal_model_train sigmoidal_model_validation sigmoidal_model_sim ...
    fullSubFolderPath

% Save all workspace variable to file:
mat_file = fullfile(fullSubFolderPath,"fitted_models.mat");
save(mat_file)
disp('已成功保存数据!')