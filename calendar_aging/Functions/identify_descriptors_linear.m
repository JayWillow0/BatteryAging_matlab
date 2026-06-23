function FitInfo = identify_descriptors_linear(x, xvars, y, CV, plotOpt)
rng('default')
lambda = [0 logspace(-6, 0, 200)]; % 在10e-6-1之间取对数间隔点200个,网格搜索最佳Lasso正则化惩罚系数λ
[BB, FitInfo] = lasso(x, y, 'Lambda', lambda, 'CV', CV, 'PredictorNames', xvars);
if strcmp(plotOpt.CVplot, 'On')
    lassoPlot(BB, FitInfo, 'PlotType', 'CV'); legend('show')
    set(gca, 'YScale', 'log')
    title("Cross-Validated MSE of Lasso Fit for Linear Model")
end
if strcmp(plotOpt.lambdaplot,'On')
    lassoPlot(BB, FitInfo, 'PlotType', 'Lambda');
end
% 1SE 规则（1-Standard Error Rule）用于选择模型复杂性较低但仍具有较好预测性能的 λ 值。
% 具体来说，它是在交叉验证中寻找使均方误差（MSE）接近最优值的最大 λ 值，其规则如下：
% 在交叉验证中，找出使 MSE 最小的 λ 值（LambdaMinMSE）
% 对应的 MSE 和标准误差（SE）可以通过交叉验证获得
% 误差阈值为 MSE_min + SE，即最小 MSE 加上对应的标准误差
% 找到误差小于等于 MSE_min + SE 的所有 λ 值
% 选择其中最大的 λ 值，作为 Lambda1SE
% SE是K个 MSE 值的标准差/sqrt(K)
% Grab results:
BB_1SE = BB(:,FitInfo.Index1SE);
nonzero_params = BB_1SE ~= 0;
p_1SE = [FitInfo.Intercept(FitInfo.Index1SE), BB_1SE(nonzero_params)'];
Xvars_1SE = ['c_0', xvars(nonzero_params)];
% Create equation string and display the chosen descriptors:
disp("Linear equation:")
for i = 1:length(Xvars_1SE)
    if i == 1
        eq_1SE = 'p(1)';
        eq_str = "c_0";
    else
        eq_1SE = strcat(eq_1SE,'+p(',num2str(i),').*',Xvars_1SE{i});
        eq_str = strcat(eq_str, " + c_", num2str(i-1), "*", Xvars_1SE{i});
    end
end
disp(eq_str)
FitInfo.Xvars_1SE = Xvars_1SE;
FitInfo.eq_1SE = eq_1SE;
FitInfo.p_1SE = p_1SE;
end
