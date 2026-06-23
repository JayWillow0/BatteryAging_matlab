function siFitInfo = identify_descriptors_SISSO(xvars, obj, type)

SS_1SE = obj.listOfCoefs{length(obj.selectedIndicesL0)}; % 系数矩阵
nonzero_params = SS_1SE ~= 0; % 非零系数
Xvars_1SE = ['c_0', xvars(nonzero_params)]; % 非零系数的符号

% p_1SE = [obj.intercept, SS_1SE(nonzero_params)]; % 系数
% 
% if length(obj.selectedIndicesCurrent) == 1 % 判断优化符号算子的维度
%     eq_1SE = 'p(1)';
%     eq_str = "c_0";
% else
%     eq_1SE = 'p(1)';
%     eq_str = "c_0";
% 
%     for i=1:length(obj.selectedIndicesCurrent)
%         eq_1SE = strcat(eq_1SE,'+p(',num2str(i+1),').*', xvars{obj.selectedIndicesCurrent(i)});
%         eq_str = strcat(eq_str, " + c_", num2str(i), "*", xvars{obj.selectedIndicesCurrent(i)});
%     end
% end

if strcmp(type, 'lin')
    % 线性
    p_1SE = [obj.intercept, SS_1SE(nonzero_params)]; % 系数
    if length(obj.selectedIndicesCurrent) == 1
        eq_1SE = 'p(1)';
        eq_str = "c_0";
    else
        eq_1SE = 'p(1)';
        eq_str = "c_0";
        for i=1:length(obj.selectedIndicesCurrent)
            eq_1SE = strcat(eq_1SE,'+p(',num2str(i+1),').*', xvars{obj.selectedIndicesCurrent(i)});
            eq_str = strcat(eq_str, " + c_", num2str(i), "*", xvars{obj.selectedIndicesCurrent(i)});
        end
    end

else
    % 乘法
    p_1SE = [exp(obj.intercept), SS_1SE(nonzero_params)]; % 系数
    if length(obj.selectedIndicesCurrent) == 1
        eq_1SE = 'p(1)';
        eq_str = "c_0";
    else
        eq_1SE = 'p(1)';
        eq_str = "c_0";
        for i=1:length(obj.selectedIndicesCurrent)
            eq_1SE = strcat(eq_1SE,'.*exp(p(',num2str(i+1),').*',xvars{obj.selectedIndicesCurrent(i)},')');
            eq_str = strcat(eq_str, "*exp(c_", num2str(i), "*", xvars{obj.selectedIndicesCurrent(i)}, ")");
        end
    end

end

disp(eq_str);
siFitInfo.Xvars_1SE = Xvars_1SE;
siFitInfo.eq_1SE = eq_1SE;
siFitInfo.p_1SE = p_1SE;
end