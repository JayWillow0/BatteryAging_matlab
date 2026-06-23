% create by liuyang
% 08/06/2024
% 提取原始循环老化数据
% 保存为table
clear; clc; close all;

%%
data = table();
files = dir('*.mat'); % 获取所有.mat文件
cell_num = 1;

for k = 1:length(files)
    file = files(k).name;
    data_temp = load(file);

    if k == 1
        Tem = 40*ones(length(data_temp.Legend_Vec),1);
        SOC = 0.5*ones(length(data_temp.Legend_Vec),1);
        DOD = 0.8*ones(length(data_temp.Legend_Vec),1);
        C_avg = [1;0.5;0.2;0.75;0.75;1.5];
    elseif k == 2
        Tem = 40*ones(length(data_temp.Legend_Vec),1);
        SOC = 0.5*ones(length(data_temp.Legend_Vec),1);
        DOD = [1;1;0.8;0.4;0.2;0.1;0.05];
        C_avg = 1*ones(length(data_temp.Legend_Vec),1);
    elseif k == 3
        Tem = [40;40;25;25];
        SOC = 0.5*ones(length(data_temp.Legend_Vec),1);
        DOD = [1;0.8;0.8;1];
        C_avg = 1*ones(length(data_temp.Legend_Vec),1);
    else
        Tem = 40*ones(length(data_temp.Legend_Vec),1);
        SOC = [0.5;0.75;0.25];
        DOD = 0.2*ones(length(data_temp.Legend_Vec),1);
        C_avg = 1*ones(length(data_temp.Legend_Vec),1);
    end


    for i = 1:length(data_temp.Legend_Vec)
        
        temp = table();
        t_data = data_temp.X_Axis_Data_Mat(:,i);
        q_data = data_temp.Y_Axis_Data_Mat(:,i);

        t_notnan = t_data(~isnan(t_data)); % 防止有NaN数据
        q_notnan = q_data(~isnan(q_data));
        temp.cellNum = cell_num.*ones(length(t_notnan),1);
        % temp.EFCs = data_temp.X_Axis_Data_Mat(:,i); % 等效全生命周期
        temp.t = t_notnan; % 等效全生命周期EFCs
        temp.tdays = t_notnan./(12/C_avg(i));
        temp.TdegC = Tem(i).*ones(length(t_notnan),1);
        temp.TdegK = (Tem(i)+273.15).*ones(length(t_notnan),1);
        temp.soc = SOC(i).*ones(length(t_notnan),1);
        temp.dod = DOD(i).*ones(length(t_notnan),1);
        temp.Cavg = C_avg(i).*ones(length(t_notnan),1);
        temp.qdis = q_notnan;
        temp.qloss = 1 - q_notnan;

        data = [data; temp];
        cell_num = cell_num+1;
    end

end

%%
save('../Naumann.mat', "data"); % 保存至上一级目录


