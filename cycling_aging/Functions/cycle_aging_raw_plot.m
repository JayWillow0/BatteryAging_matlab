function cycle_aging_raw_plot(x, x_days, y, data, cellNums, s)
figure;
gt = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
% 倍率影响 cap vs days, cellNum 1-6
nexttile
colors1 = colormap(slanCM('Greens'));
for cellNum = unique(cellNums,'stable')'
    if cellNum <= 6
        mask = cellNums == cellNum;
        Ctemp = data.Cavg(mask);
        scatter(x_days(mask), y(mask), s, ...
            'LineWidth', 2, 'MarkerEdgeColor', colors1(round(Ctemp(1)/1.5*256), :), 'MarkerFaceColor', 'None');
    end

    grid on; box on; hold on;
    xlabel('Time (days)');
    ylabel('Relative capacity');
    xlim([0 1000]);
    ylim([0.7 1]);
end
title('C-rate effect')

% 倍率影响 cap vs EFCs, cellNum 1-6
nexttile
colors1 = colormap(slanCM('Greens'));
for cellNum = unique(cellNums,'stable')'
    if cellNum <= 6
        mask = cellNums == cellNum;
        Ctemp = data.Cavg(mask);

        scatter(x(mask), y(mask), s, ...
            'LineWidth', 2, 'MarkerEdgeColor', colors1(round(Ctemp(1)/1.5*256), :), 'MarkerFaceColor','None');

        c1 = colorbar('Ticks', [0.2, 0.5, 1.5], 'Location', 'east', 'Limits', [0.2, 1.5], colormap=colors1);
        c1.Label.String = 'C_{avg} (hr^{-1})';
        c1.Label.FontSize = 8;
        c1.FontSize = 8;
    end

    grid on; box on; hold on;
    xlabel('Equivalent full cycles');
    ylabel('Relative capacity');
    xlim([0 15000]);
    ylim([0.7 1]);
end
title('C-rate effect')

% DOD影响 cap vs EFCs, cellNum 7-17
nexttile
colors2 = colormap(slanCM('RdBu'));
for cellNum = unique(cellNums,'stable')'
    if cellNum <= 17 && cellNum >= 7
        mask = cellNums == cellNum;
        dodtemp = data.dod(mask);

        scatter(x(mask), y(mask), s, ...
            'LineWidth', 2, 'MarkerEdgeColor', colors2(round(dodtemp(1)/1*254), :), 'MarkerFaceColor','None');

        c2 = colorbar('Ticks', [0.05, 0.25, 0.5, 0.75, 1], 'Location', 'east', 'Limits', [0.05, 1],  colormap=colors2);
        c2.Label.String = 'DOD';
        c2.Label.FontSize = 8;
        c2.FontSize = 8;
    end
    
    grid on; box on; hold on;
    xlabel('Equivalent full cycles');
    ylabel('Relative capacity');
    xlim([0 15000]);
    ylim([0.7 1]);
end
title('DOD effect')

% SOC影响 cap vs EFCs, cellNum 18-20
nexttile
colors3 = colormap(slanCM('Copper'));
for cellNum = unique(cellNums,'stable')'
    if cellNum >= 18
        mask = cellNums == cellNum;
        soctemp = data.soc(mask);

        scatter(x(mask), y(mask), s, ...
            'LineWidth', 2, 'MarkerEdgeColor', colors3(round(soctemp(1)/1*256), :), 'MarkerFaceColor','None');

        c3 = colorbar('Ticks', [0.25, 0.5, 0.75], 'Location', 'east', 'Limits', [0.25, 0.75],  colormap=colors3);
        c3.Label.String = 'SOC';
        c3.Label.FontSize = 8;
        c3.FontSize = 8;
    end
    
    grid on; box on; hold on;
    xlabel('Equivalent full cycles');
    ylabel('Relative capacity');
    xlim([0 15000]);
    ylim([0.7 1]);
end
title('SOC effect')

set(gcf, 'units','inches','PaperPosition', [0 0 12 10]);
print(gcf, 'raw_plot','-r600','-dpng') % 注意更改保存位置