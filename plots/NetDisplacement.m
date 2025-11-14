%====================================================================================
%% --- Plot net displacement without Overlap --- skip first window!
%====================================================================================

figure('color', 'w', 'Position', [300, 300, 1000, 400]);

x_s     = 5; 
x_total = x_s;

for i = 1:4

    epochLength = epochWindowCounts(i);
    bins = x_s:x_s + epochLength - 1;
    % === Shadow bar for last 4 windows of ALL epochs ===
    fill([bins(end-3)-0.5, bins(end)+0.5, bins(end)+0.5, bins(end-3)-0.5], ...
         [0, 0, 200, 200], ...
         [0.85 0.85 0.85], 'EdgeColor', 'none');
    hold on;

   
        fill([bins(1)-0.5, bins(4)+0.5, bins(4)+0.5, bins(1)-0.5], ...
             [0, 0, 200, 200], ...
             [0.85 0.85 0.85], 'EdgeColor', 'none');  % light gray
        hold on
   

    s1 = shadedErrorBar(x_s:x_s+epochWindowCounts(i)-1, mean(distance.(epochs{i})), ...
                           std(distance.(epochs{i}))/sqrt(numSubjects), ...
                           {'-', 'color', [0.3, 0.05, 0.4], 'LineWidth', 2}, 1, 1);
    
    hold on

    s2 = shadedErrorBar(x_s:x_s+epochWindowCounts(i)-1, mean(distance.(epochs{i+4})), ...
                           std(distance.(epochs{i+4}))/sqrt(numSubjects), ...
                           {'-', 'color',[0.1 , 0.6, 0.4], 'LineWidth', 2}, 1, 1);
    
    x_s = x_s + epochWindowCounts(i) + 10; 
    x_total = [x_total, x_s-10, x_s];

end


box off
set(gca, "FontSize", 20, 'LineWidth', .8, 'TickDir', 'out')
%yticks(70:10:120); ylim([70, 120])
xticks(x_total ); xticklabels([1, 21, 22, 51, 52, 131, 132, 161])
xlabel("Bin Number"); pl = ylabel("Net Displacement"); 
xlim([0, x_total(end-1)])
ylim([10, 25])
pl.Position = pl.Position - [3, 0, 0];

legend([s1.patch, s2.patch], {"Day 1", "Day 2"}, "Location", "NorthEast")
legend box off

%====================================================================================
%% Early and Late Barplots
%====================================================================================


figure('color', 'w', 'Position', [300, 300, 1000, 400])
b = bar(myMean, 'grouped'); hold on; 

% Color scheme
blue_light = [0.1, 0.6, 0.4];
red_light  = [0.3, 0.05, 0.4];

for k = 1:numel(epochs)
    b(1).FaceColor = 'flat';
    b(2).FaceColor = 'flat';

    if k <= 4
        
        b(1).CData(k,:) = red_light + .15*(1-red_light);
        b(2).CData(k,:) = red_light + .5*(1-red_light);
    else
        b(1).CData(k,:) = blue_light + .06*(1-blue_light);
        b(2).CData(k,:) = blue_light + .5*(1-blue_light);
    end
end

% Add scatter points (simulating error bars with spread data)
ngroups = numel(epochs);
nbars = 2;
groupwidth = min(0.8, nbars / (nbars + 1.5));
markerSize = 40;

for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    for p = 1:ngroups
        data = allSubjectData{p, i};
        n = numel(data);
        % Slight jitter to spread dots around x-position
        jitterX = x(p) + 0.08 * (rand(1, n) - 0.5);

       
        scatter(jitterX, data, markerSize, 'k', 'filled', 'MarkerFaceAlpha', 0.6)
    end
end


hold on;
h1 = bar(nan, nan, 'FaceColor', red_light, 'EdgeColor', 'none');
h2 = bar(nan, nan, 'FaceColor', red_light + .5*(1-red_light), 'EdgeColor', 'none');
h3 = bar(nan, nan, 'FaceColor', blue_light, 'EdgeColor', 'none');
h4 = bar(nan, nan, 'FaceColor', blue_light + .5*(1-blue_light), 'EdgeColor', 'none');
hold off;

% Create manual legend
legend([h1, h2, h3, h4], ...
       {'Early Day1', 'Late Day1', ...
        'Early Day2', 'Late Day2'}, ...
       'Location', 'northoutside');

xlabel('Epoch');
ylabel('Net Displacement');
xticks(1:numel(epochs))
xticklabels(myEpochs)
xlim([0.5, numel(epochs) + 0.5])
yticks(5:5:30); ylim([5, 30])
box off
set(gca, "FontSize", 20, 'LineWidth', 0.8, 'TickDir', 'non', 'XTickLabel', '')

