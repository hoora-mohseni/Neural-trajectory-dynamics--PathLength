%====================================================================================
%% --- Plot Path Length with days overlap 
%====================================================================================

figure('color', 'w', 'Position', [300, 300, 1400, 400]);

x_s     = 5; 
x_total = x_s;
red_dark   = [1 0 0];

for i = 1:4
    epochLength = epochWindowCounts(i);
    bins = x_s:x_s + epochLength - 1;
    %=== Shadow bar for last 4 windows of ALL epochs ===
    fill([bins(end-3)-0.5, bins(end)+0.5, bins(end)+0.5, bins(end-3)-0.5], ...
         [0, 0, 1000, 1000], ...
         [0.85 0.85 0.85], 'EdgeColor', 'none'); % light gray
    hold on;
   
    fill([bins(1)-0.5, bins(4)+0.5, bins(4)+0.5, bins(1)-0.5], ...
         [0, 0, 1000, 1000], ...
         [0.85 0.85 0.85], 'EdgeColor', 'none');  
    hold on


    
    s1 = shadedErrorBar(x_s:x_s+epochWindowCounts(i)-1, mean(pathLengthSliding.(epochs{i})), ...
                           std(pathLengthSliding.(epochs{i}))/sqrt(numSubjects), ...
                            {'-', 'color', [0.3, 0.05, 0.4], 'LineWidth', 2}, 1, 1);
    
    hold on
    
    s2 = shadedErrorBar(x_s:x_s+epochWindowCounts(i)-1, mean(pathLengthSliding.(epochs{i+4})), ...
                           std(pathLengthSliding.(epochs{i+4}))/sqrt(numSubjects), ...
                           {'-', 'color',  [0.1 ,  0.6   , 0.4], 'LineWidth', 1.5}, 1, 1);

    x_s = x_s + epochWindowCounts(i) + 10; 
    
    x_total = [x_total, x_s-10];
    x_total = [x_total, x_s];
   
end


box off
set(gca, "FontSize", 20, 'LineWidth', .8, 'TickDir', 'out')
%yticks(340:20:440); ylim([340, 450])
yticks(70:10:120); ylim([70, 120])
%xticks(x_total); xticklabels([1, 28, 29, 68, 67, 172, 173, 212]) %12 windows!
xticks(x_total); xticklabels([1, 21, 22, 51, 52, 131, 132, 161]) %16 windows!
%xticks(x_total); xticklabels([1, 16, 17, 40, 41, 104, 105, 128]) %20 windows!
%xticks(x_total); xticklabels([1, 10, 11, 25, 26, 65, 66, 80]) %32 windows!
%xticks(x_total); xticklabels([1, 6, 7, 16, 17, 42, 43, 52]) %48Wins
%xticks(x_total); xticklabels([1, 4, 5,11, 12, 31, 32, 41]) %64wins
xlabel("Trial Number"); pl = ylabel("Path Length"); 
xlim([0, x_total(end-1)])
pl.Position = pl.Position - [3, 0, 0];
legend([s1.patch, s2.patch], {"Day 1", "Day 2"}, "Location", "NorthEast")
legend box off

%====================================================================================
%% --- Plot Path Length without Overlap ---
%====================================================================================

figure('color', 'w', 'Position', [300, 300, 1400, 400]);

x_s     = 5; 
x_total = x_s;
red_dark   = [1 0 0];

for i = 1:8
        
    epochLength = epochWindowCounts(i);
    bins = x_s:x_s + epochLength - 1;
   
    if i <= 4
        epochClr = [0.3, 0.05, 0.4];
    else
        epochClr = [0.1,  0.6, 0.4];
    end

    
    s1 = shadedErrorBar(x_s:x_s+epochWindowCounts(i)-1, mean(pathLengthSliding.(epochs{i})), ...
                           std(pathLengthSliding.(epochs{i}))/sqrt(numSubjects), ...
                            {'-', 'color', epochClr, 'LineWidth', 1.5}, 1, 1);
    hold on;

    %angular error without overlap
    % s1 = shadedErrorBar(x_s:x_s + epochWindowCounts(i) - 1, mean(error.(epochs{i}), 'omitnan'), std(error.(epochs{i}), 'omitnan')/sqrt(32), {'-', 'color', epochClr, 'LineWidth', 1.5}, 1, 1); 
    % 
    % hold on
    
    x_s = x_s + epochWindowCounts(i) + 15; 
    x_total = [x_total, x_s-17];
    x_total = [x_total, x_s];
 
end


box off
set(gca, "FontSize", 20, 'LineWidth', .8, 'TickDir', 'out')
yticks(70:10:120); ylim([70, 120])
%yticks(0:15:45); ylim([0, 45]) %angular error without overlap
xticks(x_total); xticklabels([1, 21, 22, 51, 52, 131, 132, 161, 1, 21, 22, 51, 52, 131, 132, 161])
xlabel("Trial Number"); pl = ylabel("Path Length"); 
xlim([0, x_total(end-1)])


%====================================================================================
%%  Early & Late barplots
%====================================================================================

allSubjectData = cell(numel(epochs), 2);  % For scatter

% Compute means and collect per-subject values
for p = 1:numel(epochs)
   
    myMean(p, 1) = mean(mean(pathLength.early.(epochs{p}), 2));
    myMean(p, 2) = mean(mean(pathLength.late.(epochs{p}), 2));

    allSubjectData{p, 1} = mean(pathLength.early.(epochs{p}), 2);
    allSubjectData{p, 2} = mean(pathLength.late.(epochs{p}), 2);

    myEpochs{p} = epochs{p}(1:end-4);
end

% Plot bars
figure('Color', 'w', 'Position', [300, 300, 1000, 400])
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
% 
% allSubjectData{1,1} = zeros(32,1)
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

% Legend
h1 = bar(nan, nan, 'FaceColor', red_light, 'EdgeColor', 'none');
h2 = bar(nan, nan, 'FaceColor', red_light + .5*(1-red_light), 'EdgeColor', 'none');
h3 = bar(nan, nan, 'FaceColor', blue_light, 'EdgeColor', 'none');
h4 = bar(nan, nan, 'FaceColor', blue_light + .5*(1-blue_light), 'EdgeColor', 'none');

legend([h1, h2, h3, h4], ...
    {'Early Day1', 'Late Day1', 'Early Day2', 'Late Day2'}, ...
    'Location', 'northoutside');

% xlabel('Epoch');
% ylabel('Path Length');
% xticks(1:numel(epochs))
%xticklabels(myEpochs)
xlim([0.5, 8 + 0.5])  %numel(epochs)
yticks(60:20:140); ylim([50, 140])
box off
set(gca, "FontSize", 20, 'LineWidth', 0.8, 'TickDir', 'non', 'XTickLabel', '')
legend([h1, h2, h3, h4], {}, "Location", "northeast")
legend box off


