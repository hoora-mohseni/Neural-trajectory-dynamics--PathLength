%====================================================================================
%% Plot PathLength for each PC
%====================================================================================

for pc = 1:10
    figure('color', 'w', 'Position', [300, 300, 700, 300]);
    x_s     = 5; % Initialize x-axis offset
    x_total = x_s;
    red_dark   = [1 0 0];
    for i = 1:4

        epochLength = epochWindowCounts(i);
        bins = x_s:x_s + epochLength - 1;
        %=== Shadow bar for last 4 windows of ALL epochs ===
        fill([bins(end-3)-0.5, bins(end)+0.5, bins(end)+0.5, bins(end-3)-0.5], ...
             [0, 0, 100, 100], ...
             [0.85 0.85 0.85], 'EdgeColor', 'none');
        hold on;

        s1 = shadedErrorBar(x_s:x_s+epochWindowCounts(i)-1, mean(pathLengthPCs.(epochs{i})(:,:,pc)), ...
                               std(pathLengthPCs.(epochs{i})(:,:,pc))/sqrt(numSubjects), ...
                               {'-', 'color', [0.3, 0.05, 0.4], 'LineWidth', 2}, 1, 1);
        hold on
    
        s2 = shadedErrorBar(x_s:x_s+epochWindowCounts(i)-1, mean(pathLengthPCs.(epochs{i+4})(:,:,pc)), ...
                               std(pathLengthPCs.(epochs{i+4})(:,:,pc))/sqrt(numSubjects), ...
                               {'-', 'color',  [0.1 ,  0.6   , 0.4], 'LineWidth', 1.5}, 1, 1);
        x_s = x_s + epochWindowCounts(i) + 10; 
        
        x_total = [x_total, x_s-10];
        x_total = [x_total, x_s];
   
    end
    box off
    set(gca, "FontSize", 15)
    %yticks(10:5:30); ylim([10, 30])
    xticks(x_total); xticklabels([1, 21, 22, 51, 52, 131, 132, 161])
    xlabel("Trial Number"); pl = ylabel(sprintf("Path Length, pc %d", pc)); xlim([0, x_total(end-1)])
    pl.Position = pl.Position - [3, 0, 0];
    legend([s1.patch, s2.patch], {"Day 1", "Day 2"}, "Location", "North")
    legend box off
end

%====================================================================================
%%  Early & Late barplots
%====================================================================================

for pc = 1:10
    
    allSubjectData = cell(numel(epochs), 2);  % For scatter
    
    % Compute means and collect per-subject values
    for p = 1:numel(epochs)
       
        myMean(p, 1) = mean(mean(PCspathlength.early.(epochs{p})(:, :, pc), 2));
        myMean(p, 2) = mean(mean(PCspathlength.late.(epochs{p})(:, :, pc), 2));
    
        allSubjectData{p, 1} = mean(PCspathlength.early.(epochs{p})(:, :, pc), 2);
        allSubjectData{p, 2} = mean(PCspathlength.late.(epochs{p})(:, :, pc), 2);
    
        myEpochs{p} = epochs{p}(1:end-4);
    end
    
    % Plot bars
    figure('Color', 'w', 'Position', [300, 300, 800, 400])
    b = bar(myMean, 'grouped'); hold on;
    
    % Color scheme
    blue_light = [0.1, 0.6, 0.4];
    blue_dark  = [0 0 1];
    red_light  = [0.3, 0.05, 0.4];
    red_dark   = [1 0 0];
    
    for k = 1:numel(epochs)
        if k <= 4
            b(1).FaceColor = 'flat';
            b(1).CData(k,:) = red_light;
            b(2).FaceColor = 'flat';
           
            b(2).CData(k,:) = red_light;
            b(2).FaceAlpha = 0.5;
        else
            b(1).CData(k,:) = blue_light;
            b(2).CData(k,:) = blue_light;
            b(2).FaceAlpha = 0.5;
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
    h2 = bar(nan, nan, 'FaceColor', red_light, 'EdgeColor', 'none');
    h3 = bar(nan, nan, 'FaceColor', blue_light, 'EdgeColor', 'none');
    h4 = bar(nan, nan, 'FaceColor', blue_light, 'EdgeColor', 'none');
    
    % legend([h1, h2, h3, h4], ...
    %     {'Early Day1', 'Late Day1', 'Early Day2', 'Late Day2'}, ...
    %     'Location', 'northoutside');
    
    % xlabel('Epoch');
    % ylabel('Path Length');
    xticks(1:numel(epochs))
    xticklabels(myEpochs)
    xlim([0.5, 8 + 0.5])  %numel(epochs)
    yticks(20:10:50); ylim([20,50])
    box off
    set(gca, "FontSize", 20, 'LineWidth', 0.8, 'TickDir', 'non', 'XTickLabel', '')

end

% legend([h1, h2,], {}, "Location", "northeast")
% legend box off

%====================================================================================
%% pathLengh total and each Pc correlation- both Days
%====================================================================================

for p = 1:10
    
    figure('color', 'w', 'Position', [300, 300, 1000, 400]);
    r_path = [];
    r_err = [];
    colorIdx = 1;

    for i = [2:4, 6:8]
        allPath = mean(pathLengthSliding.(epochs{i}));
        err     = mean(pathLengthPCs.(epochs{i})(:,:, p));
    
        % allPath(err>35) = [];
        % err(err>35) = [];
        
        if i == 2 || i ==6
            sign = 'o' ;
        elseif i == 3 || i == 7
            sign = "square";
        else
            sign = '^';
        end
    
        if i == 2 || i == 3 || i == 4
            for j = 1:numel(err)
        
                a(i) = plot(err(j), allPath(j), sign, 'MarkerFaceColor', [0.3, 0.05, 0.4], 'MarkerEdgeColor', [0.3, 0.05, 0.4]); hold on
                colorIdx = colorIdx + 1;
            end
        end
        if i == 6 || i == 7 || i == 8
            for j = 1:numel(err)
        
                a(i) = plot(err(j), allPath(j), sign, 'MarkerFaceColor', [0.1 , 0.6  , 0.4], 'MarkerEdgeColor', [0.1 , 0.6  , 0.4]); hold on
                colorIdx = colorIdx + 1;
            end
        end
    
        r_err =  [r_err, err];
        r_path = [r_path, allPath];
    end
    
    % Remove NaNs before correlation
    validIdx = ~isnan(r_err) & ~isnan(r_path);
    [R, P] = corrcoef(r_err(validIdx), r_path(validIdx));
    R_value = R(1,2);
    P_value = P(1,2);
    
    % Display R and p in command window
    fprintf('R value: %.3f\n', R_value);
    fprintf('p value: %.3f\n', P_value);
    
    % Annotate R and p on the plot
    text(0.05, 0.9, sprintf('R = %.2f\np = %.3f', R_value, P_value ), ...
       'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold');
    
    
     mdl = fitlm(r_err, r_path);
     plot(r_err, mdl.Fitted, 'color', 'k', 'LineWidth', 2);
    
    %yticks(70:10:120)
    %xticks(0:15:45)
    %xlim([0, 45]); %ylim([70, 120])
    ylabel('all PCs Path Length ')
    xlabel(['PC', num2str(p)])
    box off; 
    %c.Ticks = [];        % Removes both ticks and labels
    set(gca, "FontSize", 20, 'LineWidth', .8, 'TickDir', 'out')
    axis square
    %legend([a(6), a(7), a(8)], {' ', ' ', ' '});

end

%====================================================================================
%% Show first 10 PCs variance
%====================================================================================

disp(table((1:10)', l(1:10), ex(1:10), ...
    'VariableNames', {'PC', 'Variance', 'PercentExplained'}));

cumulativeVariance = cumsum(ex);
disp(cumulativeVariance(1:10));

figure('color', 'w', 'Position', [200, 200, 400, 600]);

x = linspace(1,10,10);  % default evenly spaced
x = x * 0.8;

b = bar(x, ex(1:10), 'FaceColor', [0.3, 0.05, 0.4], ...
    'BarWidth', 0.5); 

hold on;
plot(x, cumulativeVariance(1:10), 'k.-', 'LineWidth', 2, ...
    'Marker', '.', 'MarkerSize', 15);
set(gca,'XTick', x, 'XTickLabel', 1:10);
xlabel('PC'); ylabel('Percent Explained');
title('Explained Variance & Cumulative Sum');
legend('Percent Explained','Cumulative Sum','Location','best');
box off;

%====================================================================================
%%  Correlate Each PC with Brain Map 
%====================================================================================


meanActivity = [];
r = [];

for i = 1:numel(epochs)
    phase = epochs{i};
    for subj = 1:numSubjects
        currentMatrix = eval(phase);  % Dimensions: time x regions x subjects
        meanActivity.early.(phase)(subj, :) = mean(currentMatrix(1:40,:,subj), 1); 
        meanActivity.late.(phase)(subj, :) = mean(currentMatrix(end-39:end,:,subj), 1);
        for pc_num = 1:10

            pc_vector = coeff_ref(:,pc_num);
            r.early.(phase)(subj, pc_num) = corr(pc_vector, meanActivity.early.(phase)(subj, :)');  % scalar result
            r.late.(phase)(subj, pc_num) = corr(pc_vector, meanActivity.late.(phase)(subj, :)'); 
        end
    end
end

%% 

for pc  = 1:3
    
    allSubjectData = cell(numel(epochs), 2);  % For scatter
    
    % Compute means and collect per-subject values
    for p = 1:numel(epochs)
       
        myMean(p, 1) = mean(r.early.(epochs{p})(: , pc));
        myMean(p, 2) = mean(r.late.(epochs{p})(: , pc));
    
        allSubjectData{p, 1} = r.early.(epochs{p})(: , pc);
        allSubjectData{p, 2} = r.early.(epochs{p})(: , pc);
    
        myEpochs{p} = epochs{p}(1:end-4);
    end
    
    % Plot bars
    figure('Color', 'w', 'Position', [300, 300, 800, 400])
    b = bar(myMean, 'grouped'); hold on;
    
    % Color scheme
    blue_light = [0.1, 0.6, 0.4];
    blue_dark  = [0 0 1];
    red_light  = [0.3, 0.05, 0.4];
    red_dark   = [1 0 0];
    
    for k = 1:numel(epochs)
        if k <= 4
            b(1).FaceColor = 'flat';
            b(1).CData(k,:) = red_light;
            b(2).FaceColor = 'flat';
           
            b(2).CData(k,:) = red_light;
            b(2).FaceAlpha = 0.5;
        else
            b(1).CData(k,:) = blue_light;
            b(2).CData(k,:) = blue_light;
            b(2).FaceAlpha = 0.5;
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
    h2 = bar(nan, nan, 'FaceColor', red_light, 'EdgeColor', 'none');
    h3 = bar(nan, nan, 'FaceColor', blue_light, 'EdgeColor', 'none');
    h4 = bar(nan, nan, 'FaceColor', blue_light, 'EdgeColor', 'none');
    
    % legend([h1, h2, h3, h4], ...
    %     {'Early Day1', 'Late Day1', 'Early Day2', 'Late Day2'}, ...
    %     'Location', 'northoutside');
    % 
    % xlabel('Epoch');
    % ylabel('Path Length');
    % xticks(1:numel(epochs))
    % xticklabels(myEpochs)
    % xlim([0.5, 8 + 0.5])  %numel(epochs)
    % yticks(20:10:50); ylim([20,50])
    box off
    set(gca, "FontSize", 20, 'LineWidth', 0.8, 'TickDir', 'non', 'XTickLabel', '')

end

% legend([h1, h2,], {}, "Location", "northeast")
% legend box off




