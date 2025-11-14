
%% net displacement 

W = 16; % Window size
S = 8;  % Step size

distance = struct();
epochWindowCounts = zeros(1, numel(epochs)); 

for i = 1:numel(epochs)
   
    maxEpochLength = size(pcaProjections.(epochs{i}), 1);
    maxNumWindows = max(0, floor((maxEpochLength - W) / S) + 1);
    
    distancePCs.(epochs{i}) = zeros(numSubjects, maxNumWindows);
    distance.(epochs{i}) = zeros(numSubjects, maxNumWindows);
end

for subj = 1:numSubjects
    for i = 1:numel(epochs)
        epochData = pcaProjections.(epochs{i})(:, :, subj);
        epochLength = size(epochData, 1);  

        numWindows = max(0, floor((epochLength - W) / S) + 1);
        epochWindowCounts(i) = numWindows; 
        
        for win = 1:numWindows
            
            startIdx = (win - 1) * S + 1;
            endIdx   = startIdx + W - 1;
            
            if endIdx > epochLength
                continue; % Skip if window exceeds data length
            end
            if startIdx == 1
                startIdx = startIdx+3;
            end
            if win == numWindows
                endIdx = endIdx - 4;
            end
            startPos  = epochData(startIdx, :);         
            endPos    = epochData(endIdx, :);  
            netDistance = sqrt(sum((endPos - startPos).^2));
            distance.(epochs{i})(subj, win, :) = netDistance;

            for pc = 1:10
                startPos  = epochData(startIdx, pc);         
                endPos    = epochData(endIdx, pc);  
                netDistance = sqrt(sum((endPos - startPos).^2));
                distancePCs.(epochs{i})(subj, win, pc) = netDistance;
            end
        end
    end
end

%% --- Plot net displacement without Overlap --- skip first window!

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

%% EArly and Late epoches 

displacement = struct();
for i=1:numel(epochs)
    displacement.early.(epochs{i}) = distance.(epochs{i})(: ,1:4);
    displacement.late.(epochs{i}) = distance.(epochs{i})(: ,end-3:end);
end
%%
stages = {'early', 'late'};
T = table(subjects, 'VariableNames', {'Subject'});

for p = 1:numel(epochs)
    for i =1:numel(stages)
        columnName = sprintf('%s_%s', stages{i}, epochs{p}); % Format column name
        
        % Extract 32×1 data for each subject
        dataColumn = displacement.(stages{i}).(epochs{p});
        avgDataColumn = mean(dataColumn, 2);
        % Append to table
        T.(columnName) = avgDataColumn;
    end
end

% Save table as CSV file
%writetable(T, 'displacement_EarlyLate_Table.csv');

%% PCs - EArly and Late epoches

PCsDisplacement = struct();
for i = 1:numel(epochs)
        PCsDisplacement.early.(epochs{i})= distancePCs.(epochs{i})(:,1:4,:);
        PCsDisplacement.late.(epochs{i})= distancePCs.(epochs{i})(:,end-3:end,:);
end

%%  Convert data to table 

stages = {'early', 'late'};
T = table(subjects, 'VariableNames', {'Subject'});

for pc = 1:10
    for p = 1:numel(epochs)
        for i =1:numel(stages)
            columnName = sprintf('PC%d_%s_%s', pc, stages{i}, epochs{p}); % Format column name
            
            % Extract 32×1 data for each subject
            dataColumn = PCsDisplacement.(stages{i}).(epochs{p})(:, :, pc);
            avgDataColumn = mean(dataColumn, 2);
            % Append to table
            T.(columnName) = avgDataColumn;
        end
    end
end

% Save table as CSV file
%writetable(T, 'PCsDisplacement_earlyLateTable.csv')

%% net displacement BarPlots

allSubjectData = cell(numel(epochs), 2);  % For scatter

for p = 1:numel(epochs)
    myMean(p, 1) = mean(mean(displacement.early.(epochs{p}), 2));
    myMean(p, 2) = mean(mean(displacement.late.(epochs{p}), 2));

    allSubjectData{p, 1} = mean(displacement.early.(epochs{p}), 2);
    allSubjectData{p, 2} = mean(displacement.late.(epochs{p}), 2);

    myEpochs{p} = epochs{p}(1:end-4);
end


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



%% displacement and error correlation- both Days

figure('color', 'w', 'Position', [300, 300, 1000, 400]);


r_distance = [];
r_err = [];

for i = [2:4, 6:8]
    alldistance = mean(distance.(epochs{i}));
    err     = mean(error.(epochs{i}));

    alldistance(err>35) = [];
    err(err>35) = [];
    
    if i == 2 || i ==6
        sign = 'o' ;
    elseif i == 3 || i == 7
        sign = "square";
    else
        sign = '^';
    end

    if i == 2 || i == 3 || i == 4
        for j = 1:numel(err)
    
            a(i) = plot(err(j), alldistance(j), sign, 'MarkerFaceColor', [0.3, 0.05, 0.4], 'MarkerEdgeColor', [0.3, 0.05, 0.4]); hold on
            colorIdx = colorIdx + 1;
        end
    end
    if i == 6 || i == 7 || i == 8
        for j = 1:numel(err)
    
            a(i) = plot(err(j), alldistance(j), sign, 'MarkerFaceColor', [0.1 , 0.6  , 0.4], 'MarkerEdgeColor', [0.1 , 0.6  , 0.4]); hold on
            colorIdx = colorIdx + 1;
        end
    end

    r_err =  [r_err, err];
    r_distance = [r_distance, alldistance];
end

% Remove NaNs before correlation
validIdx = ~isnan(r_err) & ~isnan(r_distance);
[R, P] = corrcoef(r_err(validIdx), r_distance(validIdx));
R_value = R(1,2);
P_value = P(1,2);

% Display R and p in command window
fprintf('R value: %.3f\n', R_value);
fprintf('p value: %.3f\n', P_value);

% Annotate R and p on the plot
text(0.05, 0.9, sprintf('R = %.2f\np = %.3f', R_value, P_value), ...
    'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold');


 mdl = fitlm(r_err, r_distance);
 plot(r_err, mdl.Fitted, 'color', 'k', 'LineWidth', 2);

%yticks(70:10:120)
%xticks(0:15:45)
%xlim([0, 45]); ylim([70, 120])
box off; 
%c.Ticks = [];        % Removes both ticks and labels
set(gca, "FontSize", 20, 'LineWidth', .8, 'TickDir', 'out')
axis square
lgd = legend([a(2), a(3), a(4)], {'Baseline', 'Learning', 'Washout'});
lgd.FontSize = 12;

%% net displacement and reaction time correlation- Day 1

figure('color', 'w', 'Position', [300, 300, 1000, 400]);

myColor = viridis(137);
allclr = [];
colorIdx = 1;
r_distance = [];
r_reactT = [];

for i = 2:4 
    alldistance = mean(distance.(epochs{i}));
    reactT     = mean(RT.(epochs{i}));

    if i == 2 
        sign = 'o' ;
    elseif i == 3 
        sign = "square";
    else
        sign = '^';
    end

    for j = 1:numel(reactT)
       clr = myColor(colorIdx, :);
        allclr = [allclr; clr];
        a(i) = plot(reactT(j), alldistance(j), sign, 'MarkerFaceColor', clr, 'MarkerEdgeColor', clr); hold on
        colorIdx = colorIdx + 1;
    end
    r_reactT =  [r_reactT, reactT];
    r_distance = [r_distance, alldistance];
end

% Remove NaNs before correlation
validIdx = ~isnan(r_reactT) & ~isnan(r_distance);
[R, P] = corrcoef(r_reactT(validIdx), r_distance(validIdx));
R_value = R(1,2);
P_value = P(1,2);

% Display R and p in command window
fprintf('R value: %.3f\n', R_value);
fprintf('p value: %.3f\n', P_value);

% Annotate R and p on the plot
text(0.05, 0.9, sprintf('R = %.2f\np = %.3f', R_value, P_value), ...
    'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold');



mdl = fitlm(r_reactT, r_distance);
plot(r_reactT, mdl.Fitted, 'color', myColor(20, :), 'LineWidth', 2);

%yticks(70:10:120)
%xticks(0:15:45)
%xlim([0, 45]); %ylim([70, 120])
box off; 
colormap(allclr);
c = colorbar;
c.Ticks = [];        % Removes both ticks and labels
set(gca, "FontSize", 20, 'LineWidth', .8, 'TickDir', 'out')
axis square
%legend([a(2), a(3), a(4)], {' ', ' ', ' '});
