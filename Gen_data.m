clc 
clear
close

%%

% Cortex data directory
data_dir = [pwd, '/schaefer400_subcortex'];
cerebellum_dir = [pwd, '/nettekoven_cerebellar']; % Cerebellum data directory

files = dir(data_dir);
cerebellum_files = dir(cerebellum_dir);

% Subject list
subjects = table2array(readtable('valid_subjects.csv'));

% Remove spike subjects
subjects([6, 21, 27]) = [];
subjects = array2table(subjects, 'VariableNames', {'sub'});

% Match cortex files with subjects
matchedFiles = {};
for i = 1:numel(files)
    for j = 1:numel(subjects)
        if contains(files(i).name, subjects{j, 1}, 'IgnoreCase', true)
            matchedFiles{end+1, 1} = files(i).name;
            break; 
        end
    end
end

% Match cerebellum files with subjects
matchedCerebellumFiles = {};
for i = 1:numel(cerebellum_files)
    for j = 1:numel(subjects)
        if contains(cerebellum_files(i).name, subjects{j, 1}, 'IgnoreCase', true)
            matchedCerebellumFiles{end+1, 1} = cerebellum_files(i).name;
            break; 
        end
    end
end

% Load data and combine cortex and cerebellum
combinedData = cell(numel(matchedFiles), 1); % Store combined data

for i = 1:numel(matchedFiles)
    cortexData = readmatrix([data_dir, '/', matchedFiles{i}], 'FileType', 'text', 'Delimiter', '\t');
    cerebellumData = readmatrix([cerebellum_dir, '/', matchedCerebellumFiles{i}], 'FileType', 'text', 'Delimiter', '\t');
    combinedData{i} = [cortexData, cerebellumData];
end

headers = readtable([data_dir, '/', matchedFiles{1}], 'FileType', 'text');
headers = headers.Properties.VariableNames;


% Combined data is now ready for further analysis
disp('Data successfully combined for all matched subjects.');


%% Extract Phases

% leftbaseline - rightbaseline - lefttransfer - rightlearning
numSubjects = numel(subjects);

% Separate rest Data
rest = {};

for i = 1:numel(matchedFiles)
    if contains(matchedFiles{i}, 'rest', 'IgnoreCase', true)
      rest{end+1, 1} = i;   
    end
end

Rest = []; 
for idx = 1:numel(rest)
    Rest(:, :, idx) = combinedData{rest{idx}}; 
end

% Separate leftbaseline Data
leftbaseline = {};

for i = 1:numel(matchedFiles)
    if contains(matchedFiles{i}, 'leftbaseline', 'IgnoreCase', true)
      leftbaseline{end+1, 1} = i;   
    end
end

LeftBaseline = []; 
for idx = 1:numel(leftbaseline)
    LeftBaseline(:, :, idx) = combinedData{leftbaseline{idx}}; 
end
LeftBaseline = LeftBaseline(6:end-4, :, :);

% Separate rightbaseline Data
rightbaseline = {};

for i = 1:numel(matchedFiles)
    if contains(matchedFiles{i}, 'rightbaseline', 'IgnoreCase', true)
      rightbaseline{end+1, 1} = i;   
    end
end

RightBaseline = []; 
for idx = 1:numel(rightbaseline)
    RightBaseline(:, :, idx) = combinedData{rightbaseline{idx}}; 
end
RightBaseline = RightBaseline(6:end-4, :, :);

% Separate lefttransfer Data
lefttransfer = {};

for i = 1:numel(matchedFiles)
    if contains(matchedFiles{i}, 'lefttransfer', 'IgnoreCase', true)
      lefttransfer{end+1, 1} = i;   
    end
end

LeftTransfer = []; 
for idx = 1:numel(lefttransfer)
    LeftTransfer(:, :, idx) = combinedData{lefttransfer{idx}}; 
end
LeftTransfer = LeftTransfer(6:end-4, :, :);

% Separate rightlearning Data
rightlearning = {};

for i = 1:numel(matchedFiles)
    if contains(matchedFiles{i}, 'rightlearning', 'IgnoreCase', true)
      rightlearning{end+1, 1} = i;   
    end
end

RightLearning = []; 
for idx = 1:numel(rightlearning)
    RightLearning(:, :, idx) = combinedData{rightlearning{idx}}; 
end
RightLearning = RightLearning(6:end-4, :, :);


phases = {'Rest', 'LeftBaseline', 'RightBaseline', 'RightLearning', 'LeftTransfer'};


%%  Smooth data using Savitzky-Golay filter 

% Initialize structure to store smoothed data 
smoothedData   = struct();

% Define the polynomial order and frame length for the Savitzky-Golay filter
polynomialOrder = 1;
frameLength     = 3;  % You can adjust this depending on your data

% Apply Savitzky-Golay filter and PCA projection for each phase
for i = 1:numel(phases)
    phase = phases{i};
    for subj = 1:numel(subjects)
        currentMatrix = eval(phase); 
        % Smooth the data using the Savitzky-Golay filter
        smoothedData.(phase)(:,:,subj) = sgolayfilt(currentMatrix(:,:,subj), polynomialOrder, frameLength, [], 1);

        % Project the smoothed data into PCA space (using the first few components)
        [coeff, score, ~] = pca(smoothedData.(phase)(:,:,subj));
        pcaProjections.(phase)(:,:,subj) = score(:, 1:10);  % Using the first 10 principal components
        pcaCoeff.(phase)(:, :, subj) = coeff(:, 1:10);
        
    end
end

%% Sliding Window Path Length and Displacement Calculation

W = 24; % Window size
S = 12;  % Step size

% Concatenate all phases into a single time series per subject
pathLengthSliding = struct();
epochWindowCounts = zeros(1, numel(phases)); % Store number of windows for each epoch

% First, preallocate matrices for all phases before looping over subjects
for i = 1:numel(phases)
    % Find the maximum possible number of windows across subjects
    maxEpochLength = size(pcaProjections.(phases{i}), 1);
    maxNumWindows = max(0, floor((maxEpochLength - W) / S) + 1);
    
    % Preallocate with correct dimensions **outside the subject loop**
    pathLengthSliding.(phases{i}) = zeros(numSubjects, maxNumWindows);
end

% Now process each subject and update the values correctly
for subj = 1:numSubjects
    for i = 1:numel(phases)
        epochData = pcaProjections.(phases{i})(:, :, subj);
        epochLength = size(epochData, 1);  % Determine actual time length of epoch

        % Dynamically determine number of windows for the given epoch
        numWindows = max(0, floor((epochLength - W) / S) + 1);
        epochWindowCounts(i) = numWindows; % Store window count
        
        % Compute metrics in sliding window fashion
        for win = 1:numWindows
            startIdx = (win - 1) * S + 1;
            endIdx   = startIdx + W - 1;
            
            if endIdx > epochLength
                continue; % Skip if window exceeds data length
            end
            
            traj       = epochData(startIdx:endIdx, :); % Extract windowed data
            diffs      = diff(traj);
            dists      = sqrt(sum(diffs.^2, 2));
            pathLength = sum(dists); % Sum over distances in window
            
            % Compute netdisplacement
            
            if startIdx == 1
                startIdx = startIdx+3;
            end
            if win == numWindows
                endIdx = endIdx - 4;
            end
            startPos  = epochData(startIdx, :);         
            endPos    = epochData(endIdx, :);  
            netDistance = sqrt(sum((endPos - startPos).^2));

            % Store results without overwriting previous phases
            displacement.(phases{i})(subj, win, :) = netDistance;
            pathLengthSliding.(phases{i})(subj, win) = pathLength;
            
        end
    end
end

%% --- Plot Path Length without Overlap 

figure('color', 'w', 'Position', [300, 300, 1400, 400]);

x_s     = 5; % Initialize x-axis offset
x_total = x_s;
red_dark   = [1 0 0];

for i = 2:5
    epochLength = epochWindowCounts(i);
    bins = x_s:x_s + epochLength - 1;
    
    s1 = shadedErrorBar(x_s:x_s+epochWindowCounts(i)-1, mean(pathLengthSliding.(phases{i})), ...
                           std(pathLengthSliding.(phases{i}))/sqrt(numSubjects), ...
                            {'-', 'color', [0.3, 0.05, 0.4], 'LineWidth', 2}, 1, 1);
    
hold on

    x_s = x_s + epochWindowCounts(i) + 10; 
    
    x_total = [x_total, x_s-10];
    x_total = [x_total, x_s];
   
end


box off
set(gca, "FontSize", 20, 'LineWidth', .8, 'TickDir', 'out')
yticks(130:10:170); ylim([130, 170])
xticks(x_total); xticklabels([0, 15, 16, 31, 32, 71, 72, 91])  %window size 24v
%xlabel("Trial Number"); pl = ylabel("Path Length"); 
xlim([0, x_total(end-1)])
box off

%% --- Plot Diceplacement without Overlap 

figure('color', 'w', 'Position', [300, 300, 1400, 400]);

x_s     = 5; % Initialize x-axis offset
x_total = x_s;
red_dark   = [1 0 0];

for i = 1:5
    epochLength = epochWindowCounts(i);
    bins = x_s:x_s + epochLength - 1;
    % %=== Shadow bar for last 4 windows of ALL phases ===
    % fill([bins(end-3)-0.5, bins(end)+0.5, bins(end)+0.5, bins(end-3)-0.5], ...
    %      [0, 0, 500, 500], ...
    %      [0.85 0.85 0.85], 'EdgeColor', 'none');
    % hold on;
    % 
    % if i ~= 1 && i ~= 5
    %     fill([bins(1)-0.5, bins(4)+0.5, bins(4)+0.5, bins(1)-0.5], ...
    %          [0, 0, 500, 500], ...
    %          [0.85 0.85 0.85], 'EdgeColor', 'none');  % light gray
    %     hold on
    % end

    
    s1 = shadedErrorBar(x_s:x_s+epochWindowCounts(i)-1, mean(displacement.(phases{i})), ...
                           std(displacement.(phases{i}))/sqrt(numSubjects), ...
                            {'-', 'color', [0.3, 0.05, 0.4], 'LineWidth', 2}, 1, 1);
    
hold on

    x_s = x_s + epochWindowCounts(i) + 10; 
    
    x_total = [x_total, x_s-10];
    x_total = [x_total, x_s];
   
end


box off
set(gca, "FontSize", 20, 'LineWidth', .8, 'TickDir', 'out')
%yticks(10:5:20); ylim([130, 170])
xticks(x_total); xticklabels([-23, 0, 1, 15, 16, 31, 32, 71, 72, 91])  %window size 24v
xlabel("Trial Number"); pl = ylabel("Net Displacement"); 
xlim([0, x_total(end-1)])
% pl.Position = pl.Position - [3, 0, 0];
% legend([s1.patch, s2.patch], {"Day 1", "Day 2"}, "Location", "NorthEast")
% legend box off


%% EArly and Late epoches 

Gen_pathLength = struct();
for i=1:numel(phases)
    Gen_pathLength.early.(phases{i}) = pathLengthSliding.(phases{i})(: ,1:4);
    Gen_pathLength.late.(phases{i}) = pathLengthSliding.(phases{i})(: ,end-3:end);
end

Gen_displacement = struct();
for i=1:numel(phases)
    Gen_displacement.early.(phases{i}) = displacement.(phases{i})(: ,1:4);
    Gen_displacement.late.(phases{i}) = displacement.(phases{i})(: ,end-3:end);
end
%% Path Length Table
stages = {'early', 'late'};
T = table((1: numel(subjects))', 'VariableNames', {'Subject'});

for p = 1:numel(phases)
    for i =1:numel(stages)
        columnName = sprintf('%s_%s', stages{i}, phases{p}); % Format column name
        
        % Extract 32×1 data for each subject
        dataColumn = Gen_pathLength.(stages{i}).(phases{p});
        avgDataColumn = mean(dataColumn, 2);
        % Append to table
        T.(columnName) = avgDataColumn;
    end
end

% Save table as CSV file
%writetable(T, 'Gen_pathLength_EarlyLate_Table.csv');


%% displacement Table


stages = {'early', 'late'};
T = table((1: numel(subjects))', 'VariableNames', {'Subject'});

for p = 1:numel(phases)
    for i =1:numel(stages)
        columnName = sprintf('%s_%s', stages{i}, phases{p}); % Format column name
        
        % Extract 32×1 data for each subject
        dataColumn = Gen_displacement.(stages{i}).(phases{p});
        avgDataColumn = mean(dataColumn, 2);
        % Append to table
        T.(columnName) = avgDataColumn;
    end
end

% Save table as CSV file
writetable(T, 'Gen_displacement_EarlyLate_Table.csv');

%%  PathLength Early & Late barplots

allSubjectData = cell(numel(phases), 2);  % For scatter

% Compute means and collect per-subject values
for p = 2:numel(phases)
   
    myMean(p, 1) = mean(mean(Gen_pathLength.early.(phases{p}), 2));
    myMean(p, 2) = mean(mean(Gen_pathLength.late.(phases{p}), 2));

    allSubjectData{p, 1} = mean(Gen_pathLength.early.(phases{p}), 2);
    allSubjectData{p, 2} = mean(Gen_pathLength.late.(phases{p}), 2);

    myphases{p} = phases{p}(1:end-4);
end


% Plot bars
figure('Color', 'w', 'Position', [300, 300, 1000, 400])
b = bar(myMean, 'grouped'); hold on;

% Color scheme
blue_light = [0.1, 0.6, 0.4];
red_light  = [0.3, 0.05, 0.4];


for k = 1:numel(phases)
    b(1).FaceColor = 'flat';
    b(2).FaceColor = 'flat';
        
    b(1).CData(k,:) = red_light + .15*(1-red_light);
    b(2).CData(k,:) = red_light + .5*(1-red_light);
  
end


% Add scatter points (simulating error bars with spread data)
ngroups = numel(phases);
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
        
        if i == 1 && (p == 1)
            continue;
        end
        scatter(jitterX, data, markerSize, 'k', 'filled', 'MarkerFaceAlpha', 0.6)
    end
end

% Legend
h1 = bar(nan, nan, 'FaceColor', red_light, 'EdgeColor', 'none');
h2 = bar(nan, nan, 'FaceColor', red_light, 'EdgeColor', 'none');
% h3 = bar(nan, nan, 'FaceColor', blue_light, 'EdgeColor', 'none');
% h4 = bar(nan, nan, 'FaceColor', blue_light, 'EdgeColor', 'none');

% legend([b(1), b(2)], ...
%     {'Early', 'Late'}, 'Location', 'northeast');

%xlabel('Phases');
%ylabel('Path Length');
%xticks(1:numel(phases))
xticklabels(myphases)
xlim([1.5, 5+ 0.5])  %numel(phases)
yticks(120:20:200); ylim([120, 200])
box off
set(gca, "FontSize", 20, 'LineWidth', 0.8, 'TickDir', 'non', 'XTickLabel', '')

% legend([h1, h2,], {}, "Location", "northeast")
% legend box off


%%  Net displacement Early & Late barplots

allSubjectData = cell(numel(phases), 2);  % For scatter

% Compute means and collect per-subject values
for p = 1:numel(phases)
   
    myMean(p, 1) = mean(mean(Gen_displacement.early.(phases{p}), 2));
    myMean(p, 2) = mean(mean(Gen_displacement.late.(phases{p}), 2));

    allSubjectData{p, 1} = mean(Gen_displacement.early.(phases{p}), 2);
    allSubjectData{p, 2} = mean(Gen_displacement.late.(phases{p}), 2);

    myphases{p} = phases{p}(1:end-4);
end

myMean(1,1) = 0;

% Plot bars
figure('Color', 'w', 'Position', [300, 300, 800, 400])
b = bar(myMean, 'grouped'); hold on;

% Color scheme
blue_light = [0.1, 0.6, 0.4];
red_light  = [0.3, 0.05, 0.4];


for k = 1:numel(phases)
    b(1).FaceColor = 'flat';
    b(2).FaceColor = 'flat';
        
    b(1).CData(k,:) = red_light + .15*(1-red_light);
    b(2).CData(k,:) = red_light + .5*(1-red_light);
  
end


% Add scatter points (simulating error bars with spread data)
ngroups = numel(phases);
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
        
        if i == 1 && (p == 1)
            continue;
        end
        scatter(jitterX, data, markerSize, 'k', 'filled', 'MarkerFaceAlpha', 0.6)
    end
end

% Legend
h1 = bar(nan, nan, 'FaceColor', red_light, 'EdgeColor', 'none');
h2 = bar(nan, nan, 'FaceColor', red_light, 'EdgeColor', 'none');
% h3 = bar(nan, nan, 'FaceColor', blue_light, 'EdgeColor', 'none');
% h4 = bar(nan, nan, 'FaceColor', blue_light, 'EdgeColor', 'none');

legend([b(1), b(2)], ...
    {'Early', 'Late'}, 'Location', 'northeast');

%xlabel('Phases');
ylabel('Net Displacement');
xticks(1:numel(phases))
xticklabels(myphases)
xlim([0.5, 5 + 0.5])  %numel(phases)
%yticks(5:20:220); 
ylim([5, 25])
box off
set(gca, "FontSize", 20, 'LineWidth', 0.8, 'TickDir', 'non', 'XTickLabel', '')

% legend([h1, h2,], {}, "Location", "northeast")
% legend box off

%% find outlier subjects!

% for subj = 1:numel(subjects)
% 
%     maxVals = max(pathLengthSliding.(phases{5}),[], 2);
%     spike_subjects = find(maxVals > 230);
% 
%     plot(pathLengthSliding.(phases{5})(subj, :));
%     hold on
% end
% display(spike_subjects);

%% behavioral Data

mydir = '/Users/hooramohseni/Desktop/RL-Gen_fMRIData';
files = dir(mydir);
 
subj_behavior = readtable([mydir,'/', 'subject_behavior'], 'TextType', 'string');
 
%%

for i = 1:numel(subjects)

    subject_data = subj_behavior(strcmp(subj_behavior.sub, subjects.sub(i)), :);

    right = subject_data(strcmp(subject_data.hand, 'Right'), :);
    ErRightBaseline(i, :) = (right(right.cursorRotation == 0, :).hitAngle_hand_good(:))';
    ErRightLearning(i, :) = (right(right.cursorRotation == 45 & right.blockNo == 3, :).hitAngle_hand_good(:) + 45)';

    left = subject_data(strcmp(subject_data.hand, 'Left'), :);
    ErLeftBaseline(i, :)  = (left(left.cursorRotation == 0, :).hitAngle_hand_good(:))';
    ErLeftLearning(i, :)  = (left(left.cursorRotation == 45, :).hitAngle_hand_good(:) + 45)';

end
       
errors = {'ErLeftBaseline', 'ErRightBaseline', 'ErRightLearning', 'ErLeftLearning'}; 

W_er = 8;
S_er = 4;

error = [];

for r = 1:numel(errors)
    error_traj  = eval(errors{r}); 
    numWindows  = floor((size(error_traj, 2) - W_er) / S_er) + 1;
    for subj = 1:numel(subjects)
        for win = 1:numWindows
            startIdx = (win - 1) * S_er + 1;
            endIdx   = startIdx + W_er - 1; 
          
            error.(phases{r+1})(subj, win)  = abs(nanmean(error_traj(subj, startIdx:endIdx))); 
        end  
    end
end

%% Plot error

figure('color', 'w', 'Position', [300, 300, 1400, 400]);

x_s     = 5; % Initialize x-axis offset
x_total = x_s;

error.(phases{1}) = nan(35, 23);

for j = 2:5
    
    numWin1 = size(error.(phases{j}),2) ;
    

    s1 = shadedErrorBar(x_s:x_s+ numWin1 - 1, mean(error.(phases{j})), std(error.(phases{j}))/sqrt(32), {'.-', 'color',[0.3, 0.05, 0.4], 'LineWidth', 2}, 1, 1); 
    hold on;
    x_s = x_s + numWin1 + 15;

    x_total = [x_total, x_s-17];
    x_total = [x_total, x_s];

end 

box off
set(gca, "FontSize", 20, 'LineWidth', .8, 'TickDir', 'out')
% yticks(-15:15:45); ylim([0, 45])

xticks(x_total); 
xticklabels([0, 15, 16, 31, 32, 71, 72, 91])  %window size 24v
%xlabel("Trial Number"); pl = ylabel("Net Displacement"); 
xlim([0, x_total(end-1)]);
%xlabel("Trial Number"); pl = ylabel("Angular Error"); 
%pl.Position = pl.Position - [3, 0, 0];


%% pathLengh and error correlation 

figure('color', 'w', 'Position', [300, 300, 1000, 400]);

myColor = viridis(137);
allclr = [];
colorIdx = 30;
r_path = [];
r_err = [];

for i = 2:5
    allPath = mean(pathLengthSliding.(phases{i}));
    err     = mean(error.(phases{i}));

    allPath(err>35) = [];
    err(err>35) = [];

    if i == 2 || (i == 3)
        sign = "square";
    else 
        sign = 'o';
    end

    for j = 1:numel(err)
       clr = myColor(colorIdx, :);
        allclr = [allclr; clr];
        a(i) = plot(err(j), allPath(j), sign, 'MarkerFaceColor', clr, 'MarkerEdgeColor', clr); hold on
        colorIdx = colorIdx + 1;
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
% text(0.05, 0.9, sprintf('R = %.2f\np = %.3f', R_value, P_value), ...
%     'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold');

mdl = fitlm(r_err, r_path);
plot(r_err, mdl.Fitted, 'color', myColor(20, :), 'LineWidth', 2);

% %yticks(70:10:120)
% xticks(0:15:45)
% xlim([0, 45]); %ylim([70, 120])
box off; 
colormap(allclr);
c = colorbar;
c.Ticks = [];        % Removes both ticks and labels
set(gca, "FontSize", 20, 'LineWidth', .8, 'TickDir', 'out')
axis square
%legend([a(2), a(3), a(4)], {' ', ' ', ' '});
