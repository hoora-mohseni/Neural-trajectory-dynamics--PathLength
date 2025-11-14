clc 
clear
close

%% Extract Phases
% Cortex data directory
data_dir = [pwd, '/../VMRD1D2_fMRIData 2/vmr_schaefer400_subcortex'];
cerebellum_dir = [pwd, '/../VMRD1D2_fMRIData 2/nettekoven_cerebellum']; % Cerebellum data directory

files = dir(data_dir);
cerebellum_files = dir(cerebellum_dir);

% Subject list
subject_dir = [pwd, '/../VMRD1D2_fMRIData 2/Data4paper_Aug2020_Bins1-3 (1)'];
subjectsRange = 'C1:C32';
subjects = readcell(subject_dir, 'Range', subjectsRange);

% Match cortex files with subjects
matchedFiles = {};
for i = 1:numel(files)
    for j = 1:numel(subjects)
        if contains(files(i).name, subjects{j}, 'IgnoreCase', true)
            matchedFiles{end+1, 1} = files(i).name;
            break; 
        end
    end
end

% Match cerebellum files with subjects
matchedCerebellumFiles = {};
for i = 1:numel(cerebellum_files)
    for j = 1:numel(subjects)
        if contains(cerebellum_files(i).name, subjects{j}, 'IgnoreCase', true)
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

disp('Data successfully combined for all matched subjects.');

% Days definition
days = {'data_ses1', 'data_ses2'};

% Separate Rotation Data
rotation_D1 = {};
rotation_D2 = {};
for i = 1:numel(matchedFiles)
    if contains(matchedFiles{i}, 'rotation', 'IgnoreCase', true)
        if contains(matchedFiles{i}, 'ses-01', 'IgnoreCase', true)
            rotation_D1{end+1, 1} = i;
        else
            rotation_D2{end+1, 1} = i; 
        end
    end
end

% Load Rotation Data
numSubjects = numel(subjects);
rotationDay1 = []; 
rotationDay2 = []; 
for idx = 1:numel(rotation_D1)
    rotationDay1(:, :, idx) = combinedData{rotation_D1{idx}}; 
end

for idx = 1:numel(rotation_D2)
    rotationDay2(:, :, idx) = combinedData{rotation_D2{idx}}; 
end

% Rotaion Analysis
WL = 240;


baselineDay1 = rotationDay1(9:248, :, :);
learningDay1 = rotationDay1(248+1:248+640, :, :);

baselineDay2 = rotationDay2(9:248, :, :);
learningDay2 = rotationDay2(248+1:248+640, :, :);

% Separate Rest Data
rest_D1 = {};
rest_D2 = {};
for i = 1:numel(matchedFiles)
    if contains(matchedFiles{i}, 'rest_run-2', 'IgnoreCase', true)
        if contains(matchedFiles{i}, 'ses-01', 'IgnoreCase', true)
            rest_D1{end+1, 1} = i; 
        else
            rest_D2{end+1, 1} = i; 
        end
    end
end

% Load Rest Data
restDay1 = []; 
restDay2 = []; 
for idx = 1:numel(rest_D1)
    restDay1(:, :, idx) = combinedData{rest_D1{idx}}; 
end

for idx = 1:numel(rest_D2)
    restDay2(:, :, idx) = combinedData{rest_D2{idx}}; 
end

restDay1 = restDay1(:, :, :);
restDay2 = restDay2(:, :, :);

% Separate Washout Data
washout_D1 = {};
washout_D2 = {};
for i = 1:numel(matchedFiles)
    if contains(matchedFiles{i}, 'washout', 'IgnoreCase', true)
        if contains(matchedFiles{i}, 'ses-01', 'IgnoreCase', true)
            washout_D1{end+1, 1} = i; 
        else
            washout_D2{end+1, 1} = i; 
        end
    end
end

% Load Washout Data
washoutDay1 = [];
washoutDay2 = []; 
for idx = 1:numel(washout_D1)
    washoutDay1(:, :, idx) = combinedData{washout_D1{idx}}; 
end

for idx = 1:numel(washout_D2)
    washoutDay2(:, :, idx) = combinedData{washout_D2{idx}};
end

washoutDay1 = washoutDay1(9:248, :, :);
washoutDay2 = washoutDay2(9:248, :, :);

epochs = {'restDay1', 'baselineDay1', 'learningDay1', 'washoutDay1', 'restDay2', 'baselineDay2', 'learningDay2', 'washoutDay2'};

%%  Smooth data using Savitzky-Golay filter 

smoothedData   = struct();
polynomialOrder = 1;
frameLength     = 3;  

% Apply Savitzky-Golay filter and PCA projection for each phase
for i = 2%1:numel(epochs)
    phase = epochs{i};
    for subj = 1:numSubjects
        currentMatrix = eval(phase); 
        smoothedData.(phase)(:,:,subj) = sgolayfilt(currentMatrix(:,:,subj), polynomialOrder, frameLength, [], 1);

        % Project the smoothed data into PCA space
        [coeff, score, latent, tsquared, explained] = pca(smoothedData.(phase)(:,:,subj));
        latentAllSubjects(1:length(latent), subj) = latent;
        explainedAllSubjects(1:length(explained), subj) = explained;
        pcaProjections.(phase)(:,:,subj) = score(:, 1:10);  % Using the first 10 principal components
        pcaCoeff.(phase)(:, :, subj) = coeff(:, 1:10);
        
    end
end

%%  Smooth data using Savitzky-Golay filter using Baseline Day 1 as refrences

% Initialize structure to store smoothed data 
smoothedData   = struct();

% Define the polynomial order and frame length for the Savitzky-Golay filter
polynomialOrder = 1;
frameLength     = 3;  % You can adjust this depending on your data

% Apply Savitzky-Golay filter and PCA projection for each phase
for i = 1:numel(epochs)
    phase = epochs{i};
    for subj = 1:numSubjects
        currentMatrix = eval(phase); 
        % Smooth the data using the Savitzky-Golay filter
        smoothedData.(phase)(:,:,subj) = sgolayfilt(currentMatrix(:,:,subj), polynomialOrder, frameLength, [], 1);
        
    end
end

% ---------- Build reference PCA on Baseline Day 1 (smoothed) ----------
Xref = [];               % time-stacked baseline data across subjects
for s = 1:numSubjects
    B = smoothedData.baselineDay1(:,:,s);   % [T x R]
    Xref = [Xref; B];                        % stack along time
end

% One PCA, shared by all subjects/phases
[coeff_ref, ~, latent_ref, ~, explained_ref, mu_ref] = pca(Xref); % centers columns by default
coeff_ref = coeff_ref(:,1:10);     % [R x 10] fixed spatial maps
mu_ref    = mu_ref(:)';            % [1 x R] column means used by PCA


%% Smooth data using Simple Moving Average filter

smaSmoothedData = struct();
smaPcaProjections = struct();
smaPcaCoeff = struct();
windowSize = 3;  

for i = 1:numel(epochs)
    phase = epochs{i};
    for subj = 1:numSubjects
        currentMatrix = eval(phase);  
        smoothedMatrix = zeros(size(currentMatrix(:,:,subj)));

        % Apply moving average along time dimension 
        for region = 1:size(currentMatrix,2)
            smoothedMatrix(:, region) = movmean(currentMatrix(:, region, subj), windowSize, 'Endpoints', 'shrink');
        end

        smoothedData.(phase)(:,:,subj) = smoothedMatrix;

        % Project into PCA space 
        [coeff, score, latent, tsquared, explained] = pca(smoothedMatrix);
        pcaProjections.(phase)(:,:,subj) = score(:, 1:10);
        pcaCoeff.(phase)(:,:,subj) = coeff(:, 1:10);
    end
end


%% Sliding Window Path Length

W = 16; % Window size
S = 8;  % Step size

% Concatenate all epochs into a single time series per subject
pathLengthSliding = struct();
nonlinearity = struct();
epochWindowCounts = zeros(1, numel(epochs)); 

% preallocate matrices for all epochs
for i = 1:numel(epochs)
    % the maximum possible number of windows across subjects
    maxEpochLength = size(pcaProjections.(epochs{i}), 1);
    maxNumWindows = max(0, floor((maxEpochLength - W) / S) + 1);
    
    % Preallocate with correct dimensions **outside the subject loop**
    pathLengthSliding.(epochs{i}) = zeros(numSubjects, maxNumWindows);
    nonlinearity.(epochs{i}) = zeros(numSubjects, maxNumWindows);
end

% process each subject and update the values 
for subj = 1:numSubjects
    for i = 1:numel(epochs)
        epochData = pcaProjections.(epochs{i})(:, :, subj);
        epochLength = size(epochData, 1); 

        % Dynamically determine number of windows for the given epoch
        numWindows = max(0, floor((epochLength - W) / S) + 1);
        epochWindowCounts(i) = numWindows; 
        
        % Compute metrics in sliding window 
        for win = 1:numWindows
            startIdx = (win - 1) * S + 1;
            endIdx   = startIdx + W - 1;
            
            if endIdx > epochLength
                continue; 
            end
            
            traj       = epochData(startIdx:endIdx, :); 
            diffs      = diff(traj);
            dists      = sqrt(sum(diffs.^2, 2));
            pathLength = sum(dists); % Sum over distances in window
            
            % Store results without overwriting previous epochs
            pathLengthSliding.(epochs{i})(subj, win) = pathLength;

        end
    end
end

%% --- Plot Path Length with days overlap 

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


%% --- Plot Path Length without Overlap ---

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


%% behavioral Data
mydir = '/Users/hooramohseni/Desktop/Manifold/data';
files = dir(mydir);
 
subj_behavior = readtable([mydir,'/', 'subject_behavior'], 'TextType', 'string');
 
subject_dir = [pwd, '/../VMRD1D2_fMRIData 2/Data4paper_Aug2020_Bins1-3 (1)'];
subjectsRange = 'C1:C32';
subjects = readcell(subject_dir, 'Range', subjectsRange);
 
%%  Extract Errors 


for i = 1:numel(subjects)
    
        subject_data = subj_behavior(strcmp(subj_behavior.Number, subjects{i}), :);

        learning = subject_data(strcmp(subject_data.Block, 'Rotation'), :);
        errorlearningD1(i, :) = learning(strcmp(learning.Day, 'Day 1'), :).Error(:);
        errorlearningD2(i, :) = learning(strcmp(learning.Day, 'Day 2'), :).Error(:);
      
       
        baseline = subject_data(strcmp(subject_data.Block, 'Baseline'), :);
        errorBaseD1(i, :) = baseline(strcmp(baseline.Day, 'Day 1'), :).Error(:);
        errorBaseD2(i, :) = baseline(strcmp(baseline.Day, 'Day 2'), :).Error(:);
       
        washout = subject_data(strcmp(subject_data.Block, 'Washout'), :);
        errorWashoutD1(i, :) = washout(strcmp(washout.Day, 'Day 1'), :).Error(:);
        errorWashoutD2(i, :) = washout(strcmp(washout.Day, 'Day 2'), :).Error(:);
   

end
       
errors = {'errorBaseD1', 'errorlearningD1','errorWashoutD1', 'errorBaseD2','errorlearningD2', 'errorWashoutD2'}; 

W_er = 8;
S_er = 4;

error = [];

for r = 1:numel(errors)
    error_traj  = eval(errors{r})*(180/pi); 
    numWindows  = floor((size(error_traj, 2) - W_er) / S_er) + 1;
    for subj = 1:numSubjects
        for win = 1:numWindows
            startIdx = (win - 1) * S_er + 1;
            endIdx   = startIdx + W_er - 1; 
            if r <= 3
            error.(epochs{r+1})(subj, win)  = abs(mean(error_traj(subj, startIdx:endIdx))); 
            else
            error.(epochs{r+2})(subj, win)  = abs(mean(error_traj(subj, startIdx:endIdx)));  
            end
        end
    end
end

%% Plot error

figure('color', 'w', 'Position', [300, 300, 1000, 400]);

x_s     = 5; %
x_total = x_s;


for j = 1:3
    
    numWin1 = size(error.(epochs{j+1}),2) ;

    s1 = shadedErrorBar(x_s:x_s + numWin1 - 1, mean(error.(epochs{j+1}), 'omitnan'), std(error.(epochs{j+1}), 'omitnan')/sqrt(32), {'color',[0.3, 0.05, 0.4], 'LineWidth', 2}, 1, 1); 
    hold on;

    s2 = shadedErrorBar(x_s:x_s + numWin1 - 1, mean(error.(epochs{j+5}), 'omitnan'), std(error.(epochs{j+5}), 'omitnan')/sqrt(32), {'color',[0.1, 0.6, 0.4], 'LineWidth', 2}, 1, 1); 

    x_s = x_s + numWin1 + 15;

    x_total = [x_total, x_s-17];
    x_total = [x_total, x_s];
    
end 
hold on;

box off;
set(gca, "FontSize", 20, 'LineWidth', .8, 'TickDir', 'out')
yticks(0:15:45); ylim([0, 45])
xticks(x_total(1:end-1)); xticklabels([1, 29, 30, 109, 110, 139])
xlim([0, x_total(end-1)])
xlabel("Trial Number"); pl = ylabel("Angular Error"); 
pl.Position = pl.Position - [3, 0, 0];
lgd = legend([s1.patch, s2.patch], {"Day 1", "Day 2"}, "Location", "NorthEast");
legend box off;
lgd.ItemTokenSize = [40,50];
lgd.IconColumnWidth = 50;


%%  Extract Reaction Times  


for i = 1:numel(subjects)
    
        subject_data = subj_behavior(strcmp(subj_behavior.Number, subjects{i}), :);

        learning = subject_data(strcmp(subject_data.Block, 'Rotation'), :);
        RTlearningD1(i, :) = learning(strcmp(learning.Day, 'Day 1'), :).RT(:);
        RTlearningD2(i, :) = learning(strcmp(learning.Day, 'Day 2'), :).RT(:);
      
       
        baseline = subject_data(strcmp(subject_data.Block, 'Baseline'), :);
        RTBaseD1(i, :) = baseline(strcmp(baseline.Day, 'Day 1'), :).RT(:);
        RTBaseD2(i, :) = baseline(strcmp(baseline.Day, 'Day 2'), :).RT(:);
       
        washout = subject_data(strcmp(subject_data.Block, 'Washout'), :);
        RTWashoutD1(i, :) = washout(strcmp(washout.Day, 'Day 1'), :).RT(:);
        RTWashoutD2(i, :) = washout(strcmp(washout.Day, 'Day 2'), :).RT(:);
   

end
       
RTs = {'RTBaseD1', 'RTlearningD1','RTWashoutD1', 'RTBaseD2','RTlearningD2', 'RTWashoutD2'}; 

W_er = 8;
S_er = 4;

for r = 1:numel(RTs)
    RT_traj  = eval(RTs{r}); 
    numWindows  = floor((size(RT_traj, 2) - W_er) / S_er) + 1;
    for subj = 1:numSubjects
        for win = 1:numWindows
            startIdx = (win - 1) * S_er + 1;
            endIdx   = startIdx + W_er - 1; 
            if r <= 3
            RT.(epochs{r+1})(subj, win)  = nanmean(RT_traj(subj, startIdx:endIdx)); 
            else
            RT.(epochs{r+2})(subj, win)  = nanmean(RT_traj(subj, startIdx:endIdx));  
            end
        end
    end
end

%% Plot reaction time

figure('color', 'w', 'Position', [300, 300, 1000, 400]);

x_s     = 5; 
x_total = x_s;


for j = 1:3
    
    numWin1 = size(RT.(epochs{j+1}),2) ;
    

    s1 = shadedErrorBar(x_s:x_s + numWin1 - 1, mean(RT.(epochs{j+1})), std(RT.(epochs{j+1}))/sqrt(32), {'.-', 'color',[0.3, 0.05, 0.4], 'LineWidth', 2}, 1, 1); 
    hold on;
    s2 = shadedErrorBar(x_s:x_s + numWin1 - 1, mean(RT.(epochs{j+5})), std(RT.(epochs{j+5}))/sqrt(32), {'-', 'color',[0.1, 0.6, 0.4], 'LineWidth', 2}, 1, 1); 

    x_s = x_s + numWin1 + 15;

    x_total = [x_total, x_s-17];
    x_total = [x_total, x_s];
    
end 
hold on;
box off;
set(gca, "FontSize", 20, 'LineWidth', .8, 'TickDir', 'out')
yticks(0:0.2:1.2); ylim([.6, 1.2])
xticks(x_total(1:end-1)); xticklabels([22, 51, 52, 131, 132, 161])
xlim([0, x_total(end-1)])
xlabel("Trial Number"); pl = ylabel("Reaction Time"); 
pl.Position = pl.Position - [3, 0, 0];
lgd = legend([s1.patch, s2.patch], {"Day 1", "Day 2"}, "Location", "NorthEast");
legend box off
lgd.ItemTokenSize = [40,50];
lgd.IconColumnWidth = 50;

%% EArly and Late epoches 

pathLength = struct();
for i=1:numel(epochs)
    pathLength.early.(epochs{i}) = pathLengthSliding.(epochs{i})(: ,1:4);
    pathLength.late.(epochs{i}) = pathLengthSliding.(epochs{i})(: ,end-3:end);
end

%%
stages = {'early', 'late'};
T = table(subjects, 'VariableNames', {'Subject'});

for p = 1:numel(epochs)
    for i =1:numel(stages)
        columnName = sprintf('%s_%s', stages{i}, epochs{p}); % Format column name
        
        % Extract 32×1 data for each subject
        dataColumn = pathLength.(stages{i}).(epochs{p});
        avgDataColumn = mean(dataColumn, 2);
        % Append to table
        T.(columnName) = avgDataColumn;
    end
end

% Save table as CSV file
%writetable(T, 'Pathlength_EarlyLate_Table.csv');

%%  Early & Late barplots

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


%% Visulize the temporal trajectory 

% Define a subject 
figure;
subj = 9; % Change this to the desired subject index
subjectCentroid = struct();
a = 1:30;
projData = pcaProjections.(epochs{2}); % baseline
% nexttile; 
hold on;

% Get the colormap
nPoints = 30;
cmap = viridis(nPoints+1);

x = projData(:, 1, subj);
y = projData(:, 2, subj);
z = projData(:, 3, subj);

dx = diff(x);
dy = diff(y);
dz = diff(z);

% Plot each segment of the line with a different color
for i = 1:nPoints
    plot3(x(i:i+1), y(i:i+1), z(i:i+1), 'Color', cmap(i,:), 'LineWidth', 2);
    plot3(x(i), y(i), z(i), '.', 'MarkerSize', 15, 'MarkerEdgeColor', cmap(i, :));
    
end

plot3(x(i+1), y(i+1), z(i+1), '.', 'MarkerSize', 15, 'MarkerEdgeColor', cmap(i+1, :));

subjectCentroid = mean(projData(:, :, subj));
scatter3(subjectCentroid(1), subjectCentroid(2), subjectCentroid(3), 100, 'k', 'LineWidth', 2, 'Marker', 'x');
    
colormap(viridis);
title([epochs(2), 'Trajectory and Centroid for Subject ', num2str(subj)]);
grid on;
view(3); axis square; hold off; 
set(gcf, 'Position', [100, 100, 600, 600]); legend('Timecourse', 'Centroid', 'Location', 'bestoutside');

xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
set(gca, 'LineWidth', .7)
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'ZTickLabel', []);


%% pathLengh and error correlation- Day 1

figure('color', 'w', 'Position', [300, 300, 1000, 400]);

myColor = viridis(137);
allclr = [];
colorIdx = 1;
r_path = [];
r_err = [];

for i = 2:4
    allPath = mean(pathLengthSliding.(epochs{i}));
    err     = nanmean(error.(epochs{i}));

    if i == 2 
        sign = 'o' ;
    elseif i == 3 
        sign = "square";
    else
        sign = '^';
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
text(0.05, 0.9, sprintf('R = %.2f\np = %.3f', R_value, P_value), ...
    'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold');

mdl = fitlm(r_err, r_path);
plot(r_err, mdl.Fitted, 'color', myColor(20, :), 'LineWidth', 2);

yticks(70:10:120)
xticks(0:15:45)
xlim([0, 45]); ylim([70, 120])
box off; 
colormap(allclr);
c = colorbar;
c.Ticks = [];       
set(gca, "FontSize", 20, 'LineWidth', .8, 'TickDir', 'out');
axis square;
lgd = legend([a(2), a(3), a(4)], {'Baseline', 'Learning', 'Washout'});
lgd.FontSize = 12;

%% pathLengh and error correlation- both Days

figure('color', 'w', 'Position', [300, 300, 1000, 400]);

r_path = [];
r_err = [];

for i = [2:4, 6:8]
    allPath = mean(pathLengthSliding.(epochs{i}));
    err     = nanmean(error.(epochs{i}));

    allPath(err>35) = [];
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
text(0.05, 0.9, sprintf('R = %.2f\np = %.3f', R_value, P_value ), ...
    'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold');

mdl = fitlm(r_err, r_path);
plot(r_err, mdl.Fitted, 'color', 'k', 'LineWidth', 2);

yticks(70:10:120)
xticks(0:15:45)
xlim([0, 45]); ylim([70, 120])
box off; 
set(gca, "FontSize", 20, 'LineWidth', .8, 'TickDir', 'out')
axis square
lgd = legend([a(2), a(3), a(4)], {'Baseline', 'Learning', 'Washout'});
lgd.FontSize = 12;

%% pathLengh and reaction time correlation- Day 1

figure('color', 'w', 'Position', [300, 300, 1000, 400]);

myColor = viridis(137);
allclr = [];
colorIdx = 1;
r_path = [];
r_reactT = [];

for i = 2:4 
    allPath = mean(pathLengthSliding.(epochs{i}));
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
        a(i) = plot(reactT(j), allPath(j), sign, 'MarkerFaceColor', clr, 'MarkerEdgeColor', clr); hold on
        colorIdx = colorIdx + 1;
    end
    r_reactT =  [r_reactT, reactT];
    r_path = [r_path, allPath];
end

% Remove NaNs before correlation
validIdx = ~isnan(r_reactT) & ~isnan(r_path);
[R, P] = corrcoef(r_reactT(validIdx), r_path(validIdx));
R_value = R(1,2);
P_value = P(1,2);

% Display R and p in command window
fprintf('R value: %.3f\n', R_value);
fprintf('p value: %.3f\n', P_value);
text(0.05, 0.9, sprintf('R = %.2f\np = %.3f', R_value, P_value), ...
    'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold');

mdl = fitlm(r_reactT, r_path);
plot(r_reactT, mdl.Fitted, 'color', myColor(20, :), 'LineWidth', 2);

box off; 
colormap(allclr);
c = colorbar;
c.Ticks = [];        % Removes both ticks and labels
set(gca, "FontSize", 20, 'LineWidth', .8, 'TickDir', 'out')
axis square
lgd = legend([a(2), a(3), a(4)], {'Baseline', 'Learning', 'Washout'});
lgd.FontSize = 12;

%% Error and reaction time correlation- Day 1


figure('color', 'w', 'Position', [300, 300, 1000, 400]);

myColor = viridis(137);
allclr = [];
colorIdx = 1;
r_err = [];
r_reactT = [];

for i = [2:4]
    err = mean(error.(epochs{i}));
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
        a(i) = plot(reactT(j), err(j), sign, 'MarkerFaceColor', clr, 'MarkerEdgeColor', clr); hold on
        colorIdx = colorIdx + 1;
    end

    r_reactT =  [r_reactT, reactT];
    r_err = [r_err, err];
end

% Remove NaNs before correlation
validIdx = ~isnan(r_reactT) & ~isnan(r_err);
[R, P] = corrcoef(r_reactT(validIdx), r_err(validIdx));
R_value = R(1,2);
P_value = P(1,2);

% Display R and p in command window
fprintf('R value: %.3f\n', R_value);
fprintf('p value: %.3f\n', P_value);
text(0.05, 0.9, sprintf('R = %.2f\np = %.3f', R_value, P_value), ...
    'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold');

mdl = fitlm(r_reactT, r_err);
plot(r_reactT, mdl.Fitted, 'color', myColor(20, :), 'LineWidth', 2);

box off; 
colormap(allclr);
c = colorbar;
c.Ticks = [];        % Removes both ticks and labels
set(gca, "FontSize", 20, 'LineWidth', .8, 'TickDir', 'out');
axis square;
lgd = legend([a(2), a(3), a(4)], {'Baseline', 'Learning', 'Washout'});
lgd.FontSize = 12;

%%  multiplt regrassion reaction time & error vith PathLength
% Remove NaNs
validIdx = ~isnan(r_err) & ~isnan(r_path) & ~isnan(r_reactT);

r_pathZ = zscore(r_path(validIdx));
r_reactTZ = zscore(r_reactT(validIdx));
r_errZ = zscore(r_err(validIdx));

% Build table
tbl = table(r_errZ', r_reactTZ', r_pathZ','VariableNames', {'Error', 'RT', 'PathLength'});

% Fit model
mdl = fitlm(tbl, 'PathLength ~ Error + RT');

% Display summary
disp(mdl)

%% plot multiplt regrassion

% Fit the model
%mdl = fitlm([Error, ReT], PathLength);

% Create grid for surface
[eGrid, rtGrid] = meshgrid(linspace(min(r_errZ), max(r_errZ), 30), ...
                           linspace(min(r_reactTZ), max(r_reactTZ), 30));

% Predict over the grid
predictedPath = mdl.Coefficients.Estimate(1) + ...
                mdl.Coefficients.Estimate(2) * eGrid + ...
                mdl.Coefficients.Estimate(3) * rtGrid;

% 3D plot
figure('Color', 'w');
hold on
surf(eGrid, rtGrid, predictedPath, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

% Scatter actual data
scatter3(r_errZ, r_reactTZ, r_pathZ, 60, r_pathZ, 'filled');

xlabel('Error');
ylabel('Reaction Time');
zlabel('Path Length');
title('Multiple Linear Regression: PathLength ~ Error + RT');
colorbar;
grid on;
view(135, 25);  % Adjust the 3D view angle
set(gca, 'FontSize', 14);
axis tight

%%  Correlate Each PC with Brain Map 

meanActivity = [];
r = [];

for i = 1:numel(epochs)
    phase = epochs{i};
    for subj = 1:numSubjects
        currentMatrix = eval(phase);  % Dimensions: time x regions x subjects
        meanActivity.early.(phase)(subj, :) = mean(currentMatrix(1:40,:,subj), 1); 
        meanActivity.late.(phase)(subj, :) = mean(currentMatrix(end-39:end,:,subj), 1);
        for pc_num = 1:10

            pc_vector = pcaCoeff.(phase)(:, pc_num, subj);
            r.early.(phase)(subj, pc_num) = corr(pc_vector, meanActivity.early.(phase)(subj, :)');  % scalar result
            r.late.(phase)(subj, pc_num) = corr(pc_vector, meanActivity.late.(phase)(subj, :)'); 
        end
    end
end

%% 

for pc  = 1:10
    
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

%%

a = [nanmean(errorlearningD1(:,1:16), 2), nanmean(errorlearningD2(:,1:16), 2)]*(180/pi);

figure('Color', 'w', 'Position', [300, 300, 200, 400])
x = [0.2 1.2];
b = bar(x, mean(a), 'grouped','LineWidth',0.8, 'BarWidth', 0.6); hold on;
y = b.XEndPoints;
errorbar(y, mean(a), std(a)/sqrt(20), 'k', 'LineStyle','none', 'LineWidth',0.8, 'CapSize',6)
% % Color scheme
blue_light = [0.1, 0.6, 0.4];
red_light  = [0.3, 0.05, 0.4];
% 
% b(1).FaceColor = red_light + .15*(1-red_light);
% hold on; 
% b(2).FaceColor = blue_light + .06*(1-blue_light);
b.FaceColor = 'flat';

% Assign colors: one row per bar
b.CData = [red_light + .15*(1-red_light);   % red
           blue_light + .06*(1-blue_light)]; % blue
yticks(0:15:45); ylim([0, 45]);
xlim([-.7, 2])
set(gca, 'XTick', [], "FontSize", 20, 'LineWidth', .8, 'TickDir', 'out');
box off
