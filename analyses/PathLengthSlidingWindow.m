clc 
clear
close
%====================================================================================
%% Extract Phases
%====================================================================================
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

%====================================================================================
%%  Smooth data using Savitzky-Golay filter 
%====================================================================================

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

%====================================================================================
%% Sliding Window Path Length
%====================================================================================

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

%====================================================================================
%% behavioral Data
%====================================================================================

mydir = '/Users/hooramohseni/Desktop/Manifold/data';
files = dir(mydir);
 
subj_behavior = readtable([mydir,'/', 'subject_behavior'], 'TextType', 'string');
 
subject_dir = [pwd, '/../VMRD1D2_fMRIData 2/Data4paper_Aug2020_Bins1-3 (1)'];
subjectsRange = 'C1:C32';
subjects = readcell(subject_dir, 'Range', subjectsRange);

%==================================================================================== 
%%  Extract Errors 
%====================================================================================

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

%====================================================================================
%%  Extract Reaction Times  
%====================================================================================

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

%====================================================================================
%% Extract Early and Late epochs 
%====================================================================================
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


