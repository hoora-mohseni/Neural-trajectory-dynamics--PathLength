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

matchedFiles = {};
for i = 1:numel(files)
    for j = 1:numel(subjects)
        if contains(files(i).name, subjects{j}, 'IgnoreCase', true)
            matchedFiles{end+1, 1} = files(i).name;
            break; 
        end
    end
end

matchedCerebellumFiles = {};
for i = 1:numel(cerebellum_files)
    for j = 1:numel(subjects)
        if contains(cerebellum_files(i).name, subjects{j}, 'IgnoreCase', true)
            matchedCerebellumFiles{end+1, 1} = cerebellum_files(i).name;
            break; 
        end
    end
end

combinedData = cell(numel(matchedFiles), 1); 

for i = 1:numel(matchedFiles)
    cortexData = readmatrix([data_dir, '/', matchedFiles{i}], 'FileType', 'text', 'Delimiter', '\t');
    cerebellumData = readmatrix([cerebellum_dir, '/', matchedCerebellumFiles{i}], 'FileType', 'text', 'Delimiter', '\t');
    combinedData{i} = [cortexData, cerebellumData];
end

headers = readtable([data_dir, '/', matchedFiles{1}], 'FileType', 'text');
headers = headers.Properties.VariableNames;

disp('Data successfully combined for all matched subjects.');

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
%

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

        %Project the smoothed data into PCA space (using the first few components)
        [coeff, score, latent, ~, explained, mu] = pca(smoothedData.(phase)(:,:,subj));
        pcaProjections.(phase)(:,:,subj) = score(:, 1:10);  % Using the first 10 principal components
        pcaCoeff.(phase)(:, :, subj) = coeff(:, 1:10);
        
    end
end

%====================================================================================
%% Sliding Window Path Length and net displacement Calculation for each PC
%====================================================================================

W = 16; % Window size
S = 8;  % Step size

% Concatenate all epochs into a single time series per subject
pathLengthPCs = struct();
nonlinearity = struct();
epochWindowCounts = zeros(1, numel(epochs)); 

for i = 1:numel(epochs)
   
    maxEpochLength = size(pcaProjections.(epochs{i}), 1);
    maxNumWindows = max(0, floor((maxEpochLength - W) / S) + 1);
    
    pathLengthPCs.(epochs{i}) = zeros(numSubjects, maxNumWindows);
    nonlinearity.(epochs{i}) = zeros(numSubjects, maxNumWindows);
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
            for pc = 1:10
                traj       = epochData(startIdx:endIdx, pc); % Extract windowed data
                diffs      = diff(traj);
                dists      = sqrt(sum(diffs.^2, 2));
                pathLength = sum(dists); % Sum over distances in window
                netDistance = sqrt(sum((endPos - startPos).^2));
                
                % Store results without overwriting previous epochs
                pathLengthPCs.(epochs{i})(subj, win, pc) = pathLength;
                netdisPlacementPCs.(epochs{i})(subj, win, pc) = netDistance;
            end
        end
    end
end

%====================================================================================
%% generating Early and Late windows
%====================================================================================
PCspathlength = struct();
PCSdisplacement = struct();
for i = 1:numel(epochs)
        PCspathlength.early.(epochs{i})= pathLengthPCs.(epochs{i})(:,1:4,:);
        PCspathlength.late.(epochs{i})= pathLengthPCs.(epochs{i})(:,end-3:end,:);
        PCSdisplacement.early.(epochs{i}) = netdisPlacementPCs.(epochs{i})(: ,1:4, :);
        PCSdisplacement.late.(epochs{i}) = netdisPlacementPCs.(epochs{i})(: ,end-3:end. :);
end

