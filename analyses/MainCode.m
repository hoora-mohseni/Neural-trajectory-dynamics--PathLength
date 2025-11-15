
clc 
clear
close

%% Set Global Variables

global restDay1 restDay2 
global baselineDay1 baselineDay2 
global learningDay1 learningDay2 
global washoutDay1 washoutDay2 
global numSubjects epochs nPC

numSubjects = 32;
nPC         = 10;

%% Read Data

[restDay1, restDay2, baselineDay1, baselineDay2, learningDay1, ...
    learningDay2, washoutDay1, washoutDay2, epochs] = my_readData();

%% Calculate Correlation Matrices and Mean Correlation Matrix across Subjects

[meanCorrMatrices, corrMatrices] = my_calculateCorrMatrix();

%% Calculate Gradients Map

[gradientsMatrices, gradientsmap] = my_calculateGM(meanCorrMatrices);

%%  Smooth data & Calculate Projections 

pcaProjections = my_calculateProjections(gradientsMatrices);


% ========================================================================

%% Functions

function [restDay1, restDay2, baselineDay1, baselineDay2, ...
    learningDay1, learningDay2, washoutDay1, washoutDay2, epochs] = ...
    my_readData()

data_dir        = [pwd, '/../VMRD1D2_fMRIData 2/vmr_schaefer400_subcortex'];
cerebellum_dir  = [pwd, '/../VMRD1D2_fMRIData 2/nettekoven_cerebellum']; 
subject_dir     = [pwd, '/../VMRD1D2_fMRIData 2/Data4paper_Aug2020_Bins1-3 (1)'];
subjectsRange   = 'C1:C32';

files = dir(data_dir);
cerebellum_files = dir(cerebellum_dir);

% Subject list
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
    cortexData = readmatrix([data_dir, '/', matchedFiles{i}], ...
        'FileType', 'text', 'Delimiter', '\t');
    cerebellumData = readmatrix([cerebellum_dir, '/', ...
        matchedCerebellumFiles{i}], 'FileType', 'text', 'Delimiter', '\t');
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

epochs = {'restDay1', 'baselineDay1', 'learningDay1', 'washoutDay1', ...
    'restDay2', 'baselineDay2', 'learningDay2', 'washoutDay2'};

end

% ========================================================================

function [meanCorrMatrices, corrMatrices] = my_calculateCorrMatrix()

global restDay1 restDay2 baselineDay1 baselineDay2 learningDay1 ...
        learningDay2 washoutDay1 washoutDay2 epochs numSubjects

    corrMatrices = struct();
    for i = 1:numel(epochs)
        currentMatrix = eval(epochs{i}); 
        subjectCorrs = cell(numSubjects, 1);
    
        for subj = 1:numSubjects
            subjectMatrix = currentMatrix(:, :, subj);
            corrMatrix = corr(subjectMatrix);
            subjectCorrs{subj} = corrMatrix;
        end
        corrMatrices.(epochs{i}) = subjectCorrs;
    end


meanCorrMatrices = struct();

for i = 1:numel(epochs)
    subjectCorrs = corrMatrices.(epochs{i});
    [numRows, numCols] = size(subjectCorrs{1}); % Assuming all correlation matrices are square
    stackedCorrs = zeros(numRows, numCols, numSubjects);
    
    for subj = 1:numSubjects
        stackedCorrs(:, :, subj) = subjectCorrs{subj};
    end
    
    meanCorrMatrix = mean(stackedCorrs, 3); % Take mean along the 3rd dimension
    meanCorrMatrices.(epochs{i}) = meanCorrMatrix;
end



end

% ========================================================================

function [gradientsMatrices, gradientsmap] = my_calculateGM(meanCorrMatrices)
    
    global nPC epochs

    gm = GradientMaps('kernel','cs','approach','pca','alignment','pa','n_components', nPC);      %cs: cosine similarity
    
    gradientsmap = struct();
    gradientsMatrices = struct();

    for i = 1:numel(epochs)
        gradientsmap.(epochs{i})      = gm.fit(meanCorrMatrices.(epochs{i}));
        gradientsMatrices.(epochs{i}) = gradientsmap.(epochs{i}).gradients{1};
    end

end

% ========================================================================

% ========================================================================

function [pcaProjection, pcaCoeff] = my_calculateProjections(gradientsMatrices)

    global restDay1 restDay2 baselineDay1 baselineDay2 learningDay1 ...
        learningDay2 washoutDay1 washoutDay2 epochs numSubjects


    smoothedData    = struct();
    polynomialOrder = 1;
    frameLength     = 3;  
    
    % Apply Savitzky-Golay filter and PCA projection for each phase
    for i = 1:numel(epochs)
        phase = epochs{i};
        for subj = 1:numSubjects
            currentMatrix                   = eval(phase); 
            smoothedData.(phase)(:,:,subj)  = sgolayfilt(currentMatrix(:,:,subj), ...
                polynomialOrder, frameLength, [], 1);
            
            % Project the smoothed data into PCA space
            for TR = 1:size(smoothedData.(phase), 1)
                pcaProjections.(phase)(TR, :, subj) = ...
                    corr(squeeze(smoothedData.(phase)(TR,:,subj))', gradientsMatrices.baselineDay1);
            end

            [coeff, score, latent, ~, explained]        = pca(smoothedData.(phase)(:,:,subj));
            latentAllSubj(1:length(latent), subj)       = latent;
            explainedAllSubj(1:length(explained), subj) = explained;
            pcaProjection.(phase)(:,:,subj)             = score(:, 1:10);  % Using the first 10 principal components
            pcaCoeff.(phase)(:, :, subj)                = coeff(:, 1:10);
            
        end
    end

end


% ========================================================================




