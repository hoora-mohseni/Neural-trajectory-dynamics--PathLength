
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


