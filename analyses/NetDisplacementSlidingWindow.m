%====================================================================================
%% net displacement 
%====================================================================================

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

%====================================================================================
%% Generating EArly and Late epoches 
%====================================================================================

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

%====================================================================================
%% generating PCs' - EArly and Late epoches
%==================================================================================== 

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

