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

