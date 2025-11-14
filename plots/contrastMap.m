clc 
clear
close

%% Brain Map Generation


earlyRestD1 = restDay1(1:16, :, :);
earlyRestD2 = restDay2(1:16, :, :);
lateRestD1 = restDay1(end-16+1:end, :, :);
lateRestD2 = restDay2(end-16+1:end, :, :);

earlyBaselineD1 = baselineDay1(1:16, :, :);
earlyBaselineD2 = baselineDay2(1:16, :, :);
lateBaselineD1 = baselineDay1(end-16+1:end, :, :);
lateBaselineD2 = baselineDay2(end-16+1:end, :, :);

earlyLearningD1 = learningDay1(1:16, :, :);
earlyLearningD2 = learningDay2(1:16, :, :);
lateLearningD1 = learningDay1(end-16+1:end, :, :);
lateLearningD2 = learningDay2(end-16+1:end, :, :);

earlyWashoutD1 = washoutDay1(1:16, :, :);
earlyWashoutD2 = washoutDay2(1:16, :, :);
lateWashoutD1 = washoutDay1(end-16+1:end, :, :);
lateWashoutD2 = washoutDay2(end-16+1:end, :, :);


%% mean activation

contrastD1 = squeeze(mean(lateLearningD1, 1) - mean(earlyLearningD1, 1));  
contrastD2 = squeeze(mean(lateLearningD2, 1) - mean(earlyLearningD2, 1));
contrastData = {contrastD1, contrastD2};
numRegions = size(contrastD1, 1);
numContrasts = numel(contrastData);

p_values = zeros(numRegions, numContrasts);
t_values = zeros(numRegions, numContrasts);

for c = 1:numContrasts
    contrast = contrastData{c};
    for i = 1:numRegions
        [~, p, ~, stats] = ttest(contrast(i, :));
        p_values(i, c) = p;
        t_values(i, c) = stats.tstat;
    end
end

% FDR correction
p_fdr = zeros(size(p_values));
for c = 1:numContrasts
    p_fdr(:, c) = mafdr(p_values(:, c), 'BHFDR', true);
end
%t_values(p_fdr > 0.05) = 0;


%% plot hemispheres 

%max_abs_t = max(abs(t_values), [], 'all'); 
clim = [-6.6, 6.6];


ParcelNumber = 400;
[surf_lh, surf_rh] = load_conte69();

% Load connectivity matrix from subset of HCP data; you can specify parcellation from [100, 200, 300, 400]
labeling = load_parcellation('schaefer',ParcelNumber);
conn_matrices = load_group_fc('schaefer',ParcelNumber);

%t_values(t_values == 0) = NaN;
writematrix(t_values, 't_values.csv')

t_map = t_values(33:end-32, :);

t_map(t_map(:, 1) == min(t_map(:, 1)), 1) = -max(t_map(:, 1));

obj = plot_hemispheres(t_map, {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400, 'views', 'lm', ...
    'labeltext', []);

cmap = RedBlueCmap();
cmap(1, :) = .9*[1, 1, 1];

% Apply colormap only to cortical areas
% colormaps(obj, cmap);  % Grey out NaNs
% set([obj.axis], 'Clim', clim);


obj.colormaps({cmap, cmap});
obj.colorlimits([clim; clim]);


%% Plot Subcortex 


subcortex = t_values(1:32,:);
new_order = [17:32, 1:16];
subcortexData = subcortex( new_order , :);
% Mymin = min(min(subcortexData));
% Mymax = max(max(subcortexData));

% Load the precomputed subcortical surface data
load('subcorticalSurf.mat', 'VL', 'VR', 'nVL', 'nVR');
cmap = RedBlueCmap();

figure('color', 'w') 
hold on;
nexttile;    
colormap(cmap);
%colorbar;
axis off;
currentSubcortex = subcortexData(:,2);
plotSubcortex(currentSubcortex, [-6.6,6.6] );  

 %% Plot Cerebellum 

cerebellumData = t_values(433:464,:);  
Mymin = min(min(cerebellumData));
Mymax = max(max(cerebellumData));

for i= 1:2
    figure('color', 'w'); 
    hold on;
    nexttile;    
    colormap(cmap);
    %colorbar;
    axis off;
    plotCerebellum(cerebellumData(:,i), [-6.6,6.6] );  

end


%%

lbl = readtable([pwd, '/schaefer_2018/schaefer400NodeNames.txt']);

sig_areas_d1_inc = lbl(t_map(:, 1) > 0, :);
sig_areas_d2_inc = lbl(t_map(:, 2) > 0, :);

sig_areas_d1_dec = lbl(t_map(:, 1) < 0, :);
sig_areas_d2_dec = lbl(t_map(:, 2) < 0, :);

%% 

W = 16; % Window size
S = 8;  % Step size


for subj = 1:numSubjects
    for i = 1:numel(epochs)
        phase = epochs{i};
        currentMatrix = eval(phase); 
        epochData = currentMatrix(:,:,subj);
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
            BOLD.(epochs{i})(subj, win, : ) = mean(epochData(startIdx:endIdx, :));
        end
    end
end

%% BOLD and PL Correlation

BOLDcatD1 = BOLD.(epochs{1});
PLcatD1   = pathLengthSliding.(epochs{1});
BOLDcatD2 = BOLD.(epochs{5});
PLcatD2   = pathLengthSliding.(epochs{5});
for i = 2:4
    BOLDcatD1 = cat(2, BOLDcatD1, BOLD.(epochs{i}));
    PLcatD1   = cat(2, PLcatD1,   pathLengthSliding.(epochs{i}));
    BOLDcatD2 = cat(2, BOLDcatD2, BOLD.(epochs{i+4}));
    PLcatD2   = cat(2, PLcatD2,   pathLengthSliding.(epochs{i+4}));
end

for subj = 1:numSubjects
    for ROI = 1:464
        
        r1(subj, ROI) = corr(BOLDcatD1(subj, :, ROI)', PLcatD1(subj, :)');
        r2(subj, ROI) = corr(BOLDcatD2(subj, :, ROI)', PLcatD2(subj, :)');
    end
end

Rmean = [mean(r1,1).',  mean(r2,1).']; 

%% plot hemispheres 


ParcelNumber = 400;
[surf_lh, surf_rh] = load_conte69();

% Load connectivity matrix from subset of HCP data; you can specify parcellation from [100, 200, 300, 400]
labeling = load_parcellation('schaefer',ParcelNumber);
conn_matrices = load_group_fc('schaefer',ParcelNumber);

R_map = Rmean(33:end-32, :);

R_map(R_map(:, 1) == min(R_map(:, 1)), 1) = -max(R_map(:, 1));

obj = plot_hemispheres(R_map, {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_400,'views', 'lm', ...
    'labeltext', []);

cmap = RedBlueCmap();
cmap(1, :) = .9*[1, 1, 1];

% Apply colormap only to cortical areas
colormaps(obj, cmap);  % Grey out NaNs

%% Plot Subcortex 


subcortex = Rmean(1:32, :);
new_order = [17:32, 1:16];
subcortexData = subcortex( new_order , :);
Mymin = min(min(subcortexData));
Mymax = max(max(subcortexData));

% Load the precomputed subcortical surface data
load('subcorticalSurf.mat', 'VL', 'VR', 'nVL', 'nVR');
cmap = RedBlueCmap();

figure('color', 'w') 
hold on;
nexttile;    
colormap(cmap);
%colorbar;
axis off;
currentSubcortex = subcortexData(:,1);
plotSubcortex(currentSubcortex, [-4,4] );  

 %% Plot Cerebellum 

cerebellumData = Rmean(433:464,:);  
Mymin = min(min(cerebellumData));
Mymax = max(max(cerebellumData));

for i= 1:2
    figure('color', 'w'); 
    hold on;
    nexttile;    
    colormap(cmap);
    %colorbar;
    axis off;
    plotCerebellum(cerebellumData(:,i), [-4,4] );  

end

%% Correlation Matrix
figure('color', 'w')


corrmatrix_tvalues_6Con = [t_values_EB_LR, tValue_EL_LB, tValue_LL_EL];
corrmatrix_tvalues_4Con = [tValue_EL_LB, tValue_LL_EL];

corrmatrix = corr(corrmatrix_tvalues_4Con);

imagesc(corrmatrix); title('Correlation Matrix'); colormap(viridis); caxis([-1, 1]); colorbar; axis square;
xlabel('contrasts map t values'); ylabel('contrasts map t values');
ax = gca;
ax.XTick = [];
ax.YTick = [];

[rows, cols] = size(corrmatrix);

for i = 1:rows
    for j = 1:cols
        R_value_str = num2str(corrmatrix(i, j), '%.2f');
        current_value = corrmatrix(i, j);
        
        if abs(current_value) < 0.25 
            text_color = [0 0 0]; % Black text for mid-range (often yellow/green)
        else
            text_color = [1 1 1]; % White text for high/low (dark purple/bright yellow)
        end
        
        text(j, i, R_value_str, ...
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle', ...
             'Color', text_color, ...
             'FontSize', 10); % Adjust font size as needed
    end
end

%% UMAP 


%Convert the Similarity (Correlation) Matrix to a Dissimilarity (Distance) Matrix
distance_matrix = 1 - corrmatrix;

% Use cmdscale (Classical MDS) for dimensionality reduction
% Input: distance_matrix
% Output: Y is the 4x2 matrix of coordinates, eigenvalue is the goodness-of-fit
% 2 is the desired number of dimensions (for a 2D plot)
[Y, eigenvalue] = cmdscale(distance_matrix, 2);

% Y is your 4x2 coordinate matrix. Plotting the rows of Y will show your map.

figure('Color', 'w', 'Position', [300, 300, 400, 400]);
hold on; 
plot(Y(1, 1), Y(1, 2),  'hexagram', 'MarkerSize', 15, 'MarkerFaceColor', [0.3, 0.05, 0.4], 'MarkerEdgeColor', [0.3, 0.05, 0.4]);   
plot(Y(2, 1), Y(2, 2),   'hexagram', 'MarkerSize', 15, 'MarkerFaceColor', [0.1 , 0.6  , 0.4], 'MarkerEdgeColor', [0.1 , 0.6  , 0.4]);   
plot(Y(3, 1), Y(3, 2),  's', 'MarkerSize', 15, 'MarkerFaceColor', [0.3, 0.05, 0.4], 'MarkerEdgeColor', [0.3, 0.05, 0.4]);   
plot(Y(4, 1), Y(4, 2),   's', 'MarkerSize', 15, 'MarkerFaceColor', [0.1 , 0.6  , 0.4], 'MarkerEdgeColor', [0.1 , 0.6  , 0.4]);
plot(Y(5, 1), Y(5, 2),  'o', 'MarkerSize', 15, 'MarkerFaceColor', [0.3, 0.05, 0.4], 'MarkerEdgeColor', [0.3, 0.05, 0.4]);   
plot(Y(6, 1), Y(6, 2),   'o', 'MarkerSize', 15, 'MarkerFaceColor', [0.1 , 0.6  , 0.4], 'MarkerEdgeColor', [0.1 , 0.6  , 0.4]);

text(Y(:,1), Y(:,2), {'EB-LR D1', 'EB-LR D2',  'EL-LB D1', 'EL-LB D2', 'LL-EL D1', 'LL-EL D2'}, ...
     'VerticalAlignment','top', 'HorizontalAlignment','center'); %6 Condition
% text(Y(3:6,1), Y(3:6,2), {'EL-LB D1', 'EL-LB D2', 'LL-EL D1', 'LL-EL D2'}, ...
%      'VerticalAlignment','top', 'HorizontalAlignment','center'); %4 condition
title('MDS');
xlabel('Dimension 1'); xticks(-1:0.5:1); xlim([-1,1]); 
ylabel('Dimension 2'); yticks(-1:0.5:1); ylim([-1,1]); 
grid off;
axis equal;
box off;

%% Averagted t-values within each network barplot

% Extracting Major Network indices
Ntwrk = {'vis', 'sommot', 'dorsattn', 'salventattn', 'limbic', 'cont', 'default'}; 
for i = 1:numel(Ntwrk)
    NtwrkIdx{i} = find(cellfun(@(x) count(lower(x), lower(Ntwrk{i})), headers) == 1);
end

NtwrkContrast = struct();
for n = 1:numel(Ntwrk)
    targetLoadings = corrmatrix_tvalues_6Con(NtwrkIdx{n}, :);
    NtwrkContrast.(Ntwrk{n}) = mean(targetLoadings, 1);
end

%%
% Define the network names (for X-axis labels)
Ntwrk = {'vis', 'sommot', 'dorsattn', 'salventattn', 'limbic', 'cont', 'default'}; 
NumNetworks = numel(Ntwrk);

% --- 1. Consolidate Data into a Matrix (7 Networks x 6 Conditions) ---
% Rows will be Networks, Columns will be the 6 Conditions
PlotData = zeros(NumNetworks, 3);

for n = 1:NumNetworks
    % Extract the 1x6 vector for the current network and place it in the matrix
    % NOTE: This assumes NtwrkContrast has been correctly built in the previous step
    PlotData(n, :) = NtwrkContrast.(Ntwrk{n})(:,[1,3,5]); % day 1
    %PlotData(n, :) = NtwrkContrast.(Ntwrk{n})(:,[2,4,6]);  % day 2
end

% --- 2. Create the Grouped Bar Plot ---
figure('color', 'w', 'Position', [300, 300, 1400, 800]);
% The bar function plots each row (network) along the x-axis, with 
% grouped bars for each column (condition).
h = bar(PlotData); 
hold on;

% --- 3. Apply Specific Colors Manually ---
% h is an array of 6 graphic objects (one for each condition/group).
% We set the color based on the condition index:


h(1).FaceColor = [0.3, 0.05, 0.4]; % Day 1 
h(2).FaceColor = [0.3, 0.05, 0.4];
h(3).FaceColor = [0.3, 0.05, 0.4];

% Day 2 
% h(1).FaceColor = [0.1 , 0.6  , 0.4]; % Day 2 
% h(2).FaceColor = [0.1 , 0.6  , 0.4];
% h(3).FaceColor = [0.1 , 0.6  , 0.4];

% --- 4. Add Labels and Formatting ---
% Set the X-axis labels to the network names
set(gca, 'XTick', 1:NumNetworks, 'XTickLabel', Ntwrk, 'FontSize', 15);
title('Averaged T-values per Network and Condition');
ylabel('Mean T-value'); %ylim([-2.5, 2.5])
xlabel('Network');

% Assume conditions are paired across A, B, C phases of your learning model
LegendLabels = {'E-Baseline vs. L-Rest Day 2', 'E-Learning vs. L-Baseline Day 2',  ...
                'L-Learning vs. E-Learning Day 2'};

% Add a legend using the bar objects (h) and custom labels
legend(h, LegendLabels, 'Location', 'northeast');
grid off;
box off;
