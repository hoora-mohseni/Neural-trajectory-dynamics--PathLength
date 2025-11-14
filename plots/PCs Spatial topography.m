%%  plot hemispheres 

ParcelNumber = 400;
[surf_lh, surf_rh] = load_conte69();

% plot hemispheres of Conte69
% plot_hemispheres(ones(64984,1),{surf_lh,surf_rh});

% Load connectivity matrix from subset of HCP data; you can specify parcellation from [100, 200, 300, 400]
labeling = load_parcellation('schaefer',ParcelNumber);
conn_matrices = load_group_fc('schaefer',ParcelNumber);

% Plot the gradients

for i = 1:numel(phases)
    plot_hemispheres(gradientsMatrices.(phases{i})(33:end-32,1:5), {surf_lh,surf_rh}, ...
           'parcellation',labeling.schaefer_400);
    title([phases{i}, ' - Gradient 1 to 5']);
    
%     plot_hemispheres(gradientsMatrices.(phases{i})(33:end-32,6:10), {surf_lh,surf_rh}, ...
%            'parcellation',labeling.schaefer_400);
%     title([phases{i}, ' - Gradient 6 to 10']);
%    
end


%% Plot Subcortex 

subcortexData = struct();
numPCs          = 10;

for i = 1:numel(phases)
    subcortex = gradientsMatrices.(phases{i})(1:32,:);
    new_order = [17:32, 1:16];
    subcortexData.(phases{i}) = subcortex( new_order , :);
    Mymin(i) = min(min(subcortexData.(phases{i})));
    Mymax(i) = max(max(subcortexData.(phases{i})));
end

Mymin = min(Mymin);
Mymax = max(Mymax);

% Load the precomputed subcortical surface data
load('subcorticalSurf.mat', 'VL', 'VR', 'nVL', 'nVR');

% Visualize the subcortex for each phase using plotSubcortex.m

for pc = 1:numPCs 
    for i = 2:numel(phases) 
        figure; 
        hold on;
        nexttile;    
        colormap(parula);
        %colorbar;
        axis off;
        currentSubcortex = subcortexData.(phases{i})(:,pc); 
        plotSubcortex(currentSubcortex, [Mymin,Mymax] );  
        title(['Subcortex - PC',pc ,phases{i}]);
    end
end

%% Plot Cerebellum

cerebellumData = struct();
numPCs          = 10;

for i = 1:numel(phases)
	cerebellumData.(phases{i}) = gradientsMatrices.(phases{i})(433:464,:);  
    Mymin(i) = min(min(cerebellumData.(phases{i})));
    Mymax(i) = max(max(cerebellumData.(phases{i})));
end

Mymin = min(Mymin);
Mymax = max(Mymax);


for pc = 1:5%:numPCs  
    for i = 7%:numel(phases) 
        figure; 
        hold on;
        nexttile;    
        colormap(parula);
        %colorbar;
        axis off;
        currentCerebellum = cerebellumData.(phases{i})(:,pc); 
        plotCerebellum(currentCerebellum, [Mymin, Mymax] );  
        title(['Cerebellum - PC', pc ,phases{i}]);
    end
end
