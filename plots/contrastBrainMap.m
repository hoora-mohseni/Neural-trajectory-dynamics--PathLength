%====================================================================================
%% plot hemispheres 
%====================================================================================

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

%====================================================================================
%% Plot Subcortex 
%====================================================================================

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

%====================================================================================
%% Plot Cerebellum 
%====================================================================================

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

box off;
