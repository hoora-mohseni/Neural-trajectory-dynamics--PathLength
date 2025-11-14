%====================================================================================
%% --- Plot Motion without Overlap ---
%====================================================================================

numSubjects = numel(subjects);
figure('color', 'w', 'Position', [300, 300, 1400, 400]);

x_s     = 5; 
x_total = x_s;
red_dark   = [1 0 0];

for i = 1:4
  
    s1 = shadedErrorBar(x_s:x_s+epochWindowCounts(i)-1, mean(dvarsSliding.(epochs{i})), ...
                           std(dvarsSliding.(epochs{i}))/sqrt(numSubjects), ...
                            {'-', 'color', [0.3, 0.05, 0.4], 'LineWidth', 2}, 1, 1);
    
    hold on
    
    s2 = shadedErrorBar(x_s:x_s+epochWindowCounts(i)-1, mean(dvarsSliding.(epochs{i+4})), ...
                           std(dvarsSliding.(epochs{i+4}))/sqrt(numSubjects), ...
                           {'-', 'color',  [0.1 ,  0.6   , 0.4], 'LineWidth', 1.5}, 1, 1);

    x_s = x_s + epochWindowCounts(i) + 10; 
    
    x_total = [x_total, x_s-10];
    x_total = [x_total, x_s];
   
end

box off
set(gca, "FontSize", 20, 'LineWidth', .8, 'TickDir', 'out')
yticks(20:2:30); ylim([20, 30])
xticks(x_total); xticklabels([1, 21, 22, 51, 52, 131, 132, 161])
xlabel("Bins"); pl = ylabel("dvars"); 
xlim([0, x_total(end-1)])
pl.Position = pl.Position - [3, 0, 0];
legend([s1.patch, s2.patch], {"Day 1", "Day 2"}, "Location", "NorthEast")
legend box off

%====================================================================================
%% pathLengh and head motion correlation- Day 1
%====================================================================================

figure('color', 'w', 'Position', [300, 300, 1000, 400]);

myColor = viridis(137);
allclr = [];
colorIdx = 1;
r_path = [];
r_Dvars = [];

for i = 2:4
    allPath = mean(pathLengthSliding.(epochs{i}));
    Dvars     = mean(hmSliding.(epochs{i}));


    if i == 2 
        sign = 'o' ;
    elseif i == 3 
        sign = "square";
    else
        sign = '^';
    end

    for j = 1:numel(Dvars)
       clr = myColor(colorIdx, :);
        allclr = [allclr; clr];
        a(i) = plot(Dvars(j), allPath(j), sign, 'MarkerFaceColor', clr, 'MarkerEdgeColor', clr); hold on
        colorIdx = colorIdx + 1;
    end
    r_Dvars =  [r_Dvars, Dvars];
    r_path = [r_path, allPath];
end

% Remove NaNs before correlation
validIdx = ~isnan(r_Dvars) & ~isnan(r_path);
[R, P] = corrcoef(r_Dvars(validIdx), r_path(validIdx));
R_value = R(1,2);
P_value = P(1,2);

% Display R and p in command window
fprintf('R value: %.3f\n', R_value);
fprintf('p value: %.3f\n', P_value);
text(0.05, 0.9, sprintf('R = %.2f\np = %.3f', R_value, P_value), ...
    'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold');

mdl = fitlm(r_Dvars, r_path);
plot(r_Dvars, mdl.Fitted, 'color', myColor(20, :), 'LineWidth', 2);

% yticks(70:10:120)
% xticks(0:15:45)
% xlim([0, 45]); %ylim([70, 120])
box off; 
colormap(allclr);
c = colorbar;
c.Ticks = [];        % Removes both ticks and labels
set(gca, "FontSize", 20, 'LineWidth', .8, 'TickDir', 'out')
axis square
%legend([a(2), a(3), a(4)], {' ', ' ', ' '});

%====================================================================================
%%
%====================================================================================

figure('color', 'w')
hmAll = [];
pthAll = [];
for i = 1:numel(epochs)
    hmAll = [hmAll, hmSliding.(epochs{i})];
    pthAll = [pthAll, pathLengthSliding.(epochs{i})];
end


for i = 1:5
    if i == 5
        a = hmAll;
        b = pthAll;
    else
        a = hmSliding.(epochs{i});
        b = pathLengthSliding.(epochs{i});
    end
    
    mycorr(1:32, i) = diag(corr(a', b'));

end

b = bar(mean(mycorr)); hold on

for i = 1:5
    plot(rand(32, 1)*.2+.78+i-1, mycorr(:, i), 'k.', 'MarkerSize', 20)
    errorbar(i+.2, mean(mycorr(:, i)), std(mycorr(:, i))/sqrt(32), 'k', 'LineWidth', 2)
end

xticks(1:5); xticklabels([epochs(1:4), 'Total'])
ylabel('Correlation Path Length - Head Motion')
xlabel('epoch')
box off

% --- mapping ---
epochs  = {'restDay1','baselineDay1','learningDay1','washoutDay1', ...
           'restDay2','baselineDay2','learningDay2','washoutDay2'};
phaseOf = containers.Map(epochs, {'rest','baseline','learning','washout', ...
                                  'rest','baseline','learning','washout'});
dayOf   = containers.Map(epochs, {'day1','day1','day1','day1', ...
                                  'day2','day2','day2','day2'});

% --- gather rows: one row per (Subject, Phase, Day, Window) ---
FD = [];  Subject = [];  Phase = {};  Day = {};  Window = [];

numSubjects = numel(subjects);
for e = 1:numel(epochs)
    M = hmSliding.(epochs{e});          % [numSubjects x numWindows_e]
    [nS, nW] = size(M);

    FD      = [FD;      M(:)];
    Subject = [Subject; repelem((1:nS)', nW, 1)];
    Window  = [Window;  repmat((1:nW)', nS, 1)];
    Phase   = [Phase;   repmat({phaseOf(epochs{e})}, nS*nW, 1)];
    Day     = [Day;     repmat({dayOf(epochs{e})},   nS*nW, 1)];
end

T = table;
T.FD      = FD;
T.Subject = categorical(Subject);
T.Phase   = categorical(Phase, {'rest','baseline','learning','washout'}); % set reference order
T.Day     = categorical(Day,   {'day1','day2'});
T.Window  = Window;   % (not used in the fixed part of the model)

%====================================================================================
%% Build a table of window-level FD
%====================================================================================

% --- mapping ---
epochs  = {'restDay1','baselineDay1','learningDay1','washoutDay1', ...
           'restDay2','baselineDay2','learningDay2','washoutDay2'};
phaseOf = containers.Map(epochs, {'rest','baseline','learning','washout', ...
                                  'rest','baseline','learning','washout'});
dayOf   = containers.Map(epochs, {'day1','day1','day1','day1', ...
                                  'day2','day2','day2','day2'});

% --- gather rows: one row per (Subject, Phase, Day, Window) ---
FD = [];  Subject = [];  Phase = {};  Day = {};  Window = [];

numSubjects = numel(subjects);
for e = 1:numel(epochs)
    %M = hmSliding.(epochs{e});          % [numSubjects x numWindows_e]
    M = dvarsSliding.(epochs{e}); 
    [nS, nW] = size(M);

    FD      = [FD;      M(:)];
    Subject = [Subject; repelem((1:nS)', nW, 1)];
    Window  = [Window;  repmat((1:nW)', nS, 1)];
    Phase   = [Phase;   repmat({phaseOf(epochs{e})}, nS*nW, 1)];
    Day     = [Day;     repmat({dayOf(epochs{e})},   nS*nW, 1)];
end

T = table;
T.FD      = FD;
T.Subject = categorical(Subject);
T.Phase   = categorical(Phase, {'rest','baseline','learning','washout'}); % set reference order
T.Day     = categorical(Day,   {'day1','day2'});
T.Window  = Window;   % (not used in the fixed part of the model)

%====================================================================================
%% plot multiplt regrassion (pathlength, error, dvars)
%====================================================================================

validIdx = ~isnan(r_err) & ~isnan(r_path) & ~isnan(r_Dvars);
Error = zscore(r_err(validIdx))';
Dvars = zscore(r_Dvars(validIdx)');
PathLength = zscore(r_path(validIdx)');

% Build table
tbl = table(Error, Dvars, PathLength,'VariableNames', {'Error', 'Dvars', 'PathLength'});

% Fit model
mdl = fitlm(tbl, 'PathLength ~ Error + Dvars');

disp(mdl);
% Create grid for surface
[eGrid, rtGrid] = meshgrid(linspace(min(Error), max(Error), 30), ...
                           linspace(min(Dvars), max(Dvars), 30));

% Predict over the grid
predictedPath = mdl.Coefficients.Estimate(1) + ...
                mdl.Coefficients.Estimate(2) * eGrid + ...
                mdl.Coefficients.Estimate(3) * rtGrid;

% 3D plot
figure('Color', 'w');
hold on
surf(eGrid, rtGrid, predictedPath, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

% Scatter actual data
scatter3(Error, Dvars, PathLength, 60, PathLength, 'filled');

xlabel('Error');
ylabel('FD');
zlabel('Path Length');
title('Multiple Linear Regression: PathLength ~ Error + FD');
colorbar;
grid on;
view(135, 25);  % Adjust the 3D view angle
set(gca, 'FontSize', 14);
axis tight



