
clc
clear
close all

%%

hm = [];
stage = {'rest_run-2', 'rotation',  'washout'};
phases = {'rest', 'rotation',  'washout'};
session = {'day1', 'day2'};

%% Subject list
subject_dir = [pwd, '/../VMRD1D2_fMRIData 2/Data4paper_Aug2020_Bins1-3 (1)'];
subjectsRange = 'C1:C32';
subjects = readcell(subject_dir, 'Range', subjectsRange);

for subj = 1:numel(subjects)
    for ses = 1:numel(session)
        for i = 1:length(stage)
            pth = ['/Volumes/Raw/VMR-Learning-Complete/derivatives/2020/fmriprep/',subjects{subj},'/ses-0', num2str(ses),'/func/', subjects{subj},'_ses-0', ...
                num2str(ses),'_task-', stage{i},'_desc-confounds_regressors.tsv'];
            
            T = readtable(pth, 'FileType', 'text', 'Delimiter', '\t');
            T = T(:, {'framewise_displacement', 'dvars'});
            hm.(session{ses}).(phases{i})(1:size(T), subj) = T.framewise_displacement;
            dvars.(session{ses}).(phases{i})(1:size(T), subj) = T.dvars;
        end
    end 
end

%% save
save('hm_data.mat', 'hm', '-v7.3');  % -v7.3 if it’s big
save('dvars_data.mat', 'dvars', '-v7.3');  % -v7.3 if it’s big

%% load
load('hm_data.mat', 'hm');
load('dvars_data.mat', 'dvars');

%% 

epochs = {'restDay1', 'baselineDay1', 'learningDay1', 'washoutDay1', 'restDay2', 'baselineDay2', 'learningDay2', 'washoutDay2'};
headM = [];
headM2 = [];

headM.(epochs{1}) = hm.(session{1}).(phases{1})(:, :);
headM.(epochs{2}) = hm.(session{1}).(phases{2})(9:248, :);
headM.(epochs{3}) = hm.(session{1}).(phases{2})(248+1:248+640, :);
headM.(epochs{4}) = hm.(session{1}).(phases{3})(9:248, :);

headM.(epochs{5}) = hm.(session{2}).(phases{1})(:, :);
headM.(epochs{6}) = hm.(session{2}).(phases{2})(9:248, :);
headM.(epochs{7}) = hm.(session{2}).(phases{2})(248+1:248+640, :);
headM.(epochs{8}) = hm.(session{2}).(phases{3})(9:248, :);

headM2.(epochs{1}) = dvars.(session{1}).(phases{1})(:, :);
headM2.(epochs{2}) = dvars.(session{1}).(phases{2})(9:248, :);
headM2.(epochs{3}) = dvars.(session{1}).(phases{2})(248+1:248+640, :);
headM2.(epochs{4}) = dvars.(session{1}).(phases{3})(9:248, :);

headM2.(epochs{5}) = dvars.(session{2}).(phases{1})(:, :);
headM2.(epochs{6}) = dvars.(session{2}).(phases{2})(9:248, :);
headM2.(epochs{7}) = dvars.(session{2}).(phases{2})(248+1:248+640, :);
headM2.(epochs{8}) = dvars.(session{2}).(phases{3})(9:248, :);


%%

W = 16; % Window size
S = 8;  % Step size

hmSliding = [];
epochWindowCounts = zeros(1, numel(epochs)); 

for subj = 1:numel(subjects)
    for i = 1:numel(epochs)
        epochData = headM.(epochs{i})(:, subj);
        pochData2 = headM2.(epochs{i})(:, subj);
        epochLength = size(epochData, 1);  % Determine actual time length of epoch

        numWindows = max(0, floor((epochLength - W) / S) + 1);
        epochWindowCounts(i) = numWindows; 
        
        for win = 1:numWindows
            startIdx = (win - 1) * S + 1;
            endIdx   = startIdx + W - 1;
            
            if endIdx > epochLength
                continue; 
            end
            
            hmSliding.(epochs{i})(subj, win) = mean(epochData(startIdx:endIdx, :), 'omitnan');
            dvarsSliding.(epochs{i})(subj, win) = mean(pochData2 (startIdx:endIdx, :), 'omitnan');

        end
    end
end


%% --- Plot Motion without Overlap ---

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

%% pathLengh and head motion correlation- Day 1

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


%%

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

%% Build a table of window-level FD

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

%% Fit M0: FD ~ Phase * Day + (1 | Subject)

lme_M0 = fitlme(T, 'FD ~ Phase*Day + (1|Subject)');  % random intercept per subject

% Inference on fixed effects (Phase, Day, Phase×Day):
anova(lme_M0, 'DFMethod','Satterthwaite')

%% Get adjusted means (Phase × Day) for plotting

% Reference grid for fixed-effects predictions (random effects off)
ph = categories(T.Phase); dy = categories(T.Day);
[PH, DY] = ndgrid(ph, dy);
G = table;
G.Phase   = categorical(PH(:), ph);
G.Day     = categorical(DY(:), dy);
G.Subject = T.Subject(1);              % dummy; ignored when Conditional=false

[pred, predCI] = predict(lme_M0, G, 'Conditional', false);  % fixed-effects means
means = reshape(pred, [numel(ph) numel(dy)]);
loCI  = reshape(predCI(:,1), [numel(ph) numel(dy)]);
hiCI  = reshape(predCI(:,2), [numel(ph) numel(dy)]);

% quick bar plot (Phase on x, bars = Day)
figure('Color','w'); hold on
b = bar(means);  % 4 phases × 2 days
ng = size(means,1); nb = size(means,2);
for ii = 1:ng
    for jj = 1:nb
        x = b(jj).XEndPoints(ii);
        plot([x x], [loCI(ii,jj) hiCI(ii,jj)], 'k-', 'LineWidth', 1.5);
    end
end
set(gca,'XTickLabel', ph); ylabel('Predicted FD (mm)');
legend(dy,'Location','northwest'); box off
title('FD by Phase × Day (fixed-effects estimates)')


%% ---------- Build table: one row per (Subject, Phase, Day, Window) ----------

epochs  = {'restDay1','baselineDay1','learningDay1','washoutDay1', ...
           'restDay2','baselineDay2','learningDay2','washoutDay2'};
phaseOf = containers.Map(epochs, {'rest','baseline','learning','washout', ...
                                  'rest','baseline','learning','washout'});
dayOf   = containers.Map(epochs, {'day1','day1','day1','day1', ...
                                  'day2','day2','day2','day2'});

PL = []; FD = []; DV = [];
Subj = []; Phase = {}; Day = {}; Win = [];

numSubjects = numel(subjects);
for e = 1:numel(epochs)
    PL_e = pathLengthSliding.(epochs{e});         % [S × W_e]
    FD_e = hmSliding.(epochs{e});                 % [S × W_e]
    DV_e = dvarsSliding.(epochs{e});              % [S × W_e]
    [S, W] = size(PL_e);

    PL   = [PL;   PL_e(:)];
    FD   = [FD;   FD_e(:)];
    DV   = [DV;   DV_e(:)];
    Subj = [Subj; repelem((1:S)', W, 1)];
    Win  = [Win;  repmat((1:W)', S, 1)];
    Phase = [Phase; repmat({phaseOf(epochs{e})}, S*W, 1)];
    Day   = [Day;   repmat({dayOf(epochs{e})},   S*W, 1)];
end

T = table;
T.PathLen = PL;
T.FD      = FD;
T.DVARS   = DV;
T.Subj    = categorical(Subj);
T.Phase   = categorical(Phase, {'rest','baseline','learning','washout'}); % set order
T.Day     = categorical(Day,   {'day1','day2'});
T.Window  = Win;

% Optional: remove any rows with missing values (should be none)
T = rmmissing(T);

% (Optional but recommended) Standardize continuous predictors within-subject for interpretability
Tz = T;
for s = categories(T.Subj)'
    idx = (T.Subj==s);
    Tz.FD(idx)    = (T.FD(idx) - mean(T.FD(idx),'omitnan')) ./ std(T.FD(idx),[],'omitnan');
    Tz.DVARS(idx) = (T.DVARS(idx) - mean(T.DVARS(idx),'omitnan')) ./ std(T.DVARS(idx),[],'omitnan');
end

%% M2 Core control

% FD-only adjustment
m2_fd = fitlme(T, 'PathLen ~ Phase*Day + FD + (1|Subj)');
disp(anova(m2_fd, 'DFMethod','Satterthwaite'));      % tests Phase/Day/Interaction and FD jointly
[beta,~,FE] = fixedEffects(m2_fd,'DFMethod','Satterthwaite');
disp(table(m2_fd.CoefficientNames', beta, FE.SE, FE.DF, FE.pValue, ...
    'VariableNames',{'Term','Beta','SE','DF','p'}));

%% FD + DVARS 

% Add DVARS to avoid omitted-variable bias
m2_both = fitlme(T, 'PathLen ~ Phase*Day + FD + DVARS + (1|Subj)');
disp(anova(m2_both, 'DFMethod','Satterthwaite'));
[beta,~,FE] = fixedEffects(m2_both,'DFMethod','Satterthwaite');
disp(table(m2_both.CoefficientNames', beta, FE.SE, FE.DF, FE.pValue, ...
    'VariableNames',{'Term','Beta','SE','DF','p'}));

%%  Plot adjusted Phase×Day means (U-shape after control)

% Reference grid for predictions (random effects off)
ph = categories(T.Phase); dy = categories(T.Day);
[PH, DY] = ndgrid(ph, dy);
G = table;
G.Phase = categorical(PH(:), ph);
G.Day   = categorical(DY(:), dy);
G.Subj = repmat(T.Subj(1), height(G), 1);   % repeat the first subject for every row
  % dummy; ignored when Conditional=false
G.FD    = repmat(mean(T.FD,'omitnan'),    height(G), 1);
G.DVARS = repmat(mean(T.DVARS,'omitnan'), height(G), 1);

[pred, ci] = predict(m2_both, G, 'Conditional', false);  % fixed-effects means
means = reshape(pred, [numel(ph) numel(dy)]);
lo    = reshape(ci(:,1), [numel(ph) numel(dy)]);
hi    = reshape(ci(:,2), [numel(ph) numel(dy)]);

figure('Color','w'); b = bar(means); hold on
for i=1:size(means,1), for j=1:size(means,2)
        x = b(j).XEndPoints(i); plot([x x],[lo(i,j) hi(i,j)],'k-','LineWidth',1.5);
    end, end
set(gca,'XTickLabel',ph); ylabel('Predicted Path Length'); legend(dy,'Location','northwest'); box off
title('Path Length by Phase × Day (adjusted for FD + DVARS)')


%% ---------- Compute within/between components ----------
% Step 3 — M3 Within- vs between-subject motion
% I split motion into subject means (between) and deviations (within). Did this for FD and (optionally) DVARS.
T_M3 = T;  % or Tz if you preferred z-scales
% FD
FD_between   = zeros(height(T_M3),1);
FD_within    = zeros(height(T_M3),1);
% DVARS
DV_between   = zeros(height(T_M3),1);
DV_within    = zeros(height(T_M3),1);

for s = categories(T_M3.Subj)'
    idx = (T_M3.Subj==s);
    FD_mean = mean(T_M3.FD(idx),'omitnan');
    DV_mean = mean(T_M3.DVARS(idx),'omitnan');
    FD_between(idx) = FD_mean;
    DV_between(idx) = DV_mean;
    FD_within(idx)  = T_M3.FD(idx)    - FD_mean;
    DV_within(idx)  = T_M3.DVARS(idx) - DV_mean;
end

T_M3.FD_within  = FD_within;
T_M3.FD_between = FD_between;
T_M3.DV_within  = DV_within;
T_M3.DV_between = DV_between;

%% ---------- Fit M3 ----------
% FD only
m3_fd = fitlme(T_M3, 'PathLen ~ Phase*Day + FD_within + FD_between + (1|Subj)');
disp(anova(m3_fd, 'DFMethod','Satterthwaite'));
[beta,~,FE] = fixedEffects(m3_fd,'DFMethod','Satterthwaite');
disp(table(m3_fd.CoefficientNames', beta, FE.SE, FE.DF, FE.pValue, ...
    'VariableNames',{'Term','Beta','SE','DF','p'}));

% FD + DVARS (recommended)
m3_both = fitlme(T_M3, ['PathLen ~ Phase*Day + FD_within + FD_between ' ...
                        '+ DV_within + DV_between + (1|Subj)']);
disp(anova(m3_both, 'DFMethod','Satterthwaite'));
[beta,~,FE] = fixedEffects(m3_both,'DFMethod','Satterthwaite');
disp(table(m3_both.CoefficientNames', beta, FE.SE, FE.DF, FE.pValue, ...
    'VariableNames',{'Term','Beta','SE','DF','p'}));

%% plot multiplt regrassion (pathlength, error, dvars)
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



