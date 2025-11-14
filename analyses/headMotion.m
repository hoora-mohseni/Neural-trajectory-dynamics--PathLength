clc
clear
close all
%====================================================================================
%% Extracting HeadMotion(framewise_displacement) and Dvars
%====================================================================================

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

%====================================================================================
%% save
%====================================================================================
save('hm_data.mat', 'hm', '-v7.3');  % -v7.3 if it’s big
save('dvars_data.mat', 'dvars', '-v7.3');  % -v7.3 if it’s big

%====================================================================================
%% load
%====================================================================================
load('hm_data.mat', 'hm');
load('dvars_data.mat', 'dvars');

%====================================================================================
%% 
%====================================================================================

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

%====================================================================================
%% Sliding wondows
%====================================================================================
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
%% Fit M0: FD ~ Phase * Day + (1 | Subject)
%====================================================================================

lme_M0 = fitlme(T, 'FD ~ Phase*Day + (1|Subject)');  % random intercept per subject

% Inference on fixed effects (Phase, Day, Phase×Day):
anova(lme_M0, 'DFMethod','Satterthwaite')

%====================================================================================
%% Get adjusted means (Phase × Day) for plotting
%====================================================================================

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

%====================================================================================
%% ---------- Build table: one row per (Subject, Phase, Day, Window) ----------
%====================================================================================

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

