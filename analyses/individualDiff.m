
DPL_LB_EL1 = mean(pathLength.early.(epochs{3}), 2, 'omitnan') - mean(pathLength.late.(epochs{2}), 2, 'omitnan');
DPL_EL_LL1 = mean(pathLength.late.(epochs{3}), 2, 'omitnan')  - mean(pathLength.early.(epochs{3}), 2, 'omitnan');
DPL_LL_EW1 = mean(pathLength.early.(epochs{4}), 2, 'omitnan') - mean(pathLength.late.(epochs{3}), 2, 'omitnan');
DPL_LB_EL2 = mean(pathLength.early.(epochs{7}), 2, 'omitnan') - mean(pathLength.late.(epochs{6}), 2, 'omitnan');
DPL_EL_LL2 = mean(pathLength.late.(epochs{7}), 2, 'omitnan')  - mean(pathLength.early.(epochs{7}), 2, 'omitnan');
DPL_LL_EW2 = mean(pathLength.early.(epochs{8}), 2, 'omitnan') - mean(pathLength.late.(epochs{7}), 2, 'omitnan');


Derr_LB_EL1 = mean(error.(epochs{3})(: ,1:4), 2, 'omitnan');
Derr_EL_LL1 = mean(error.(epochs{3})(: ,end-3:end), 2, 'omitnan');
Derr_LL_EW1 = mean(error.(epochs{4})(: ,1:4), 2, 'omitnan');
Derr_LB_EL2 = mean(error.(epochs{7})(: ,1:4), 2, 'omitnan');
Derr_EL_LL2 = mean(error.(epochs{7})(: ,end-3:end), 2, 'omitnan');
Derr_LL_EW2 = mean(error.(epochs{8})(: ,1:4), 2, 'omitnan');

%%
% List of variable pairs
DPL_vars = {DPL_LB_EL1, DPL_EL_LL1, DPL_LB_EL2, DPL_EL_LL2};
Derr_vars = {Derr_LB_EL1, Derr_EL_LL1, Derr_LB_EL2, Derr_EL_LL2};
pairNames = {'EL-LB', 'LL-EL', 'EL-LB', 'LL-EL'};

% Create figure with subplots
figure('Color', 'w', 'Position', [400, 400, 800, 500]);
for i = 1:length(DPL_vars)
    nexttile;
    
    % Flatten data (subjects × windows → vector)
    x = DPL_vars{i}(:);
    y = Derr_vars{i}(:);
    
    if i < 3
        scatter(x, y, 40, red_light + .15*(1-red_light), 'filled'); hold on
    else
        scatter(x, y, 40, blue_light + .05*(1-blue_light), 'filled'); hold on
    end

    % Fit linear regression
    p = polyfit(x, y, 1);   % p(1) = slope, p(2) = intercept
    yfit = polyval(p, x);
    
    % Plot regression line in red
    plot(x, yfit, 'k-' , 'LineWidth', 2);

    % Compute correlation
    [r, pval] = corr(x(:), y(:), 'Type', 'Pearson');

    % Add text with r and p on figure
    txt = sprintf('r = %.3f, p = %.3g', r, pval);
    xpos = min(x) + 0.05*(max(x)-min(x));
    ypos = max(y) - 0.1*(max(y)-min(y));
    %text(xpos, ypos, txt, 'FontSize', 12, 'Color', 'r');


    % xlabel('Δ Path Length');
    % ylabel('Error');
    % title(pairNames{i});
    grid on;
    
    % Add regression line
    hold on;
    coeffs = polyfit(x, y, 1);
    xfit = linspace(min(x), max(x), 100);
    yfit = polyval(coeffs, xfit);
    plot(xfit, yfit, 'k-', 'LineWidth', 1.5);
    hold off;
    
    % Show correlation
    R = corrcoef(x, y);
    if numel(R) > 1
        text(mean(x), mean(y), sprintf('r = %.2f', R(1,2)), ...
            'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w');
    end
end

% sgtitle('Scatter plots: ΔPath Length vs Error');


%%


figure('Color', 'w', 'Position', [400, 400, 500, 500]),
x = 1:2;

DPL_vars = [DPL_LB_EL1, DPL_EL_LL1, DPL_LB_EL2, DPL_EL_LL2];
bar(mean(DPL_vars(:, 1:2)),'FaceColor',red_light + .15*(1-red_light),'EdgeColor',red_light,'LineWidth',1, 'BarWidth', 0.4); hold on
bar(x+2, mean(DPL_vars(:, 3:4)),'FaceColor',blue_light + .06*(1-blue_light),'EdgeColor',blue_light,'LineWidth',1, 'BarWidth', 0.5); hold on

for i = 1:4
    errorbar(i+.1, mean(DPL_vars(:, i)), std(DPL_vars(:, i))/sqrt(32), 'k', 'LineWidth', 1)
    plot((i-.1)*ones(1, 32), DPL_vars(:, i), 'k.', 'MarkerSize', 10)
end

set(gca, 'FontSize', 15); box off
xticks(1:4);

% xticklabels(pairNames);

% xlabel('Contrast');
% ylabel('Δ Path Length (DPL)');

xline(2.5, '--')
xlim([0.5, 4.5])
yticks(-60:10:60); ylim([-60, 60])
box off
% text(1.5, -40, 'Day 1', 'HorizontalAlignment','center', 'FontSize', 20)
% text(3.5, -40, 'Day 2', 'HorizontalAlignment','center', 'FontSize', 20)


%%


figure('Color', 'w', 'Position', [300, 500, 1000, 700]);

x = 1:2;
Derr_vars = [Derr_LB_EL1, Derr_EL_LL1, Derr_LB_EL2, Derr_EL_LL2];
bar(mean(Derr_vars(:, 1:2)),'FaceColor',red_light + .5*(1-red_light),'EdgeColor',red_light,'LineWidth',1.5, 'BarWidth', 0.5); hold on
bar(x+2, mean(Derr_vars(:, 3:4)),'FaceColor',blue_light + .5*(1-blue_light),'EdgeColor',blue_light,'LineWidth',1.5, 'BarWidth', 0.5); hold on
for i = 1:4
    errorbar(i+.1, mean(Derr_vars(:, i)), std(Derr_vars(:, i))/sqrt(32), 'k', 'LineWidth',0.5)
    plot((i-.1)*ones(1, 32), Derr_vars(:, i), 'k.' , 'MarkerSize', 10)
end

set(gca, 'FontSize', 15); box off
xticks(1:4);
xticklabels(pairNames);

xlabel('Contrast');
ylabel('Δ Error (Derr)');

xline(2.5, '--')
text(2, 40, 'Day 1', 'HorizontalAlignment','center', 'FontSize', 20)
text(4, 40, 'Day 2', 'HorizontalAlignment','center', 'FontSize', 20)





%%  Mixed-effects models (population effects + individual differences)

% --- mapping ---
epochs  = {'restDay1','baselineDay1','learningDay1','washoutDay1', ...
           'restDay2','baselineDay2','learningDay2','washoutDay2'};
phaseOf = containers.Map(epochs, {'rest','baseline','learning','washout', ...
                                  'rest','baseline','learning','washout'});
dayOf   = containers.Map(epochs, {'day1','day1','day1','day1', ...
                                  'day2','day2','day2','day2'});

% --- gather rows: one row per (Subject, Phase, Day, Window) ---
PL = [];  Subject = [];  Phase = {};  Day = {};  Window = [];

numSubjects = numel(subjects);
for e = 1:numel(epochs)
    %M = hmSliding.(epochs{e});          % [numSubjects x numWindows_e]
    M = pathLengthSliding.(epochs{e}); 
    [nS, nW] = size(M);

    PL      = [PL;      M(:)];
    Subject = [Subject; repelem((1:nS)', nW, 1)];
    Window  = [Window;  repmat((1:nW)', nS, 1)];
    Phase   = [Phase;   repmat({phaseOf(epochs{e})}, nS*nW, 1)];
    Day     = [Day;     repmat({dayOf(epochs{e})},   nS*nW, 1)];
end

T = table;
T.PL      = PL;
T.Subject = categorical(Subject);
T.Phase   = categorical(Phase, {'rest','baseline','learning','washout'}); % set reference order
T.Day     = categorical(Day,   {'day1','day2'});
T.Window  = Window;   % (not used in the fixed part of the model)

lmePL1  = fitlme(T, 'PL ~ Day*Phase + (1 + Phase | Subject)', 'FitMethod','REML');
lmePL2  = fitlme(T, 'PL ~ Day*Phase + (1 | Subject)', 'FitMethod','REML');

cmp = compare(lmePL2, lmePL1,'CheckNesting',true);   % likelihood-ratio test
disp(cmp)

%lmeErr = fitlme(error,  'logErr ~ day*phase + speed + (1 + phase | subject)');

disp(anova(lmePL1));   % fixed-effect tests
%disp(anova(lmeErr));

% Extract subject BLUPs (random effects) for individual differences
rePL  = randomEffects(lmePL1);   % table with levels for each subject/coef
%reErr = randomEffects(lmeErr);

rePL = reshape(rePL, 32, 4);

figure
bar(mean(rePL)); hold on
for i = 1:4
    errorbar(i+.1, mean(rePL(:, i)), std(rePL(:, i))/sqrt(32), 'k', 'LineWidth', 1.5)
    plot((i-.1)*ones(1, 32), rePL(:, i), 'k.', 'MarkerSize', 10)
end


%% ---------- Prep: turn your structs into a long table T ----------
% Required input:
%   error:      1x1 struct with fields baselineDay1, learningDay1, washoutDay1,
%               baselineDay2, learningDay2, washoutDay2; each is S x T(double)
%   pathlength: same fields/shapes as `error`

pathlength = rmfield(pathLengthSliding, {'restDay1', 'restDay2'});

assert(isstruct(error) && isstruct(pathlength), 'Provide `error` and `pathlength` structs.');

ef = fieldnames(error);
pf = fieldnames(pathlength);
assert(isequal(sort(ef), sort(pf)), 'error and pathlength must have the same fields.');

% Helper to parse "baselineDay1" -> phase="baseline", day=1
parseField = @(f) deal( ...
    regexp(f, '^(baseline|learning|washout)', 'tokens','once'), ...
    regexp(f, 'Day(\d+)$', 'tokens','once'));

all_subject   = [];
all_day       = [];
all_phase     = {};
all_trial     = [];
all_error     = [];
all_pathlen   = [];

for i = 1:numel(ef)
    f  = ef{i};
    [phaseTok, dayTok] = parseField(f);
    assert(~isempty(phaseTok) && ~isempty(dayTok), ...
        'Field name "%s" must look like baselineDay1/learningDay2/etc.', f);

    phase = string(phaseTok{1});
    day   = str2double(dayTok{1});

    E = error.(f);           % S x T
    PL = pathlength.(f);     % S x T
    assert(ismatrix(E) && isnumeric(E) && isequal(size(E), size(PL)), ...
        'Sizes of error.%s and pathlength.%s must match.', f, f);

    [S, T] = size(E);

    % Build long vectors for this field
    subj   = repelem((1:S).', T, 1);      % S*T x 1
    triali = repmat((1:T).', S, 1);       % S*T x 1
    dayv   = repmat(day, S*T, 1);
    phasev = repmat(phase, S*T, 1);

    all_subject = [all_subject; subj];
    all_day     = [all_day;     dayv];
    all_phase   = [all_phase;   cellstr(phasev)];
    all_trial   = [all_trial;   triali];
    all_error   = [all_error;   E(:)];
    all_pathlen = [all_pathlen; PL(:)];
end

T = table;
T.subject   = categorical(all_subject);                    % subject ID as categorical
T.day       = categorical(all_day, [1 2], {'day1','day2'}); % ordered levels
T.phase     = categorical(all_phase, {'baseline','learning','washout'}); % ordered levels
T.trial     = all_trial;                                   % within-epoch window index
T.error     = all_error;
T.PL        = all_pathlen;

% Optional: drop missing rows if present
T = rmmissing(T);

% Quick sanity checks
% head(T)
% summary(T)

%% ---------- Mixed-effects models (no log, no speed) ----------
% Model: y ~ day*phase + (1 + phase | subject)
% Random effects: subject-specific intercepts and phase slopes

lmePL  = fitlme(T,  'PL    ~ day*phase + (1 + phase | subject)',  'FitMethod','REML');
lmeErr = fitlme(T,  'error ~ day*phase + (1 + phase | subject)',  'FitMethod','REML');

%% ---------- Results ----------
disp('Fixed-effect tests (ANOVA) for PL:');
disp(anova(lmePL));

disp('Fixed-effect tests (ANOVA) for error:');
disp(anova(lmeErr));

disp('Fixed-effect estimates for PL:');
disp(lmePL.Coefficients);

disp('Fixed-effect estimates for error:');
disp(lmeErr.Coefficients);

% Subject-level random effects (BLUPs)
rePL  = randomEffects(lmePL);   % table of subject random intercept/slope adjustments
reErr = randomEffects(lmeErr);

% Example: extract random slopes for phase-learning vs baseline
% (names depend on internal coding; inspect once)
% rePL(strcmp(rePL.Name,'phase_learning'), :)
% reErr(strcmp(reErr.Name,'phase_learning'), :)
