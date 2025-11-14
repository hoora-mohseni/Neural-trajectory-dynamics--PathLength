%====================================================================================
%% FINAL NARRATIVE ANALYSIS SCRIPT (v_Flexible_Control)
%====================================================================================
%
% This script executes the full 3-step analysis plan + a flexible
% control analysis (Analysis 4).
%
% *** NEW: You can now set the 'num_bins_to_skip_in_control'
%     variable below to test the robustness of Analysis 4. ***
%
% --- 3 Key Metrics (calculated per subject) ---
% 1. Avg_Early_Error: Behavioral performance (avg. of D1 & D2 early learning).
% 2. Avg_Compression: The neural "quench" strategy (avg. of D1 & D2).
% 3. Overall_Neural_Error_Link: The neural "consequence" (r-value from
%    all concatenated task blocks).
%
% --- 3-Step Narrative (run for each neural metric) ---
% 1. (Null) Test: Does Compression predict Performance?
% 2. (Null) Test: Does the Neural Link predict Performance?
% 3. (Discovery) Test: Do the Neural Strategy and Consequence relate?
%
% --- 1 Control Analysis ---
% 4. (Control) Test: Does Strategy predict Consequence *even after
%    removing the overlapping data bins*?
%
% Requires the Statistics and Machine Learning Toolbox.

clear;
clc;
close all;

%====================================================================================
%% --- 1. Load Data ---
%====================================================================================
try
    loaded_data = load('mydata.mat');
    mydata = loaded_data.mydata;
    num_subjects = size(mydata.pathLength.baselineDay1, 1);
    fprintf('Data loaded successfully for %d subjects.\n', num_subjects);
catch ME
    if strcmp(ME.identifier, 'MATLAB:load:couldNotFindFile')
        fprintf('Error: Could not find "Hoora_PathLength_Data.mat".\n');
        fprintf('Please make sure the file is in the same directory as this script.\n');
    else
        fprintf('An error occurred loading the data:\n%s\n', ME.message);
    end
    return;
end

%====================================================================================
%% --- 2. Calculate All Behavioral & Neural Metrics ---
%====================================================================================

fprintf('\n--- Calculating All Metrics for All Analyses ---\n');
epoch_bins = 4; % 40 TRs / 8 TR step size
fprintf('Using epoch_bins = %d (for 40 TRs)\n', epoch_bins);

% *** NEW PARAMETER FOR CONTROL ANALYSIS 4 ***
% Set how many early learning bins to skip for the control analysis.
% Default is epoch_bins (4) to replicate the original finding.
% Try setting this to 8, 10, etc. to test robustness.
num_bins_to_skip_in_control = 10; %setting to 4 skips the first four windows (which would be the 40TR window)
fprintf('Using num_bins_to_skip_in_control = %d\n', num_bins_to_skip_in_control);


% --- Behavioral Metrics ---
% Metric A: Avg_Early_Error
early_err_d1 = mean(mydata.error.learningDay1(:, 1:epoch_bins), 2, 'omitnan');
early_err_d2 = mean(mydata.error.learningDay2(:, 1:epoch_bins), 2, 'omitnan');
avg_early_error = mean([early_err_d1, early_err_d2], 2, 'omitnan');
fprintf('Calculated: Avg_Early_Error\n');

% Metric B: Avg_Learning_Error
all_learning_error = [mydata.error.learningDay1, mydata.error.learningDay2];
avg_learning_error = mean(all_learning_error, 2, 'omitnan');
fprintf('Calculated: Avg_Learning_Error\n');


% --- PathLength Neural Metrics ---
late_bl1_pl = mean(mydata.pathLength.baselineDay1(:, end-epoch_bins+1:end), 2, 'omitnan');
early_l1_pl = mean(mydata.pathLength.learningDay1(:, 1:epoch_bins), 2, 'omitnan');
compression_d1_pl = late_bl1_pl - early_l1_pl;
late_bl2_pl = mean(mydata.pathLength.baselineDay2(:, end-epoch_bins+1:end), 2, 'omitnan');
early_l2_pl = mean(mydata.pathLength.learningDay2(:, 1:epoch_bins), 2, 'omitnan');
compression_d2_pl = late_bl2_pl - early_l2_pl;
average_compression_pl = mean([compression_d1_pl, compression_d2_pl], 2, 'omitnan');
fprintf('Calculated: Avg_Compression (PathLength)\n');

pl_all = [mydata.pathLength.baselineDay1, mydata.pathLength.learningDay1, mydata.pathLength.washoutDay1, ...
          mydata.pathLength.baselineDay2, mydata.pathLength.learningDay2, mydata.pathLength.washoutDay2];
err_all = [mydata.error.baselineDay1, mydata.error.learningDay1, mydata.error.washoutDay1, ...
           mydata.error.baselineDay2, mydata.error.learningDay2, mydata.error.washoutDay2];
overall_pl_error_link_r = NaN(num_subjects, 1);
for i = 1:num_subjects
    [r_val, ~] = corr(pl_all(i, :)', err_all(i, :)', 'rows', 'complete');
    overall_pl_error_link_r(i) = r_val;
end
fprintf('Calculated: Overall_PL_Error_Link\n');

% --- NetDisplacement Neural Metrics ---
late_bl1_disp = mean(mydata.netDisplacement.baselineDay1(:, end-epoch_bins+1:end), 2, 'omitnan');
early_l1_disp = mean(mydata.netDisplacement.learningDay1(:, 1:epoch_bins), 2, 'omitnan');
compression_d1_disp = late_bl1_disp - early_l1_disp;
late_bl2_disp = mean(mydata.netDisplacement.baselineDay2(:, end-epoch_bins+1:end), 2, 'omitnan');
early_l2_disp = mean(mydata.netDisplacement.learningDay2(:, 1:epoch_bins), 2, 'omitnan');
compression_d2_disp = late_bl2_disp - early_l2_disp;
average_compression_disp = mean([compression_d1_disp, compression_d2_disp], 2, 'omitnan');
fprintf('Calculated: Avg_Compression (NetDisplacement)\n');

disp_all = [mydata.netDisplacement.baselineDay1, mydata.netDisplacement.learningDay1, mydata.netDisplacement.washoutDay1, ...
            mydata.netDisplacement.baselineDay2, mydata.netDisplacement.learningDay2, mydata.netDisplacement.washoutDay2];
overall_disp_error_link_r = NaN(num_subjects, 1);
for i = 1:num_subjects
    [r_val, ~] = corr(disp_all(i, :)', err_all(i, :)', 'rows', 'complete');
    overall_disp_error_link_r(i) = r_val;
end
fprintf('Calculated: Overall_Disp_Error_Link\n');

% Store results for summary
results_early_err = struct();
results_avg_err = struct();

%====================================================================================
%% ---        START ANALYSIS PIPELINE (using Avg_Early_Error)          ---
%====================================================================================

fprintf('\n======================================================\n');
fprintf('          RUNNING 3-STEP NARRATIVE (using Avg_Early_Error)\n');
fprintf('======================================================\n');

% --- PathLength ---
fprintf('\n--- PathLength Analyses (vs. Avg_Early_Error) ---\n');
[r1_pl, p1_pl] = corr(average_compression_pl, avg_early_error, 'rows', 'complete');
fprintf('  Analysis 1 (Compression vs. Performance): r=%.3f, p=%.4f\n', r1_pl, p1_pl);
plot_correlation(average_compression_pl, avg_early_error, ...
    'PL Analysis 1: Compression vs. Avg Early Error', ...
    'Average Compression (PL)', 'Average Early Error', ...
    'final_pl_analysis_1_early_error.png', r1_pl, p1_pl);
results_early_err.pl1 = [r1_pl, p1_pl];

[r2_pl, p2_pl] = corr(overall_pl_error_link_r, avg_early_error, 'rows', 'complete');
fprintf('  Analysis 2 (Neural Link vs. Performance): r=%.3f, p=%.4f\n', r2_pl, p2_pl);
plot_correlation(overall_pl_error_link_r, avg_early_error, ...
    'PL Analysis 2: Neural Link vs. Avg Early Error', ...
    'Overall PL-Error Link (r-value)', 'Average Early Error', ...
    'final_pl_analysis_2_early_error.png', r2_pl, p2_pl);
results_early_err.pl2 = [r2_pl, p2_pl];

[r3_pl, p3_pl] = corr(average_compression_pl, overall_pl_error_link_r, 'rows', 'complete');
fprintf('  Analysis 3 (Compression vs. Neural Link): r=%.3f, p=%.4f\n', r3_pl, p3_pl);
plot_correlation(average_compression_pl, overall_pl_error_link_r, ...
    'PL Analysis 3 (Discovery): Neural Strategy vs. Consequence', ...
    'Average Compression (PL)', 'Overall PL-Error Link (r-value)', ...
    'final_pl_analysis_3_discovery.png', r3_pl, p3_pl);
results_early_err.pl3 = [r3_pl, p3_pl];

% --- Analysis 4: CONTROL ANALYSIS (PathLength) ---
fprintf('\n--- PathLength Analysis 4 (Control): Removing Overlap ---\n');
pl_all_modified = [mydata.pathLength.baselineDay1, mydata.pathLength.learningDay1(:, num_bins_to_skip_in_control+1:end), mydata.pathLength.washoutDay1, ...
                   mydata.pathLength.baselineDay2, mydata.pathLength.learningDay2(:, num_bins_to_skip_in_control+1:end), mydata.pathLength.washoutDay2];
err_all_modified = [mydata.error.baselineDay1, mydata.error.learningDay1(:, num_bins_to_skip_in_control+1:end), mydata.error.washoutDay1, ...
                    mydata.error.baselineDay2, mydata.error.learningDay2(:, num_bins_to_skip_in_control+1:end), mydata.error.washoutDay2];
overall_pl_error_link_r_modified = NaN(num_subjects, 1);
for i = 1:num_subjects
    [r_val, ~] = corr(pl_all_modified(i, :)', err_all_modified(i, :)', 'rows', 'complete');
    overall_pl_error_link_r_modified(i) = r_val;
end
fprintf('Calculated: Modified_Overall_PL_Error_Link (skipping %d bins)\n', num_bins_to_skip_in_control);

[r4_pl, p4_pl] = corr(average_compression_pl, overall_pl_error_link_r_modified, 'rows', 'complete');
fprintf('  Correlation (Avg_Compression vs. *Modified* PL-Error Link):\n');
fprintf('  Pearson''s r = %.3f, p = %.4f\n', r4_pl, p4_pl);
results_early_err.pl4 = [r4_pl, p4_pl];


% --- NetDisplacement ---
fprintf('\n--- NetDisplacement Analyses (vs. Avg_Early_Error) ---\n');
[r1_disp, p1_disp] = corr(average_compression_disp, avg_early_error, 'rows', 'complete');
fprintf('  Analysis 1 (Compression vs. Performance): r=%.3f, p=%.4f\n', r1_disp, p1_disp);
plot_correlation(average_compression_disp, avg_early_error, ...
    'NetDisp Analysis 1: Compression vs. Avg Early Error', ...
    'Average Compression (NetDisp)', 'Average Early Error', ...
    'final_disp_analysis_1_early_error.png', r1_disp, p1_disp);
results_early_err.disp1 = [r1_disp, p1_disp];

[r2_disp, p2_disp] = corr(overall_disp_error_link_r, avg_early_error, 'rows', 'complete');
fprintf('  Analysis 2 (Neural Link vs. Performance): r=%.3f, p=%.4f\n', r2_disp, p2_disp);
plot_correlation(overall_disp_error_link_r, avg_early_error, ...
    'NetDisp Analysis 2: Neural Link vs. Avg Early Error', ...
    'Overall Disp-Error Link (r-value)', 'Average Early Error', ...
    'final_disp_analysis_2_early_error.png', r2_disp, p2_disp);
results_early_err.disp2 = [r2_disp, p2_disp];

[r3_disp, p3_disp] = corr(average_compression_disp, overall_disp_error_link_r, 'rows', 'complete');
fprintf('  Analysis 3 (Compression vs. Neural Link): r=%.3f, p=%.4f\n', r3_disp, p3_disp);
plot_correlation(average_compression_disp, overall_disp_error_link_r, ...
    'NetDisp Analysis 3 (Discovery): Neural Strategy vs. Consequence', ...
    'Average Compression (NetDisp)', 'Overall Disp-Error Link (r-value)', ...
    'final_disp_analysis_3_discovery.png', r3_disp, p3_disp);
results_early_err.disp3 = [r3_disp, p3_disp];

% --- Analysis 4: CONTROL ANALYSIS (NetDisplacement) ---
fprintf('\n--- NetDisp Analysis 4 (Control): Removing Overlap ---\n');
disp_all_modified = [mydata.netDisplacement.baselineDay1, mydata.netDisplacement.learningDay1(:, num_bins_to_skip_in_control+1:end), mydata.netDisplacement.washoutDay1, ...
                     mydata.netDisplacement.baselineDay2, mydata.netDisplacement.learningDay2(:, num_bins_to_skip_in_control+1:end), mydata.netDisplacement.washoutDay2];
overall_disp_error_link_r_modified = NaN(num_subjects, 1);
for i = 1:num_subjects
    [r_val, ~] = corr(disp_all_modified(i, :)', err_all_modified(i, :)', 'rows', 'complete');
    overall_disp_error_link_r_modified(i) = r_val;
end
fprintf('Calculated: Modified_Overall_Disp_Error_Link (skipping %d bins)\n', num_bins_to_skip_in_control);

[r4_disp, p4_disp] = corr(average_compression_disp, overall_disp_error_link_r_modified, 'rows', 'complete');
fprintf('  Correlation (Avg_Compression vs. *Modified* Disp-Error Link):\n');
fprintf('  Pearson''s r = %.3f, p = %.4f\n', r4_disp, p4_disp);
results_early_err.disp4 = [r4_disp, p4_disp];

%====================================================================================
%% ---       START ANALYSIS PIPELINE (using Avg_Learning_Error)        ---
%====================================================================================

fprintf('\n======================================================\n');
fprintf('          RUNNING 3-STEP NARRATIVE (using Avg_Learning_Error)\n');
fprintf('======================================================\n');
results_avg_err = struct();

% --- PathLength ---
fprintf('\n--- PathLength Analyses (vs. Avg_Learning_Error) ---\n');
[r1_pl, p1_pl] = corr(average_compression_pl, avg_learning_error, 'rows', 'complete');
fprintf('  Analysis 1 (Compression vs. Performance): r=%.3f, p=%.4f\n', r1_pl, p1_pl);
plot_correlation(average_compression_pl, avg_learning_error, ...
    'PL Analysis 1: Compression vs. Avg Learning Error', ...
    'Average Compression (PL)', 'Average Learning Error (All Bins)', ...
    'final_pl_analysis_1_avg_learning_error.png', r1_pl, p1_pl);
results_avg_err.pl1 = [r1_pl, p1_pl];

[r2_pl, p2_pl] = corr(overall_pl_error_link_r, avg_learning_error, 'rows', 'complete');
fprintf('  Analysis 2 (Neural Link vs. Performance): r=%.3f, p=%.4f\n', r2_pl, p2_pl);
plot_correlation(overall_pl_error_link_r, avg_learning_error, ...
    'PL Analysis 2: Neural Link vs. Avg Learning Error', ...
    'Overall PL-Error Link (r-value)', 'Average Learning Error (All Bins)', ...
    'final_pl_analysis_2_avg_learning_error.png', r2_pl, p2_pl);
results_avg_err.pl2 = [r2_pl, p2_pl];

[r3_pl, p3_pl] = corr(average_compression_pl, overall_pl_error_link_r, 'rows', 'complete');
fprintf('  Analysis 3 (Compression vs. Neural Link): r=%.3f, p=%.4f\n', r3_pl, p3_pl);
copyfile('final_pl_analysis_3_discovery.png', 'final_pl_analysis_3_avg_learning_error.png');
fprintf('Saved plot to ''final_pl_analysis_3_avg_learning_error.png''\n');
results_avg_err.pl3 = [r3_pl, p3_pl];

% --- Analysis 4: CONTROL (PathLength, for this metric) ---
fprintf('\n--- PathLength Analysis 4 (Control): Removing Overlap ---\n');
[r4_pl, p4_pl] = corr(average_compression_pl, overall_pl_error_link_r_modified, 'rows', 'complete');
fprintf('  Correlation (Avg_Compression vs. *Modified* PL-Error Link):\n');
fprintf('  (Using %d skipped bins)\n', num_bins_to_skip_in_control);
fprintf('  Pearson''s r = %.3f, p = %.4f\n', r4_pl, p4_pl);
results_avg_err.pl4 = [r4_pl, p4_pl];


% --- NetDisplacement ---
fprintf('\n--- NetDisplacement Analyses (vs. Avg_Learning_Error) ---\n');
[r1_disp, p1_disp] = corr(average_compression_disp, avg_learning_error, 'rows', 'complete');
fprintf('  Analysis 1 (Compression vs. Performance): r=%.3f, p=%.4f\n', r1_disp, p1_disp);
plot_correlation(average_compression_disp, avg_learning_error, ...
    'NetDisp Analysis 1: Compression vs. Avg Learning Error', ...
    'Average Compression (NetDisp)', 'Average Learning Error (All Bins)', ...
    'final_disp_analysis_1_avg_learning_error.png', r1_disp, p1_disp);
results_avg_err.disp1 = [r1_disp, p1_disp];

[r2_disp, p2_disp] = corr(overall_disp_error_link_r, avg_learning_error, 'rows', 'complete');
fprintf('  Analysis 2 (Neural Link vs. Performance): r=%.3f, p=%.4f\n', r2_disp, p2_disp);
plot_correlation(overall_disp_error_link_r, avg_learning_error, ...
    'NetDisp Analysis 2: Neural Link vs. Avg Learning Error', ...
    'Overall Disp-Error Link (r-value)', 'Average Learning Error (All Bins)', ...
    'final_disp_analysis_2_avg_learning_error.png', r2_disp, p2_disp);
results_avg_err.disp2 = [r2_disp, p2_disp];

[r3_disp, p3_disp] = corr(average_compression_disp, overall_disp_error_link_r, 'rows', 'complete');
fprintf('  Analysis 3 (Compression vs. Neural Link): r=%.3f, p=%.4f\n', r3_disp, p3_disp);
copyfile('final_disp_analysis_3_discovery.png', 'final_disp_analysis_3_avg_learning_error.png');
fprintf('Saved plot to ''final_disp_analysis_3_avg_learning_error.png''\n');
results_avg_err.disp3 = [r3_disp, p3_disp];

% --- Analysis 4: CONTROL (NetDisplacement, for this metric) ---
fprintf('\n--- NetDisp Analysis 4 (Control): Removing Overlap ---\n');
[r4_disp, p4_disp] = corr(average_compression_disp, overall_disp_error_link_r_modified, 'rows', 'complete');
fprintf('  Correlation (Avg_Compression vs. *Modified* Disp-Error Link):\n');
fprintf('  (Using %d skipped bins)\n', num_bins_to_skip_in_control);
fprintf('  Pearson''s r = %.3f, p = %.4f\n', r4_disp, p4_disp);
results_avg_err.disp4 = [r4_disp, p4_disp];


fprintf('\n--- All analyses complete ---\n');

%====================================================================================
%% --- 8. Final Results Summary (for copy/paste) ---
%====================================================================================

fprintf('\n======================================================\n');
fprintf('     FINAL RESULTS SUMMARY (N=%d) - Using Avg_Early_Error\n', num_subjects);
fprintf('======================================================\n');
fprintf('--- PathLength Analyses ---\n');
fprintf('Analysis 1 (Compression vs. Performance):       r = %.3f, p = %.4f\n', results_early_err.pl1(1), results_early_err.pl1(2));
fprintf('Analysis 2 (Neural Link vs. Performance):     r = %.3f, p = %.4f\n', results_early_err.pl2(1), results_early_err.pl2(2));
fprintf('Analysis 3 (Compression vs. Neural Link):     r = %.3f, p = %.4f\n', results_early_err.pl3(1), results_early_err.pl3(2));
fprintf('Analysis 4 (Control: Comp vs. Mod. Link, %d skipped): r = %.3f, p = %.4f\n', num_bins_to_skip_in_control, results_early_err.pl4(1), results_early_err.pl4(2));
fprintf('\n--- NetDisplacement Analyses (Robustness Check) ---\n');
fprintf('Analysis 1 (Compression vs. Performance):       r = %.3f, p = %.4f\n', results_early_err.disp1(1), results_early_err.disp1(2));
fprintf('Analysis 2 (Neural Link vs. Performance):     r = %.3f, p = %.4f\n', results_early_err.disp2(1), results_early_err.disp2(2));
fprintf('Analysis 3 (Compression vs. Neural Link):     r = %.3f, p = %.4f\n', results_early_err.disp3(1), results_early_err.disp3(2));
fprintf('Analysis 4 (Control: Comp vs. Mod. Link, %d skipped): r = %.3f, p = %.4f\n', num_bins_to_skip_in_control, results_early_err.disp4(1), results_early_err.disp4(2));
fprintf('======================================================\n');

fprintf('\n======================================================\n');
fprintf('     FINAL RESULTS SUMMARY (N=%d) - Using Avg_Learning_Error\n', num_subjects);
fprintf('======================================================\n');
fprintf('--- PathLength Analyses ---\n');
fprintf('Analysis 1 (Compression vs. Performance):       r = %.3f, p = %.4f\n', results_avg_err.pl1(1), results_avg_err.pl1(2));
fprintf('Analysis 2 (Neural Link vs. Performance):     r = %.3f, p = %.4f\n', results_avg_err.pl2(1), results_avg_err.pl2(2));
fprintf('Analysis 3 (Compression vs. Neural Link):     r = %.3f, p = %.4f\n', results_avg_err.pl3(1), results_avg_err.pl3(2));
fprintf('Analysis 4 (Control: Comp vs. Mod. Link, %d skipped): r = %.3f, p = %.4f\n', num_bins_to_skip_in_control, results_avg_err.pl4(1), results_avg_err.pl4(2));
fprintf('\n--- NetDisplacement Analyses (Robustness Check) ---\n');
fprintf('Analysis 1 (Compression vs. Performance):       r = %.3f, p = %.4f\n', results_avg_err.disp1(1), results_avg_err.disp1(2));
fprintf('Analysis 2 (Neural Link vs. Performance):     r = %.3f, p = %.4f\n', results_avg_err.disp2(1), results_avg_err.disp2(2));
fprintf('Analysis 3 (Compression vs. Neural Link):     r = %.3f, p = %.4f\n', results_avg_err.disp3(1), results_avg_err.disp3(2));
fprintf('Analysis 4 (Control: Comp vs. Mod. Link, %d skipped): r = %.3f, p = %.4f\n', num_bins_to_skip_in_control, results_avg_err.disp4(1), results_avg_err.disp4(2));
fprintf('======================================================\n');

%====================================================================================
%% --- Helper Function for Plotting (with r/p values) ---
%====================================================================================
function plot_correlation(x_data, y_data, plot_title, x_label_text, y_label_text, filename, r_val, p_val)
    % Creates a scatter plot, adds a regression line, adds stats, and saves
    
    figure('color', 'w', 'Position', [300, 300, 1000, 400]);
    
    % Clean NaNs for plotting and regression
    valid_mask = ~isnan(x_data) & ~isnan(y_data);
    x_clean = x_data(valid_mask);
    y_clean = y_data(valid_mask);
    
    % Scatter plot
    %scatter(x_clean, y_clean, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
    scatter(x_clean, y_clean, 40, [0.1, 0.6, 0.4] + .15*(1-[0.1, 0.6, 0.4]), 'filled'); % 'k' for black
    hold on;
   
    % Add regression line
    coeffs = polyfit(x_clean, y_clean, 1);
    x_fit = linspace(min(x_clean), max(x_clean), 100);
    y_fit = polyval(coeffs, x_fit);
    plot(x_fit, y_fit, 'k-', 'LineWidth', 2);
    viridis(137)
    % --- Add r and p values to the plot ---
    ax = gca;
    xlims = ax.XLim;
    ylims = ax.YLim;
    % Position text
    x_pos_ratio = 0.05; % Default left
    y_pos_ratio = 0.9;  % Default top
    
    % Adjust text position to avoid overlapping with the regression line
    if r_val < 0 && (coeffs(1) < 0) % Negative slope
        x_pos_ratio = 0.60; % Move text to the right
    elseif r_val > 0 && (coeffs(1) > 0) % Positive slope
        y_pos_ratio = 0.20; % Move text to the bottom-left
    elseif r_val < 0 && (coeffs(1) > 0) % Edge case (rare)
        x_pos_ratio = 0.05;
        y_pos_ratio = 0.9;
    elseif r_val > 0 && (coeffs(1) < 0) % Edge case (rare)
        x_pos_ratio = 0.60;
        y_pos_ratio = 0.9;
    end

    x_pos = xlims(1) + x_pos_ratio * (xlims(2) - xlims(1));
    y_pos = ylims(1) + y_pos_ratio * (ylims(2) - ylims(1));
    
    text_str = sprintf('r = %.3f\np = %.4f', r_val, p_val);
    % Added BackgroundColor and EdgeColor to make text pop
    text(x_pos, y_pos, text_str, 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', 'k');
    % --- End of new code ---
    
    % Add labels and title
    title(plot_title);
    xlabel(x_label_text);
    ylabel(y_label_text);
    grid off;
    box off;
    hold off;
    %set(gca, "FontSize", 20, 'LineWidth', .8, 'TickDir', 'out')
    axis square;

    % Save the figure
    saveas(gcf, filename);
    fprintf('Saved plot to ''%s''\n', filename);
end





