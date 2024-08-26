%% code is to 1) Pair Images with Corresponding Memorability Values
%% 2) Perform Comparisons 
clear;clc;

% Path
respath = '/Users/zhuang/Documents/Projects/Project2_MEN/MEM_MEG_plus_exp/results/';
dspath = '/Users/zhuang/Documents/Projects/Project2_MEN/ratings/';

% Initialize arrays to hold the results
subject_ids = [];
correlation_coeffs = [];
p_values = [];
remembered_rate_memorable_all=[];
remembered_rate_forgettable_all=[];
mean_memorable_confidence_all=[];
mean_forgettable_confidence_all=[];

% loop across subjects
for i=1:4
    datafile = fullfile([respath,sprintf('merged_trials_sub0%d.mat',i)]);
    load(datafile)
    % Extract the image numbers into an array
    category_numbers = [merged_trials.category_nr];

    % Extract image paths from merged trials data
    image_paths = {merged_trials.image_path}; % Assuming image_path is a string or cell array of strings
    
    seen=[merged_trials.seen];
    trial_type={merged_trials.trial_type};
    image_nr=[merged_trials.image_nr];
    category_nr=[merged_trials.category_nr];

    confidence_values = [merged_trials.confidence];
    
    % Remove the 'images_meg/' prefix from image_paths
    prefix_to_remove_img = 'images_meg/';
    image_paths = cellfun(@(p) strrep(p, prefix_to_remove_img, ''), image_paths, 'UniformOutput', false);

    % Load the memorability scores (CR values) from the CSV file
    memorability_data = readtable([dspath,'THINGS_Memorability_Scores.csv']);

    % Extract image paths and CR values from memorability data
    memorability_image_paths = memorability_data.file_path; % Assuming image_path is a string or cell array of strings
    cr_values = memorability_data.cr;
    
    % Remove the 'images/' prefix from memorability_image_paths
    prefix_to_remove = 'images/';
    memorability_image_paths = cellfun(@(p) strrep(p, prefix_to_remove, ''), memorability_image_paths, 'UniformOutput', false);

    % Initialize an array to hold matched CR values
    matched_cr_values = zeros(size(confidence_values));

    % Loop through each trial and find the corresponding CR value
    for j = 1:length(image_paths)
        % Find the row in memorability_data corresponding to the current image_path
        idx = find(strcmp(memorability_image_paths, image_paths{j}));
    
        if ~isempty(idx)
            matched_cr_values(j) = cr_values(idx);
        end
    end

    % Create a table with the relevant data
    results_table = table(image_paths', confidence_values', seen', trial_type', image_nr', category_nr', matched_cr_values', ...
                          'VariableNames', {'image_paths', 'confidence', 'seen', 'TrialType', 'ImageNr', 'CategoryNr', 'MatchedCRValues'});

    % Save the results table to a .mat file
    save(fullfile(respath, sprintf('results_sub0%d.mat', i)), 'results_table');
    
    %% To determine whether people remember memorable items more than forgettable items, 
    % you need to analyze whether the CR values (memorability scores) correlate with 
    % the probability of an item being remembered.
    % Define a threshold for memorable vs. forgettable (e.g., CR value > median CR)
    cr_threshold = median(matched_cr_values(~isnan(matched_cr_values)));
   
    % Determine whether items are memorable or forgettable
    is_memorable = matched_cr_values > cr_threshold;

    % Compute remembering rates
    remembered_rate_memorable = nanmean(seen(is_memorable));
    remembered_rate_forgettable = nanmean(seen(~is_memorable & ~isnan(matched_cr_values)));
    
    % store the results
    remembered_rate_memorable_all=[remembered_rate_memorable_all;remembered_rate_memorable];
    remembered_rate_forgettable_all=[remembered_rate_forgettable_all;remembered_rate_forgettable];

    % Display results
    fprintf('Subject %d:\n', i);
    fprintf('Remembered Rate for Memorable Items: %.2f\n', remembered_rate_memorable);
    fprintf('Remembered Rate for Forgettable Items: %.2f\n', remembered_rate_forgettable);

    % Optionally, perform a statistical test to compare the proportions
    [h, p] = prop_test([nansum(seen(is_memorable)), nansum(seen(~is_memorable & ~isnan(matched_cr_values)))], ...
                        [nansum(is_memorable), nansum(~is_memorable & ~isnan(matched_cr_values))]);
    fprintf('Proportion Test result: h = %d, p = %.4f\n', h, p);
    %% Analysis: Compare confidence for memorable vs. forgettable items 
    % Separate the confidence values based on CR values
    memorable_confidence = confidence_values(matched_cr_values > cr_threshold);
    forgettable_confidence = confidence_values(matched_cr_values <= cr_threshold & ~isnan(matched_cr_values));
    
    % Calculate mean confidence values
    mean_memorable_confidence = nanmean(memorable_confidence);
    mean_forgettable_confidence = nanmean(forgettable_confidence);
    
    % store the results
    mean_memorable_confidence_all=[mean_memorable_confidence_all;mean_memorable_confidence];
    mean_forgettable_confidence_all=[mean_forgettable_confidence_all;mean_forgettable_confidence];

    % Display results
    %fprintf('Subject %d:\n', i);
    fprintf('Mean Confidence for Memorable Items: %.2f\n', mean_memorable_confidence);
    fprintf('Mean Confidence for Forgettable Items: %.2f\n', mean_forgettable_confidence);
    
    % Optionally, perform a statistical test to compare means
    [h, p] = ttest2(memorable_confidence, forgettable_confidence);
    fprintf('T-test result: h = %d, p = %.4f\n', h, p);

    %% correlation
    % Remove NaN values for correlation calculation
    valid_indices = ~isnan(confidence_values);
    valid_confidence = confidence_values(valid_indices);
    valid_cr_values = matched_cr_values(valid_indices);

    % Calculate Pearson correlation coefficient
    [R, P] = corr(valid_confidence', valid_cr_values', 'Type', 'Pearson');

    % Store the results
    subject_ids = [subject_ids; i];
    correlation_coeffs = [correlation_coeffs; R];
    p_values = [p_values; P];

    % Display results
    %fprintf('Subject %d:\n', i);
    fprintf('Pearson Correlation Coefficient between Confidence and Memorability: %.2f\n', R);
    fprintf('P-value of the correlation: %.4f\n', P);

end
% Create a summary table including additional metrics
summary_table = table(subject_ids, correlation_coeffs, p_values, ...
                      remembered_rate_memorable_all, remembered_rate_forgettable_all, ...
                      mean_memorable_confidence_all, mean_forgettable_confidence_all, ...
                      'VariableNames', {'SubjectID', 'CorrelationBetween Confidence and Memorability', 'PValue', ...
                                         'RememberedRateMemorable', 'RememberedRateForgettable', ...
                                         'MeanMemorableConfidence', 'MeanForgettableConfidence'});

% Display the summary table
disp(summary_table);

% Save the summary table to a .mat file
save(fullfile(respath, 'correlation_summary.mat'), 'summary_table');

function [h, p] = prop_test(x, n)
    % Perform a two-sample test for proportions
    p1 = x(1) / n(1);
    p2 = x(2) / n(2);
    p_pool = (x(1) + x(2)) / (n(1) + n(2));
    se = sqrt(p_pool * (1 - p_pool) * (1 / n(1) + 1 / n(2)));
    z = (p1 - p2) / se;
    p = 2 * (1 - normcdf(abs(z))); % Two-tailed p-value
    h = p < 0.05; % Hypothesis test result (0: fail to reject H0, 1: reject H0)
end