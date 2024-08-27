clear;clc;

% Path
respath = '/Users/zhuang/Documents/Projects/Project2_MEN/MEM_MEG_plus_exp/results/';
dspath = '/Users/zhuang/Documents/Projects/Project2_MEN/ratings/';

% Initialize arrays to hold the results
subject_ids = [];
number_remembered_items = [];
number_forgotten_items = [];
mean_remembered_cr_array = [];
mean_forgotten_cr_array = [];
ttest_p_values = [];
correlation_coefficients = [];
correlation_p_values = [];
confidence_correlation_coefficients = [];
confidence_correlation_p_values = [];

% loop across subjects
for i = 1:4
    datafile = fullfile([respath, sprintf('merged_trials_sub0%d.mat', i)]);
    load(datafile)

    % Extract image paths, seen values, and original confidence values
    image_paths = {merged_trials.image_path};
    seen = [merged_trials.seen];
    confidence_values = [merged_trials.confidence];
    trial_type = {merged_trials.trial_type};

    % Assign new confidence values based on seen and unseen conditions
    assigned_confidence_values = zeros(size(confidence_values)); % Initialize assigned confidence array
    
    % Assign confidence values according to the conditions
    assigned_confidence_values(seen == 1 & confidence_values == 2) = 1.5;    % High confidence seen
    assigned_confidence_values(seen == 1 & confidence_values == 1) = 0.5;     % Low confidence seen
    assigned_confidence_values(seen == 0 & confidence_values == 1) = -0.5;    % Low confidence not seen
    assigned_confidence_values(seen == 0 & confidence_values == 2) = -1.5;   % High confidence not seen

    % Remove the 'images_meg/' prefix from image_paths
    prefix_to_remove_img = 'images_meg/';
    image_paths = cellfun(@(p) strrep(p, prefix_to_remove_img, ''), image_paths, 'UniformOutput', false);

    % Load memorability scores (CR values)
    memorability_data = readtable([dspath, 'THINGS_Memorability_Scores.csv']);
    memorability_image_paths = memorability_data.file_path;
    cr_values = memorability_data.cr;

    % Remove the 'images/' prefix from memorability_image_paths
    prefix_to_remove = 'images/';
    memorability_image_paths = cellfun(@(p) strrep(p, prefix_to_remove, ''), memorability_image_paths, 'UniformOutput', false);

    % Match CR values to the trials
    matched_cr_values = zeros(size(assigned_confidence_values));
    for j = 1:length(image_paths)
        idx = find(strcmp(memorability_image_paths, image_paths{j}));
        if ~isempty(idx)
            matched_cr_values(j) = cr_values(idx);
        end
    end
    
    % Identify remembered and forgotten items based on seen==1 and trial_type=='target'
    remembered = (seen == 1) & strcmp(trial_type, 'target');
    forgotten = (seen == 0) & strcmp(trial_type, 'target');
    
    % memorability values of remembered/ forgotten
    remembered_cr=matched_cr_values(remembered);
    forgotten_cr=matched_cr_values(forgotten);
    % mean
    mean_remembered_cr=mean(remembered_cr);
    mean_forgotten_cr=mean(forgotten_cr);
    mean_remembered_cr_array = [mean_remembered_cr_array; mean_remembered_cr];
    mean_forgotten_cr_array = [mean_forgotten_cr_array; mean_forgotten_cr];
    % t-test
    [h, p] = ttest2(remembered_cr, forgotten_cr);
    ttest_p_values = [ttest_p_values; p];

    % Combine remembered and forgotten into a single binary variable
    remembered_forgotten = remembered - forgotten;

    % sum remembered to caculate how many items has been remembered when
    % being seen
    number_remembered = sum(remembered);
    number_remembered_items = [number_remembered_items; number_remembered];
    % sum remembered to caculate how many items has been forgotten when
    % being seen
    number_forgotten = sum(forgotten);
    number_forgotten_items = [number_forgotten_items; number_forgotten];

    % remove 0 [0 refers to 'foil']
    valid_remembered_forgotten=remembered_forgotten;
    valid_remembered_forgotten(valid_remembered_forgotten == 0) = [];
    
    % remembered =1; forgotten =0;
    valid_remembered_forgotten(valid_remembered_forgotten==-1)=0;

    remembered_forgotten(remembered_forgotten==-1)=1;
    logicalremembered_forgotten = remembered_forgotten ~= 0;
    
    % Run the correlation between the remembered/forgotten binary variable and memorability
    valid_cr_values = matched_cr_values(logicalremembered_forgotten);

    [R_cr, P_cr] = corr(valid_remembered_forgotten', valid_cr_values', 'Type', 'Pearson');
    correlation_coefficients = [correlation_coefficients; R_cr];
    correlation_p_values = [correlation_p_values; P_cr];

    % Display
    fprintf('Subject %d:\n', i);
    fprintf('numbers of remembered items: %.2f\n', number_remembered);
    fprintf('numbers of forgotten items: %.2f\n', number_forgotten);
    fprintf('average of memorability for remembered items: %.2f\n', mean_remembered_cr);
    fprintf('average of memorability for of forgotten items: %.2f\n', mean_forgotten_cr);
    fprintf('difference of memorability between remembered and forgotten (t-test): %.2f\n', p);

    fprintf('Pearson Correlation between remembered (1) / forgotten (0) and memorability: %.2f\n', R_cr);
    fprintf('P-value of the correlation between remembered / forgotten and memorability: %.4f\n', P_cr);

    %% Calculate Pearson correlation coefficient between assigned confidence and memorability
    valid_indices = ~isnan(assigned_confidence_values);
    valid_assigned_confidence = assigned_confidence_values(valid_indices);
    valid_cr_values = matched_cr_values(valid_indices);
    [R, P] = corr(valid_assigned_confidence', valid_cr_values', 'Type', 'Pearson');
    confidence_correlation_coefficients = [confidence_correlation_coefficients; R];
    confidence_correlation_p_values = [confidence_correlation_p_values; P];

    % Store subject IDs
    subject_ids = [subject_ids; i];

    % Display results
    fprintf('Pearson Correlation between the confidence(1.5,0.5,-0.5,-1.5) and memorability: %.2f\n', R);
    fprintf('P-value of the correlation between the confidence and memorability: %.4f\n', P);

    
end

% Create a table to save the results
results_table = table(subject_ids, number_remembered_items, number_forgotten_items, ...
                      mean_remembered_cr_array, mean_forgotten_cr_array, ...
                      ttest_p_values, correlation_coefficients, correlation_p_values, ...
                      confidence_correlation_coefficients, confidence_correlation_p_values,...
                      'VariableNames', {'SubjectID', 'number_remembered_items', 'number_forgotten_items', ...
                                         'avg_memorability_of_remembered', 'avg_memorability_of_forgotten', ...
                                         'P values of difference_mem_between_remembered&forgotten', ...
                                         'r between remembered (1) / forgotten (0) and memorability',...
                                         'P of the correlation','r between the confidence(1.5,0.5,-0.5,-1.5) and memorability',...
                                         'P of the correlation between the confidence and memorability'});
results_table1 = table(subject_ids, number_remembered_items, number_forgotten_items, ...
                      mean_remembered_cr_array, mean_forgotten_cr_array, ...
                      ttest_p_values,...
                      'VariableNames', {'SubjectID', 'number_remembered_items', 'number_forgotten_items', ...
                                         'avg_memorability_of_remembered', 'avg_memorability_of_forgotten', ...
                                         'P values of difference_mem_between_remembered&forgotten'});
% Create a table to save the results
results_table2 = table(subject_ids, correlation_coefficients, correlation_p_values, ...
                      confidence_correlation_coefficients, confidence_correlation_p_values,...
                      'VariableNames', {'SubjectID',...
                                         'r between remembered (1) / forgotten (0) and memorability',...
                                         'P of the correlation','r between the confidence(1.5,0.5,-0.5,-1.5) and memorability',...
                                         'P of the correlation between the confidence and memorability'});


display(results_table1)
display(results_table2)

% Save the summary table to a .mat file
save(fullfile(respath, 'results_table.mat'), 'results_table');
