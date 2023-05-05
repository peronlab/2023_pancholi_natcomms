%Script to plot Figure 4f: population dynamics during training
%
%To use this script, first set path to session_response_data.mat file
%
%Next, define which ROIs you want to use (set roi_selection)
%   - 1: Opsin-expressing (O+)
%   - 2: Opsin non-expressing (O-)
%   - 3: Opsin-ambiguous (O?)
%
%Next, define your response pool threshold. This is the percentile threshold for responses (i.e. 0 = all ROIs, 95 = top 5%, 99 = top 1%).

%% Parameters
[all_dat settings] = get_all_data_pstim; 
response_data_path = [settings.base_dir filesep];
roi_selection = 1; %Which ROIs to plot. (1) Opsin-expressing ROIs (O+), (2) Opsin non-expressing ROIs (O-), (3) Opsin-ambiguous ROIs (O?)
response_pool_threshold = 0; %Percentile threshold for response. ROIs with response > percentile response threshold will be included.

%% Get data

%Load session response data
try
    load([response_data_path, 'session_response_data.mat'])
catch
    error('Error loading session response data. Check file path and try again.')
end

%Select stim value to observe across sessions
stim_value_across_sessions = 9; %Will get number of ROIs responding to this stim value across sessions
selected_stim_index = find(possible_stim_values == stim_value_across_sessions, 1); %Get selected stim value

%Get ROIs to use
roi_logicals = cell(1, length(animal_IDs)); %Preallocate
for animal_index = 1:length(animal_IDs) %Loop over animals
    switch roi_selection
        case 1 %Red
            roi_logicals{animal_index} = ismember(response_struct.all_roi_IDs{animal_index}, response_struct.red_roi_IDs{animal_index});
        case 2 %Green
            roi_logicals{animal_index} = ismember(response_struct.all_roi_IDs{animal_index}, response_struct.green_roi_IDs{animal_index});
        case 3 %Ambig
            roi_logicals{animal_index} = ismember(response_struct.all_roi_IDs{animal_index}, response_struct.ambig_roi_IDs{animal_index});        
        otherwise
            error('Invalid value for roi_selection. Use either 1 for opsin-expressing, 2 for opsin non-expressing, and 3 for opsin-ambiguous ROIs.')
    end
end

%% Plot 

%Figure parameters
animal_coloring_scheme = {[255, 31, 91]/255, [0, 154, 222]/255, [175, 88, 186]/255, [255, 198, 30]/255, [242, 133, 34]/255, [0 255 0]/255}; %Colors for different animals
figure_position = [100, 300, 1500, 500]; %[left bottom width height] in pixels. Value vector for figure 'Position' property to control figure window position and size.

%% Z-scored DFF (top row in panel)

%Plot
figure
summary_matrix = NaN(3, length(animal_IDs)); %Matrix to hold summary data fpr ROI responses
for animal_index = 1:length(animal_IDs)
    session_groups = {1:3, round(length(response_struct.session_date{animal_index})/2)-1:round(length(response_struct.session_date{animal_index})/2)+1, length(response_struct.session_date{animal_index})-2:length(response_struct.session_date{animal_index})}; %Session grouping indices for early, mid, and late sessions
    response_matrix = response_struct.normalized_response{animal_index, selected_stim_index};
    response_threshold = prctile(response_matrix, response_pool_threshold, 'all'); %Get response threshold
    response_logicals = nansum(response_matrix > response_threshold, 2) > 0; %Logicals for ROIs with at least one session with response > response threshold
    response_logicals = response_logicals & roi_logicals{animal_index}; %Apply ROI logicals
    response_matrix = response_matrix(response_logicals, :); %Apply response threshold
    for group_index = 1:length(session_groups)
        summary_matrix(group_index, animal_index) = nanmedian(nanmean(response_matrix(:, session_groups{group_index}), 2)); %Get median response over ROIs in each session group
    end
    plot(summary_matrix(:, animal_index), 'o', 'Color', animal_coloring_scheme{animal_index}, 'MarkerSize', 6) %Plot marker for each animal
    hold on
    plot(summary_matrix(:, animal_index), '-', 'Color', animal_coloring_scheme{animal_index}, 'LineWidth', 2) %Plot line for each animal
end
errorbar(nanmean(summary_matrix, 2), nanstd(summary_matrix, 0, 2)/sqrt(length(animal_IDs)), 'o', 'Color', 'k', 'MarkerSize', 10, 'CapSize', 0) %Plot grand mean
xlim([0.5 3.5]); ylim([0 8])
xticks([1 2 3]); xticklabels({'Early', 'Middle', 'Late'});
xlabel('Sessions'); ylabel('Median z-scored \DeltaF/F')
title({['z-scored \DeltaF/F to ', num2str(stim_value_across_sessions), ' stimulus pulses']})
set(gcf, 'Position', [300, 300, 350, 200]) %Set figure window size to be consistent

%Stats
[~, t_test_p_em] = ttest(summary_matrix(1, :), summary_matrix(2, :)); %Paired t-test between early and middle
[~, t_test_p_ml] = ttest(summary_matrix(2, :), summary_matrix(3, :)); %Paired t-test between middle and late
corrected_p_DFF = [t_test_p_em, t_test_p_ml] * 2; %Paired t-tests corrected for multiple comparisons (Bonferroni)

%% Probability of response (middle row in panel)

%Plot
figure
summary_matrix = NaN(3, length(animal_IDs)); %Matrix to hold summary data fpr ROI responses
for animal_index = 1:length(animal_IDs)
    session_groups = {1:3, round(length(response_struct.session_date{animal_index})/2)-1:round(length(response_struct.session_date{animal_index})/2)+1, length(response_struct.session_date{animal_index})-2:length(response_struct.session_date{animal_index})}; %Session grouping indices for early, mid, and late sessions
    prob_matrix = response_struct.prob_response{animal_index, selected_stim_index}; %Get P(response)
    response_matrix = response_struct.normalized_response{animal_index, selected_stim_index}; %Get response
    response_threshold = prctile(response_matrix, response_pool_threshold, 'all'); %Get response threshold
    response_logicals = nansum(response_matrix > response_threshold, 2) > 0; %Logicals for ROIs with at least one session with response > response threshold
    response_logicals = response_logicals & roi_logicals{animal_index}; %Apply ROI logicals
    prob_matrix = prob_matrix(response_logicals, :); %Apply response threshold
    for group_index = 1:length(session_groups)
        summary_matrix(group_index, animal_index) = nanmedian(nanmean(prob_matrix(:, session_groups{group_index}), 2)); %Get median P(response) over ROIs in each session group
    end
    plot(summary_matrix(:, animal_index), 'o', 'Color', animal_coloring_scheme{animal_index}, 'MarkerSize', 6) %Plot marker for each animal
    hold on
    plot(summary_matrix(:, animal_index), '-', 'Color', animal_coloring_scheme{animal_index}, 'LineWidth', 2) %Plot line for each animal
end
errorbar(nanmean(summary_matrix, 2), nanstd(summary_matrix, 0, 2)/sqrt(length(animal_IDs)), 'o', 'Color', 'k', 'MarkerSize', 10, 'CapSize', 0) %Plot grand mean
xlim([0.5 3.5]); ylim([0 1])
xticks([1 2 3]); xticklabels({'Early', 'Middle', 'Late'});
xlabel('Sessions'); ylabel('P(response)')
title({['P(response) to ', num2str(stim_value_across_sessions), ' stimulus pulses']})
set(gcf, 'Position', [300, 300, 350, 200]) %Set figure window size to be consistent

%Stats
[~, t_test_p_em] = ttest(summary_matrix(1, :), summary_matrix(2, :)); %Paired t-test between early and middle
[~, t_test_p_ml] = ttest(summary_matrix(2, :), summary_matrix(3, :)); %Paired t-test between middle and late
corrected_p_prob = [t_test_p_em, t_test_p_ml] * 2; %Paired t-tests corrected for multiple comparisons (Bonferroni)

%% Neuron count (bottom row in panel)

%Plot number of ROIs with response > threshold (pool size)
figure
pool_size = NaN(length(animal_IDs), 1); %Matrix to hold pool size for each animal
for animal_index = 1:length(animal_IDs)
    response_matrix = response_struct.normalized_response{animal_index, selected_stim_index};
    response_threshold = prctile(response_matrix, response_pool_threshold, 'all'); %Get response threshold
    response_logicals = nansum(response_matrix > response_threshold, 2) > 0; %Logicals for ROIs with at least one session with response > response threshold
    response_logicals = response_logicals & roi_logicals{animal_index}; %Apply ROI logicals
    pool_size(animal_index) = sum(response_logicals); %Get pool size
    plot(pool_size(animal_index), 'o', 'Color', animal_coloring_scheme{animal_index}, 'MarkerSize', 6) %Plot marker for each animal    
    hold on
end
errorbar(nanmean(pool_size), nanstd(pool_size)/sqrt(length(animal_IDs)), 'o', 'Color', 'k', 'MarkerSize', 10, 'CapSize', 0) %Plot mean
ylim([0 3000]); ylabel('Neuron count')
title({['Number of neurons responding to ', num2str(stim_value_across_sessions), ' stimulus pulses']})
set(gcf, 'Position', [300, 300, 350, 200]) %Set figure window size to be consistent

