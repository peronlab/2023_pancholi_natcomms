%Script to plot Figure 2f-h: photostimulation response increases with pulse count
%
%To use this script, first set path to session_response_data.mat file
%
%Next, define whether you want to plot opsin-expressing ROIs or opsin non-expressing ROIs (set variable use_opsin_expressing_rois)
%   - 0: Use opsin non-expressing ROIs
%   - 1: Use opsin-expressing ROIs

%% Parameters
[all_dat settings] = get_all_data_pstim; 
response_data_path = [settings.base_dir filesep];
use_opsin_expressing_rois = 0; %Plot responses for opsin-expressing rois (1) or opsin non-expressing rois (0)

%% Get data

%Load session response data
try
    load([response_data_path, 'session_response_data.mat'])
catch
    error('Error loading session response data. Check file path and try again.')
end

%Get ROIs to use
roi_logicals = cell(1, length(animal_IDs)); %Preallocate
for animal_index = 1:length(animal_IDs) %Loop over animals
    if use_opsin_expressing_rois == 1
        roi_logicals{animal_index} = ismember(response_struct.all_roi_IDs{animal_index}, response_struct.red_roi_IDs{animal_index});
    elseif use_opsin_expressing_rois == 0
        roi_logicals{animal_index} = ismember(response_struct.all_roi_IDs{animal_index}, response_struct.green_roi_IDs{animal_index});
    else
        error('Invalid value for use_opsin_expressing_rois. Use either 1 for opsin-expressing or 0 for opsin non-expressing.')
    end
end

%% Plot 

%Figure parameters
animal_coloring_scheme = {[255, 31, 91]/255, [0, 154, 222]/255, [175, 88, 186]/255, [255, 198, 30]/255, [242, 133, 34]/255, [0 255 0]/255}; %Colors for different animals
figure_position = [100, 300, 1500, 500]; %[left bottom width height] in pixels. Value vector for figure 'Position' property to control figure window position and size.

%% Fraction of ROIs (panel 2f)

%Boxplot
figure
example_animal_ID = 'an295578'; %Example animal ID
example_session = '2021_01_10'; %Example session date
animal_index = find(cellfun(@(x) strcmp(x, example_animal_ID), animal_IDs), 1); %Get index of example animal
session_index = find(cellfun(@(x) strcmp(x, example_session), response_struct.session_date{animal_index}), 1); %Get index of example session
frac_resp_vector = []; grouping_vector = []; %Preallocate
for stim_index = 1:length(possible_stim_values)
    stim_logicals = response_struct.stim_values{animal_index}(response_struct.stim_session_indices{animal_index} == session_index) == possible_stim_values(stim_index); %Get logicals for this stim value in this session
    logical_matrix = response_struct.response_logicals{animal_index}{session_index}(stim_logicals, roi_logicals{animal_index})'; %Get logical matrix
    frac_resp_trial = nansum(logical_matrix) ./ sum(~isnan(logical_matrix)); %Fraction of logical ROIs responsive to this stim value on this session
    frac_resp_vector = [frac_resp_vector; frac_resp_trial']; %Add to vector
    grouping_vector = [grouping_vector; repmat(stim_index, length(frac_resp_trial), 1)]; %Add stim index to grouping vector
end
boxplot(frac_resp_vector, grouping_vector, 'Widths', 0.3, 'Symbol', '') %Boxplot
ylim([0 0.5])
xticks((1:5)); xticklabels({'1', '3', '5', '7', '9'})
xlabel('Pulse count'); ylabel('Fraction of ROIs')
title({'Fraction of ROIs responding to each pulse count', [example_animal_ID, ', ', replace(example_session, '_', '-')]})

%Plot median fraction of ROIs responding to each pulse count across full pulse sessions for each animal
figure
frac_rois_resp = NaN(length(full_pulse_stim_values), length(animal_IDs)); %Preallocate
for animal_index = 1:length(animal_IDs)
    full_sess_indices = find(response_struct.stim_phase{animal_index} == 5); %Get indices for full pulse sessions
    for stim_index = 1:length(full_pulse_stim_values) %Loop over stim pulses in full pulse sessions
        frac_resp_vector = [];
        for session_index = 1:length(full_sess_indices)
            stim_logicals = response_struct.stim_values{animal_index}(response_struct.stim_session_indices{animal_index} == full_sess_indices(session_index)) == full_pulse_stim_values(stim_index); %Get logicals for this stim value in this session
            logical_matrix = response_struct.response_logicals{animal_index}{full_sess_indices(session_index)}(stim_logicals, roi_logicals{animal_index})'; %Get logical matrix
            frac_resp_trial = nansum(logical_matrix) ./ sum(~isnan(logical_matrix)); %Fraction of ROIs responsive 
            frac_resp_vector = [frac_resp_vector; frac_resp_trial']; %Add to vector
        end        
        frac_rois_resp(stim_index, animal_index) = nanmedian(frac_resp_vector); %Get median fraction of ROIs (of this type) responding to this stim value across full pulse sessions
    end
    plot(frac_rois_resp(:, animal_index), 'o', 'Color', animal_coloring_scheme{animal_index}, 'MarkerSize', 6) %Plot median frac resp
    hold on
    plot(frac_rois_resp(:, animal_index), '-', 'Color', animal_coloring_scheme{animal_index}, 'LineWidth', 1) %Plot median frac resp
end
errorbar(mean(frac_rois_resp, 2), std(frac_rois_resp, 0, 2)/sqrt(length(animal_IDs)), 'o', 'Color', 'k', 'MarkerSize', 10, 'CapSize', 0) %Plot grand mean +/- SEM across animals
xlim([0.5 5.5]); ylim([0 0.5])
xticks(1:length(full_pulse_stim_values)); xticklabels({'1', '3', '5', '7', '9'})
title({'Median fraction of ROIs responding to each pulse count', 'Full pulse count sessions'})
xlabel('Pulse count'); ylabel('Fraction of ROIs')

%% Z-scored DFF (panel 2g)

%Boxplot for example animal
figure
example_animal_ID = 'an295578'; %Example animal ID
example_session = '2021_01_10'; %Example session date
animal_index = find(cellfun(@(x) strcmp(x, example_animal_ID), animal_IDs), 1); %Get index of example animal
session_index = find(cellfun(@(x) strcmp(x, example_session), response_struct.session_date{animal_index}), 1); %Get index of example animal
plotting_matrix = NaN(sum(roi_logicals{animal_index}), length(full_pulse_stim_values)); %Preallocate
for stim_index = 2:length(possible_stim_values)
    plotting_matrix(:, stim_index-1) = response_struct.normalized_response{animal_index, stim_index}(roi_logicals{animal_index}, session_index);
end
boxplot(plotting_matrix, 'Widths', 0.3, 'Symbol', '')
ylim([-2 4])
xticks((1:5)); xticklabels({'1', '3', '5', '7', '9'})
xlabel('Pulse count'); ylabel('Normalized \DeltaF/F')
title({'Normalized \DeltaF/F of all ROIs to each pulse count', [example_animal_ID, ', ', replace(example_session, '_', '-')]})

%%% --- Plot median normalized response to each pulse count across full pulse sessions for each animal --- %%%
figure
median_norm_response = NaN(length(full_pulse_stim_values), length(animal_IDs)); %Preallocate
for animal_index = 1:length(animal_IDs)
    full_sess_logicals = response_struct.stim_phase{animal_index} == 5; %Get logicals for full pulse sessions
    for stim_index = 1:length(full_pulse_stim_values) %Loop over stim pulses in full pulse sessions    
        full_pulse_stim_index = find(possible_stim_values == full_pulse_stim_values(stim_index), 1); %Find index of full pulse session stim value in possible stim values       
        full_sess_norm_resp = response_struct.normalized_response{animal_index, full_pulse_stim_index}(roi_logicals{animal_index}, full_sess_logicals); %Get normalized response of ROIs for each full pulse session
        median_norm_response(stim_index, animal_index) = nanmedian(nanmean(full_sess_norm_resp, 2)); %Get mean of each ROIs response across full pulse sessions and then take median over all ROIs
    end
    %Plot
    plot(median_norm_response(:, animal_index), 'o', 'Color', animal_coloring_scheme{animal_index}, 'MarkerSize', 6) %Plot mean normalized response for this animal
    hold on
    plot(median_norm_response(:, animal_index), '-', 'Color', animal_coloring_scheme{animal_index}, 'LineWidth', 1) %Plot mean normalized response for this animal
end
errorbar(mean(median_norm_response, 2), std(median_norm_response, 0, 2)/sqrt(length(animal_IDs)), 'o', 'Color', 'k', 'MarkerSize', 10, 'CapSize', 0) %Plot grand mean across animals
xlim([0.5 5.5]); ylim([0 2])
xticks(1:length(full_pulse_stim_values)); xticklabels({'1', '3', '5', '7', '9'})
title({'Median normalized \DeltaF/F of ROIs responding to each pulse count', 'Full pulse count sessions'})
xlabel('Pulse count'); ylabel('Normalized \DeltaF/F')

%% Probability of response (panel 2h)

%Boxplot for example animal
figure
example_animal_ID = 'an295578'; %Example animal ID
example_session = '2021_01_10'; %Example session date
animal_index = find(cellfun(@(x) strcmp(x, example_animal_ID), animal_IDs), 1); %Get index of example animal
session_index = find(cellfun(@(x) strcmp(x, example_session), response_struct.session_date{animal_index}), 1); %Get index of example animal
plotting_matrix = NaN(sum(roi_logicals{animal_index}), length(full_pulse_stim_values)); %Preallocate
for stim_index = 2:length(possible_stim_values)
    plotting_matrix(:, stim_index-1) = response_struct.prob_response{animal_index, stim_index}(roi_logicals{animal_index}, session_index);
end
boxplot(plotting_matrix, 'Widths', 0.3, 'Symbol', '')
ylim([0 1])
xticks((1:5)); xticklabels({'1', '3', '5', '7', '9'})
xlabel('Pulse count'); ylabel('P(response)')
title({'P(response) of all ROIs to each pulse count', [example_animal_ID, ', ', replace(example_session, '_', '-')]})

%Plot median P(response) to each pulse count across full pulse sessions for each animal
figure
median_prob_response = NaN(length(full_pulse_stim_values), length(animal_IDs)); %Preallocate
for animal_index = 1:length(animal_IDs)
    full_sess_logicals = response_struct.stim_phase{animal_index} == 5; %Get logicals for full pulse sessions
    for stim_index = 1:length(full_pulse_stim_values) %Loop over stim pulses in full pulse sessions    
        full_pulse_stim_index = find(possible_stim_values == full_pulse_stim_values(stim_index), 1); %Find index of full pulse session stim value in possible stim values       
        full_sess_prob_resp = response_struct.prob_response{animal_index, full_pulse_stim_index}(roi_logicals{animal_index}, full_sess_logicals); %Get prob response of ROIs for each full pulse session
        median_prob_response(stim_index, animal_index) = nanmedian(nanmean(full_sess_prob_resp, 2)); %Get mean of each ROIs response across full pulse sessions and then take median over all ROIs
    end
    %Plot
    plot(median_prob_response(:, animal_index), 'o', 'Color', animal_coloring_scheme{animal_index}, 'MarkerSize', 6) %Plot median prob resp for this animal
    hold on
    plot(median_prob_response(:, animal_index), '-', 'Color', animal_coloring_scheme{animal_index}, 'LineWidth', 1) %Plot median prob resp for this animal
end
errorbar(mean(median_prob_response, 2), std(median_prob_response, 0, 2)/sqrt(length(animal_IDs)), 'o', 'Color', 'k', 'MarkerSize', 10, 'CapSize', 0) %Plot grand mean across animals
xlim([0.5 5.5]); ylim([0 0.5])
xticks(1:length(full_pulse_stim_values)); xticklabels({'1', '3', '5', '7', '9'})
title({'Median P(response) of ROIs responding to each pulse count', 'Full pulse count sessions'})
xlabel('Pulse count'); ylabel('P(response)')
