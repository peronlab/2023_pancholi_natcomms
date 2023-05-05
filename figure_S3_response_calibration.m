%Script to plot Extended Data Figure 3: neural responsiveness varies as a function of stimulation parameters
%
%To use this script, first set path to response_calibration_data.mat file
%
%Next, define which calibration type you want to plot (set variable calibration type). The options are as follows:
%   - 1: Interpulse interval (panel b)
%   - 2: Pulse count vs. Pulse amplitude (panels a and d)
%   - 3: Pulse count vs. Pulse duration (panel c)
%
%Next, define whether you want to plot opsin-expressing ROIs or opsin non-expressing ROIs (set variable use_opsin_expressing_rois)
%   - 0: Use opsin non-expressing ROIs
%   - 1: Use opsin-expressing ROIs

%% Parameters
[all_dat settings] = get_all_data_pstim; 
calibration_data_path = [settings.base_dir filesep];
calibration_type = 2; %Define calibration type to plot. (1) = interpulse interval; (2) = Pulse count vs. Pulse amplitude; (3) Pulse count vs. Pulse duration
use_opsin_expressing_rois = 1; %Plot responses for opsin-expressing rois (1) or opsin non-expressing rois (0)

%% Load data

%Load animal behavioral data
try
    load([calibration_data_path, 'response_calibration_data.mat'])
catch
    error('Error loading response calibration data. Check file path and try again.')
end

%Define response and trial group structs
switch calibration_type
    case 1 %Interpulse interval
        trial_response_struct = response_calibration_data.interval_responses;
        valid_trial_groups_struct = response_calibration_data.interval_stim;
    case 2 %Count vs. amplitude
        trial_response_struct = response_calibration_data.count_vs_amplitude_responses;
        valid_trial_groups_struct = response_calibration_data.count_vs_amplitude_stim;    
    case 3 %Count vs. duration
        trial_response_struct = response_calibration_data.count_vs_duration_responses;
        valid_trial_groups_struct = response_calibration_data.count_vs_duration_stim;        
end

%% Get ROI responses

%Get ROIs to use
sess_roi_IDs = trial_response_struct.roi_ID; %ROIs in this session
if use_opsin_expressing_rois == 1
    roi_logicals = ismember(sess_roi_IDs, response_calibration_data.opsin_expressing_rois); %Opsin-expressing ROI logicals
elseif use_opsin_expressing_rois == 0
    roi_logicals = ismember(sess_roi_IDs, response_calibration_data.opsin_non_expressing_rois); %Opsin non-expressing ROI logicals 
else
    error('Invalid value for use_opsin_expressing_rois. Use either 1 for opsin-expressing or 0 for opsin non-expressing.')
end

%Get stimulus values
full_stim_value_matrix = valid_trial_groups_struct.parameter_value_matrix; %Col: Pulse count | duration | amplitude | interval
switch calibration_type
    case 1
        stim_values = full_stim_value_matrix(:, 4); %Get pulse interval
    case 2
        stim_values = full_stim_value_matrix(:, [1,3]); %Get pulse count and amplitude
    case 3
        stim_values = full_stim_value_matrix(:, [1,2]); %Get pulse count and duration
end

%Get stim combinations
stim_combinations = horzcat(valid_trial_groups_struct.unique_stim_combinations{:})';
total_stim_combinations = size(stim_combinations, 1);

%Preallocate
response_struct.frac_resp = cell(1, total_stim_combinations); %Fraction of ROIs responding to each stim combination
response_struct.prob_resp = NaN(length(sess_roi_IDs), total_stim_combinations); %Probability of response of each ROI for each stim combination
response_struct.norm_resp = cell(1, total_stim_combinations); %z-scored DFF of each ROI for each stim combination

%Loop over stim combinations
strongest_stim_combination = NaN; %Preallocate variable to hold index for stim combination with the strongest median normalized response
max_median_response = -Inf;
for stim_index = 1:total_stim_combinations
    %Get trials with this stim combination
    trial_logicals = ismember(stim_values, stim_combinations(stim_index, :), 'rows');
    
    %Get responses to this stim combination
    response_struct.frac_resp{stim_index} = nansum(trial_response_struct.logical_matrix(trial_logicals, roi_logicals), 2) ./ sum(~isnan(trial_response_struct.logical_matrix(trial_logicals, roi_logicals)), 2); %Fraction of ROIs responding
    response_struct.prob_resp(:, stim_index) =  nansum(trial_response_struct.logical_matrix(trial_logicals, :)) ./ sum(~isnan(trial_response_struct.logical_matrix(trial_logicals, :))); %Probability of response
    response_struct.norm_resp{stim_index} = trial_response_struct.normalized_response_matrix(trial_logicals, :); %z-scored DFF 
    
    %Find strongest stim combination
    median_response = nanmedian(nanmean(response_struct.norm_resp{stim_index}));
    if median_response > max_median_response
        strongest_stim_combination = stim_index;
        max_median_response = median_response;
    end
end

%% Plot

%Parameters
pulse_count_colors = {[0 0 0], [0 0 255]/255, [102 153 204]/255, [179 179 179]/255, [255 138 127]/255, [255 0 0]/255}; %Colors for pulse counts [0 1 3 5 7 9]
lower_percentile = 25; %Lower percentile when plotting percentiles
upper_percentile = 75; %Upper percentile when plotting percentiles
if use_opsin_expressing_rois; roi_type = 'Opsin+'; else; roi_type = 'Opsin-'; end %Set string for ROI color
figure_position = [100, 300, 1500, 500]; %[left bottom width height] in pixels. Value vector for figure 'Position' property to control figure window position and size.

%%
%%% --- Plot fraction of ROIs responding to each stim combination --- %%%

%Line plot
if calibration_type == 1 %Pulse interval
    
    %Get median and percentiles for each stim combination
    median_combination = cellfun(@median, response_struct.frac_resp); %Median
    lower_percentile_combination = cellfun(@(x) prctile(x, lower_percentile), response_struct.frac_resp); %Lower percentile
    upper_percentile_combination = cellfun(@(x) prctile(x, upper_percentile), response_struct.frac_resp); %Upper percentile
    
    %Plot
    figure
    errorbar(1:length(stim_combinations), median_combination, median_combination-lower_percentile_combination, upper_percentile_combination-median_combination, 'o', 'Color', 'k', 'MarkerSize', 10, 'CapSize', 0) %Plot median and IQR
    xlim([0.5, length(stim_combinations)+0.5]); ylim([0 0.3])
    xticks(1:length(stim_combinations)); xticklabels(cellstr(num2str(stim_combinations)));
    title({'Pulse interval', ['Fraction of ', roi_type, ' ROIs responding']})
    xlabel('Pulse interval (s)'); ylabel('Fraction of ROIs responding')
    
else %Pulse amplitude or duration
        
    %Plot fraction responding to max amplitude/duration
    figure
    unique_pulse_counts = unique(stim_combinations(:, 1)); %Get unique pulse counts   
    max_responses = NaN(3, length(unique_pulse_counts)); %Preallocate matrix to hold response median and percentiles for max amplitude/duration for each pulse count
    for count_index = 1:length(unique_pulse_counts) 
        %Get logicals for this pulse count among stim combinations
        count_logicals = stim_combinations(:, 1) == unique_pulse_counts(count_index); 

        %Get median and IQR for each stim combination
        median_combination = cellfun(@median, response_struct.frac_resp(count_logicals));
        lower_percentile_combination = cellfun(@(x) prctile(x, lower_percentile), response_struct.frac_resp(count_logicals));
        upper_percentile_combination = cellfun(@(x) prctile(x, upper_percentile), response_struct.frac_resp(count_logicals));
        max_responses(:, count_index) = [median_combination(end); lower_percentile_combination(end); upper_percentile_combination(end)]; %Add to max responses matrix        
    end    
    errorbar(1:length(unique_pulse_counts), max_responses(1, :), max_responses(1, :)-max_responses(2, :), max_responses(3, :)-max_responses(1, :), 'o', 'Color', 'k', 'MarkerSize', 10, 'CapSize', 0) %Plot marker
    xlim([0.5 length(unique_pulse_counts)+0.5]); ylim([0 0.3])
    xticks(1:length(unique_pulse_counts)); xticklabels(cellstr(num2str(unique_pulse_counts)));
    if calibration_type == 2
        title({'Pulse count vs. Pulse amplitude', ['Fraction of ', roi_type, ' ROIs responding to max amplitude']})
        xlabel('Pulse count'); ylabel('Fraction of ROIs responding')
    else
        title({'Pulse count vs. Pulse duration', ['Fraction of ', roi_type, ' ROIs responding to max duration']})
        xlabel('Pulse count'); ylabel('Fraction of ROIs responding')
    end
end

%Heatplot
if calibration_type > 1 %Do this only for amplitude or duration
    
    %Build matrix
    unique_pulse_counts = unique(stim_combinations(:, 1)); %Get unique pulse counts
    unique_pulse_var = unique(stim_combinations(:, 2)); %Get unique pulse amplitudes or durations
    response_matrix = NaN(length(unique_pulse_counts), length(unique_pulse_var)); %Preallocate
    for count_index = 1:length(unique_pulse_counts) %Loop over pulse counts
        %Get logicals for this pulse count among stim combinations
        count_logicals = stim_combinations(:, 1) == unique_pulse_counts(count_index); 

        %Get median of each stim combination
        median_combination = cellfun(@median, response_struct.frac_resp(count_logicals));

        %Add to matrix
        if unique_pulse_counts(count_index) == 0 %Zero pulse count
            response_matrix(:, 1) = repmat(median_combination, length(unique_pulse_counts), 1);
            response_matrix(1, :) = repmat(median_combination, length(unique_pulse_var), 1);
        else %Other combination
            response_matrix(count_index, 2:end) = median_combination;
        end 
    end
    
    %Plot
    figure
    custom_colormap = colormap_human; colormap(custom_colormap); %Use colormap 'human' defined in function colormap_human.m
    imagesc(response_matrix); %Heatmap of responses
    caxis([0 0.25])
    xticks(1:length(unique_pulse_var)); xticklabels(cellstr(num2str(unique_pulse_var)));
    yticks(1:length(unique_pulse_counts)); yticklabels(cellstr(num2str(unique_pulse_counts)));
    if calibration_type == 2
        title({'Pulse count vs. Pulse amplitude', ['Fraction of ', roi_type, ' ROIs responding']})
        xlabel('Pulse amplitude'); ylabel('Pulse count')
    else
        title({'Pulse count vs. Pulse duration', ['Fraction of ', roi_type, ' ROIs responding']})
        xlabel('Pulse duration'); ylabel('Pulse count')
    end
    colorbar
end
%%
%%% --- Plot prob response of ROIs responding to each stim combination --- %%%

%Get response logicals
response_matrix = nanmean(response_struct.norm_resp{strongest_stim_combination})'; %Get mean response to strongest stim combination for each ROI

%Line plot
if calibration_type == 1 %Pulse interval
    
    %Get median and IQR for each stim combination
    median_combination = nanmedian(response_struct.prob_resp(roi_logicals, :));
    lower_percentile_combination = prctile(response_struct.prob_resp(roi_logicals, :), lower_percentile);
    upper_percentile_combination = prctile(response_struct.prob_resp(roi_logicals, :), upper_percentile);
    
    %Plot
    figure
    errorbar(1:length(stim_combinations), median_combination, median_combination-lower_percentile_combination, upper_percentile_combination-median_combination, 'o', 'Color', 'k', 'MarkerSize', 10, 'CapSize', 0) %Plot median and IQR
    xlim([0.5, length(stim_combinations)+0.5]); ylim([0 0.4])
    xticks(1:length(stim_combinations)); xticklabels(cellstr(num2str(stim_combinations)));
    title({'Pulse interval', ['P(response) of ', roi_type, ' ROIs']})
    xlabel('Pulse interval (s)'); ylabel('P(response)')
    
else %Pulse amplitude or duration
    
    %Plot responses to max amplitude/duration
    figure
    unique_pulse_counts = unique(stim_combinations(:, 1)); %Get unique pulse counts
    max_responses = NaN(3, length(unique_pulse_counts)); %Preallocate matrix to hold response median and percentiles for max amplitude/duration for each pulse count
    for count_index = 1:length(unique_pulse_counts) 
        %Get logicals for this pulse count among stim combinations
        count_logicals = stim_combinations(:, 1) == unique_pulse_counts(count_index); 

        %Get median and IQR for each stim combination
        median_combination = nanmedian(response_struct.prob_resp(roi_logicals, count_logicals));
        lower_percentile_combination = prctile(response_struct.prob_resp(roi_logicals, count_logicals), lower_percentile);
        upper_percentile_combination = prctile(response_struct.prob_resp(roi_logicals, count_logicals), upper_percentile);
        max_responses(:, count_index) = [median_combination(end); lower_percentile_combination(end); upper_percentile_combination(end)]; %Add to max responses matrix
    end
    errorbar(1:length(unique_pulse_counts), max_responses(1, :), max_responses(1, :)-max_responses(2, :), max_responses(3, :)-max_responses(1, :), 'o', 'Color', 'k', 'MarkerSize', 10, 'CapSize', 0) %Plot marker
    xlim([0.5 length(unique_pulse_counts)+0.5]); ylim([0 0.4])
    xticks(1:length(unique_pulse_counts)); xticklabels(cellstr(num2str(unique_pulse_counts)));
    if calibration_type == 2
        title({'Pulse count vs. Pulse amplitude', ['P(response) of ', roi_type, ' ROIs to max amplitude']})
        xlabel('Pulse count'); ylabel('P(response)')
    else
        title({'Pulse count vs. Pulse duration', ['P(response) of ', roi_type, ' ROIs to max duration']})
        xlabel('Pulse count'); ylabel('P(response)')
    end
end

%Heatplot
if calibration_type > 1 %Do this only for amplitude or duration
    
    %Build matrix
    unique_pulse_counts = unique(stim_combinations(:, 1)); %Get unique pulse counts
    unique_pulse_var = unique(stim_combinations(:, 2)); %Get unique pulse amplitudes or durations
    response_matrix = NaN(length(unique_pulse_counts), length(unique_pulse_var)); %Preallocate
    for count_index = 1:length(unique_pulse_counts) %Loop over pulse counts
        %Get logicals for this pulse count among stim combinations
        count_logicals = stim_combinations(:, 1) == unique_pulse_counts(count_index); 

        %Get median of each stim combination
        median_combination = median(response_struct.prob_resp(roi_logicals, count_logicals));

        %Add to matrix
        if unique_pulse_counts(count_index) == 0 %Zero pulse count
            response_matrix(:, 1) = repmat(median_combination, length(unique_pulse_counts), 1);
            response_matrix(1, :) = repmat(median_combination, length(unique_pulse_var), 1);
        else %Other combination
            response_matrix(count_index, 2:end) = median_combination;
        end 
    end
    
    %Plot
    figure
    custom_colormap = colormap_human; colormap(custom_colormap); %Use colormap 'human' defined in function colormap_human.m
    imagesc(response_matrix); %Heatmap of responses
    caxis([0 0.2])
    xticks(1:length(unique_pulse_var)); xticklabels(cellstr(num2str(unique_pulse_var)));
    yticks(1:length(unique_pulse_counts)); yticklabels(cellstr(num2str(unique_pulse_counts)));
    if calibration_type == 2
        title({'Pulse count vs. Pulse amplitude', ['P(response) of ', roi_type, ' ROIs']})
        xlabel('Pulse amplitude'); ylabel('Pulse count')
    else
        title({'Pulse count vs. Pulse duration', ['P(response) of ', roi_type, ' ROIs']})
        xlabel('Pulse duration'); ylabel('Pulse count')
    end
    colorbar
end
%%
%%% --- Plot z-scored DFF response of ROIs responding to each stim combination --- %%%

%Get response array
response_matrix = nanmean(response_struct.norm_resp{strongest_stim_combination})'; %Get mean response to strongest stim combination
response_array = cellfun(@(x) x(:, roi_logicals), response_struct.norm_resp, 'UniformOutput', false); %Apply roi logicals

%Line plot
if calibration_type == 1 %Pulse interval
  
    %Get median and IQR for each stim combination
    median_combination = cellfun(@(x) nanmedian(nanmean(x)), response_array);
    lower_percentile_combination = cellfun(@(x) prctile(nanmean(x), lower_percentile), response_array);
    upper_percentile_combination = cellfun(@(x) prctile(nanmean(x), upper_percentile), response_array);
    
    %Plot
    figure
    errorbar(1:length(stim_combinations), median_combination, median_combination-lower_percentile_combination, upper_percentile_combination-median_combination, 'o', 'Color', 'k', 'MarkerSize', 10, 'CapSize', 0) %Plot median and IQR
    xlim([0.5, length(stim_combinations)+0.5]); ylim([0 2])
    xticks(1:length(stim_combinations)); xticklabels(cellstr(num2str(stim_combinations)));
    title({'Pulse interval', ['z-scored \DeltaF/F of ', roi_type, ' ROIs']})
    xlabel('Pulse interval (s)'); ylabel('z-scored \DeltaF/F')
 
else %Pulse amplitude or duration

    %Plot responses to max amplitude/duration
    figure
    unique_pulse_counts = unique(stim_combinations(:, 1)); %Get unique pulse counts
    max_responses = NaN(3, length(unique_pulse_counts)); %Preallocate matrix to hold response median and percentiles for max amplitude/duration for each pulse count
    for count_index = 1:length(unique_pulse_counts) 
        %Get logicals for this pulse count among stim combinations
        count_logicals = stim_combinations(:, 1) == unique_pulse_counts(count_index); 

        %Get median and IQR for each stim combination
        median_combination = cellfun(@(x) nanmedian(nanmean(x)), response_array(count_logicals));
        lower_percentile_combination = cellfun(@(x) prctile(nanmean(x), lower_percentile), response_array(count_logicals));
        upper_percentile_combination = cellfun(@(x) prctile(nanmean(x), upper_percentile), response_array(count_logicals));
        max_responses(:, count_index) = [median_combination(end); lower_percentile_combination(end); upper_percentile_combination(end)]; %Add to max responses matrix
        
    end
    errorbar(1:length(unique_pulse_counts), max_responses(1, :), max_responses(1, :)-max_responses(2, :), max_responses(3, :)-max_responses(1, :), 'o', 'Color', 'k', 'MarkerSize', 10, 'CapSize', 0) %Plot marker
    xlim([0.5 length(unique_pulse_counts)+0.5]); ylim([-0.5 2])
    xticks(1:length(unique_pulse_counts)); xticklabels(cellstr(num2str(unique_pulse_counts)));
    if calibration_type == 2
        title({'Pulse count vs. Pulse amplitude', ['z-scored \DeltaF/F of ', roi_type, ' ROIs to max amplitude']})
        xlabel('Pulse count'); ylabel('z-scored \DeltaF/F')
    else
        title({'Pulse count vs. Pulse duration', ['z-scored \DeltaF/F of ', roi_type, ' ROIs to max duration']})
        xlabel('Pulse count'); ylabel('z-scored \DeltaF/F')
    end
    
end

%Heatplot
if calibration_type > 1 %Do this only for amplitude or duration
    
    %Build matrix
    unique_pulse_counts = unique(stim_combinations(:, 1)); %Get unique pulse counts
    unique_pulse_var = unique(stim_combinations(:, 2)); %Get unique pulse amplitudes or durations
    response_matrix = NaN(length(unique_pulse_counts), length(unique_pulse_var)); %Preallocate
    for count_index = 1:length(unique_pulse_counts) %Loop over pulse counts
        %Get logicals for this pulse count among stim combinations
        count_logicals = stim_combinations(:, 1) == unique_pulse_counts(count_index); 

        %Get mean of each ROI's response to a specific stim combination and then median over all ROIs
        median_combination = cellfun(@(x) nanmedian(nanmean(x)), response_array(count_logicals));

        %Add to matrix
        if unique_pulse_counts(count_index) == 0 %Zero pulse count
            response_matrix(:, 1) = repmat(median_combination, length(unique_pulse_counts), 1);
            response_matrix(1, :) = repmat(median_combination, length(unique_pulse_var), 1);
        else %Other combination
            response_matrix(count_index, 2:end) = median_combination;
        end 
    end
    
    %Plot
    figure
    custom_colormap = colormap_human; colormap(custom_colormap); %Use colormap 'human' defined in function colormap_human.m
    imagesc(response_matrix); %Heatmap of responsesmax_median_response
    caxis([0 1])
    xticks(1:length(unique_pulse_var)); xticklabels(cellstr(num2str(unique_pulse_var)));
    yticks(1:length(unique_pulse_counts)); yticklabels(cellstr(num2str(unique_pulse_counts)));
    if calibration_type == 2
        title({'Pulse count vs. Pulse amplitude', ['z-scored \DeltaF/F of ', roi_type, ' ROIs']})
        xlabel('Pulse amplitude'); ylabel('Pulse count')
    else
        title({'Pulse count vs. Pulse duration', ['z-scored \DeltaF/F of ', roi_type, ' ROIs']})
        xlabel('Pulse duration'); ylabel('Pulse count')
    end
    colorbar
end

