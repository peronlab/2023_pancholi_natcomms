%Script to plot Extended Data Figure 5: behavioral performance segregated by trial type
%
%To use this script, set path to animal_behavioral_data.mat file

%% Parameters
[all_dat settings] = get_all_data_pstim; 
behavioral_data_path = [settings.base_dir filesep];

%% Get data

%Load animal behavioral data
try
    load([behavioral_data_path, 'animal_behavioral_data.mat'])
catch
    error('Error loading animal behavioral data. Check file path and try again.')
end

%Compute hit, miss, false alarm, and correct reject rates for each animal
signal_detection_array = build_signal_detection_array(animal_behavioral_data.pre_lesion);

%Get session stimulus types (i.e. unique pulse counts presented in each session)
stim_type_array = build_stimulus_type_array(animal_behavioral_data.pre_lesion); %1 = [0 9], 2 = [0 1 9]/[1 9], 3 = [1 7 9]/[1 3 9]/[1 3 7 9], 4 = [1 3 5 7 9]

%% Plotting

animal_coloring_scheme = {[255, 31, 91]/255, [0, 154, 222]/255, [175, 88, 186]/255, [255, 198, 30]/255, [242, 133, 34]/255, [0, 255, 0]/255}; %Colors for different animals
num_animals = length(animal_behavioral_data.pre_lesion); %Set num animals to be the number of pre-lesion animals

%%% --- Plot hit and false alarm rates of all animals across sessions --- %%%
figure_position = [300, 300, 800, 400]; %[left bottom width height] in pixels. Value vector for figure 'Position' property to control figure window position and size.
figure
subplot(2, 1, 1) %Plot hit rate
for animal_index = 1:num_animals
    plot(signal_detection_array{animal_index}(:, 1)*100, '-', 'Color', animal_coloring_scheme{animal_index}, 'LineWidth', 4); %Plot line for hit rate for this animal
    hold on
end
xlim([1 max(cellfun(@(x) size(x, 1), signal_detection_array))])
ylim([40 100]); yticks([50 75 100]); ylabel({'Percent correct', 'on high pulse count trials'})
title('Percent correct by trial type')
subplot(2, 1, 2) %Plot false alarm rate
for animal_index = 1:num_animals
    plot(signal_detection_array{animal_index}(:, 3)*100, '-', 'Color', animal_coloring_scheme{animal_index}, 'LineWidth', 4); %Plot line for hit rate for this animal
    hold on
end
xlim([1 max(cellfun(@(x) size(x, 1), signal_detection_array))])
ylim([0 70]); yticks([20 40 60]); ylabel({'Percent incorrect', 'on low pulse count trials'})
xlabel('Session number'); 
set(gcf, 'Position', figure_position) %Set figure window size to be consistent

%%% --- Plot change in mean hit and false alarm rate upon introduction of 1 pulse stim --- %%%
figure_position = [300, 300, 800, 500]; %[left bottom width height] in pixels. Value vector for figure 'Position' property to control figure window position and size.
figure
num_sess_to_average = 1; %Number of sessions before/after lesion to average together. Uses N sessions before lesion and last N sessions after lesion.
mean_hit_pre = NaN(num_animals, 1); mean_hit_post = NaN(num_animals, 1); %Preallocate
mean_FA_pre = NaN(num_animals, 1); mean_FA_post = NaN(num_animals, 1);
all_hit_pre = NaN(num_animals*num_sess_to_average, 1); all_hit_post = NaN(num_animals*num_sess_to_average, 1);
all_FA_pre = NaN(num_animals*num_sess_to_average, 1); all_FA_post = NaN(num_animals*num_sess_to_average, 1);
for animal_index = 1:num_animals
    %Identify indices of task progresspion
    first_progress = find(stim_type_array{animal_index} == 2, 1); %Find index of switch to 9-1 task
    pre_indices = first_progress-num_sess_to_average:first_progress-1; %Pre-switch indices
    post_indices = first_progress:first_progress+num_sess_to_average-1; %Post-switch indices
    
    %Get mean of pre and post sessions
    mean_hit_pre(animal_index) = mean(signal_detection_array{animal_index}(pre_indices, 1))*100;
    mean_hit_post(animal_index) = mean(signal_detection_array{animal_index}(post_indices, 1))*100;
    mean_FA_pre(animal_index) = mean(signal_detection_array{animal_index}(pre_indices, 3))*100;
    mean_FA_post(animal_index) = mean(signal_detection_array{animal_index}(post_indices, 3))*100;
    
    %Get all pre and post values
    array_indices = animal_index*num_sess_to_average-(num_sess_to_average-1):animal_index*num_sess_to_average; %Get array indices
    all_hit_pre(array_indices) = signal_detection_array{animal_index}(pre_indices, 1)*100;
    all_hit_post(array_indices) = signal_detection_array{animal_index}(post_indices, 1)*100;
    all_FA_pre(array_indices) = signal_detection_array{animal_index}(pre_indices, 3)*100;
    all_FA_post(array_indices) = signal_detection_array{animal_index}(post_indices, 3)*100;
    
    %Plot
    hit_plot = subplot(1, 2, 1); %Plot hits
    plot(1, mean_hit_pre(animal_index), '.', 'Color', animal_coloring_scheme{animal_index}, 'MarkerSize', 20) %Plot pre-lesion average
    hold on
    plot(2, mean_hit_post(animal_index), '.', 'Color', animal_coloring_scheme{animal_index}, 'MarkerSize', 20) %Plot post-lesion average
    plot([1 2], [mean_hit_pre(animal_index) mean_hit_post(animal_index)], '-', 'Color', animal_coloring_scheme{animal_index}, 'LineWidth', 2) %Plot line connecting means to show change  
    FA_plot = subplot(1, 2, 2); %Plot FAs
    plot(1, mean_FA_pre(animal_index), '.', 'Color', animal_coloring_scheme{animal_index}, 'MarkerSize', 20) %Plot pre-lesion average
    hold on
    plot(2, mean_FA_post(animal_index), '.', 'Color', animal_coloring_scheme{animal_index}, 'MarkerSize', 20) %Plot post-lesion average
    plot([1 2], [mean_FA_pre(animal_index) mean_FA_post(animal_index)], '-', 'Color', animal_coloring_scheme{animal_index}, 'LineWidth', 2) %Plot line connecting means to show change  
end
errorbar(hit_plot, [1 2], [mean(all_hit_pre), mean(all_hit_post)], [std(all_hit_pre)/sqrt((num_animals-1)), std(all_hit_post)/sqrt((num_animals-1))], '.', 'Color', 'k', 'MarkerSize', 25) %Plot hit mean +/- SE
errorbar(FA_plot, [1 2], [mean(all_FA_pre), mean(all_FA_post)], [std(all_FA_pre)/sqrt((num_animals-1)), std(all_FA_post)/sqrt((num_animals-1))], '.', 'Color', 'k', 'MarkerSize', 25) %Plot FA mean +/- SE
plot(hit_plot, [1 2], [mean(all_hit_pre), mean(all_hit_post)], '-', 'Color', 'k', 'MarkerSize', 20) %Plot lines
plot(FA_plot, [1 2], [mean(all_FA_pre), mean(all_FA_post)], '-', 'Color', 'k', 'MarkerSize', 20)
xlim(hit_plot, [0 3]); xlim(FA_plot, [0 3]); ylim(hit_plot, [65 90]); ylim(FA_plot, [0 30])
xticks(hit_plot, [1 2]); xticklabels(hit_plot, {'Before', 'After'}); xticks(FA_plot, [1 2]); xticklabels(FA_plot, {'Before', 'After'})
yticks(hit_plot, [70 80 90]); yticks(FA_plot, [0 10 20 30])
ylabel(hit_plot, {'Percent correct', 'on high pulse count trials'})
ylabel(FA_plot, {'Percent incorrect', 'on low pulse count trials'})
suptitle('Change in performance before/after addition of 1 pulse stimulus')
set(gcf, 'Position', figure_position) %Set figure window size to be consistent
[~, p_hit,~,~] = ttest(mean_hit_pre, mean_hit_post); %Perform paired-sample t-test for hit rate
[~, p_FA,~,~] = ttest(mean_FA_pre, mean_FA_post); %Perform paired-sample t-test for FA rate

%% Helper functions

%Helper function to build signal detection arrays
function signal_detection_array = build_signal_detection_array(behav_data_structs)
    
    %Initialize
    session_mask_trials = 50; %Number of trials to remove at beginning and end of session when computing performance

    signal_detection_array = cell(1, length(behav_data_structs)); %Cell array to hold signal detection matrices for each animal
    for animal_index = 1:length(behav_data_structs) %Loop over animals
        
        %Preallocate
        signal_detection_matrix = NaN(length(behav_data_structs{animal_index}.session_date), 4); %Preallocate matrix (col1: hit rate, col2: miss rate, col3: FA rate, col4: CR rate)
        
        %Get logicals
        special_mode_logicals = (cellfun(@(x) ~strcmp(x, 'Normal'), {behav_data_structs{animal_index}.trial_settings(:).specialMode})' & ...
            cellfun(@(x) ~strcmp(x, 'Brutal'), {behav_data_structs{animal_index}.trial_settings(:).specialMode})') | ...
            cellfun(@(x) strcmp(x, 'Assist'), {behav_data_structs{animal_index}.trial_settings(:).lickportPrepositionState})'; %Get "special trials" (i.e. licking, timing, alternate) to be discounted from calculations
        correct_trials = [behav_data_structs{animal_index}.trial_outcomes(:).isCorrect]'; %Get correct trial logicals
        incorrect_trials = [behav_data_structs{animal_index}.trial_outcomes(:).isError]'; %Get incorrect trial logicals
        ignore_trials = [behav_data_structs{animal_index}.trial_outcomes(:).isIgnore]'; %Ignore trial logicals      
        right_trials = cellfun(@(x) strcmp(x, 'right'), {behav_data_structs{animal_index}.trial_settings(:).trialType})'; %Get right trials for this session
        
        %Loop over sessions        
        for session_index = 1:length(behav_data_structs{animal_index}.session_date) 
            
            %Get logicals for this session
            this_session_logicals = behav_data_structs{animal_index}.session_index == session_index; %Get logicals for trials in this session        
            last_non_ignore_index = find(this_session_logicals & ~ignore_trials, 1, 'last'); %Last index in logicals for this session that was not an ignore (to exclude spontaneous and LED stability tests)
            if sum(this_session_logicals) < (2*session_mask_trials + sum(this_session_logicals) - last_non_ignore_index) %If not enough trials, continue
                continue
            end
            session_mask = [zeros(find(this_session_logicals, 1)+session_mask_trials-1, 1); ...
                ones(sum(this_session_logicals) - 2*session_mask_trials - (find(this_session_logicals, 1, 'last') - last_non_ignore_index), 1); ...
                zeros(length(this_session_logicals)-find(this_session_logicals, 1, 'last')+(find(this_session_logicals, 1, 'last')-last_non_ignore_index)+session_mask_trials, 1)]; %Build session mask logicals that removes trials at beginning and end of session
            this_session_logicals = this_session_logicals & session_mask & ~special_mode_logicals; %Apply mask and remove special trials

            %Compute hit, miss, FA, and CR rates
            hits = correct_trials & this_session_logicals & right_trials; %Get hits for this session (hit = correct right trial)
            misses = incorrect_trials & this_session_logicals & right_trials; %Get misses for this session (miss = incorrect right trial)
            false_alarms = incorrect_trials & this_session_logicals & ~right_trials; %Get false alarms for this session (false alarm = incorrect left trial)
            correct_rejects = correct_trials & this_session_logicals & ~right_trials; %Get correct rejects for this session (correct rejects = correct left trial)
            hit_rate = sum(hits) / sum(hits | misses); %Hit rate
            miss_rate = 1 - hit_rate; %Miss rate
            false_alarm_rate = sum(false_alarms) / sum(false_alarms | correct_rejects); %False alarm rate
            correct_reject_rate = 1 - false_alarm_rate; %Correct reject rate
            
            %Add to matrix
            signal_detection_matrix(session_index, :) = [hit_rate, miss_rate, false_alarm_rate, correct_reject_rate];

        end
        signal_detection_array{animal_index} = signal_detection_matrix; %Add signal detection matrix to cell array
    end
end

%Helper function to get session stimulus sets (i.e. number of pulses presented)
function stim_type_array = build_stimulus_type_array(behav_data_structs)

    stim_type_array = cell(length(behav_data_structs), 1); %Cell array to hold stimulus type in each session for each animal
    for animal_index = 1:length(behav_data_structs) %Loop over animals
        
        %Initialize
        stimulus_type_vector = zeros(length(behav_data_structs{animal_index}.session_date), 1); %Preallocate vector to hold stimulus set for each session
        special_mode_logicals = (cellfun(@(x) ~strcmp(x, 'Normal'), {behav_data_structs{animal_index}.trial_settings(:).specialMode})' & ...
            cellfun(@(x) ~strcmp(x, 'Brutal'), {behav_data_structs{animal_index}.trial_settings(:).specialMode})') | ...
            cellfun(@(x) strcmp(x, 'Assist'), {behav_data_structs{animal_index}.trial_settings(:).lickportPrepositionState})'; %Get "special trials" (i.e. licking, timing, alternate) to be discounted from calculations
        
        %Get logicals
        pulse_counts = [behav_data_structs{animal_index}.trial_settings(:).optogenStimCount]'; %Get pulse counts for each trial
        pulse_amplitudes = cellfun(@mean, {behav_data_structs{animal_index}.trial_settings(:).optogenStimAmplitude})'; %Get pulse amplitudes for each trial
        pulse_counts = pulse_counts .* pulse_amplitudes; %Change 0 amplitude trials to be 0 pulse count
        
        %Loop over sessions
        for session_index = 1:length(behav_data_structs{animal_index}.session_date) 
            valid_trial_logicals = behav_data_structs{animal_index}.session_index == session_index; %Get logicals for trials in this session
            valid_trial_logicals = valid_trial_logicals & ~special_mode_logicals; %Remove special trials
            unique_pulse_counts = unique(pulse_counts(valid_trial_logicals)); %Get unique pulse counts in this session
            if length(unique_pulse_counts) == 2 && all(unique_pulse_counts == [0;9]) %Set stimulus type values. (1 = [0 9], 2 = [0 1 9]/[1 9], 3 = [1 7 9]/[1 3 9]/[1 3 7 9], 4 = [1 3 5 7 9])
                stimulus_type_vector(session_index) = 1;
            elseif (length(unique_pulse_counts) == 2 && all(unique_pulse_counts == [1;9])) || (length(unique_pulse_counts) == 3 && all(unique_pulse_counts == [0;1;9]))
                stimulus_type_vector(session_index) = 2;
            elseif (length(unique_pulse_counts) == 3 && all(unique_pulse_counts == [1;7;9])) || (length(unique_pulse_counts) == 3 && all(unique_pulse_counts == [1;3;9])) || ...
                (length(unique_pulse_counts) == 4 && all(unique_pulse_counts == [1;3;7;9]))
                stimulus_type_vector(session_index) = 3;
            elseif length(unique_pulse_counts) == 5 && all(unique_pulse_counts == [1;3;5;7;9])
                stimulus_type_vector(session_index) = 4;
            else
                warning(['Session contains a pulse count set that does not match intended pulse counts. Check animal index: ', num2str(animal_index), ...
                    ', session index: ', num2str(session_index)]); 
            end
        end
        stim_type_array{animal_index} = stimulus_type_vector; %Add stimulus vector to cell array
    end
end 
