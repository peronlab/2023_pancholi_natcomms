%Script to plot Figure 1h: change in animal performance (percent correct) after lesion of opsin-expressing area

%% Parameters

%Data loading parameters
[all_dat settings] = get_all_data_pstim; 
behavioral_data_path = settings.base_dir;

%% Get data

%Load animal behavioral data
try
    load([behavioral_data_path filesep 'animal_behavioral_data.mat'])
catch
    error('Error loading animal behavioral data. Check file path and try again.')
end

%Compute percent correct for each session for each animal
percent_correct_array_pre_lesion = build_percent_correct_array(animal_behavioral_data.pre_lesion);
percent_correct_array_post_lesion = build_percent_correct_array(animal_behavioral_data.post_lesion);

%% Plotting

%%% --- Plot change in percent correct before and after lesion --- %%%

figure_position = [300, 300, 400, 500]; %[left bottom width height] in pixels. Value vector for figure 'Position' property to control figure window position and size.
figure
num_sess_to_average = 3; %Number of sessions before/after lesion to average together. Uses N sessions before lesion and last N sessions after lesion.
num_animals = length(animal_behavioral_data.post_lesion); %Set num animals to be the number of post-lesion animals (only first 5 animals were lesioned)
mean_pre_lesion = NaN(num_animals, 1); mean_post_lesion = NaN(num_animals, 1); %Preallocate
for animal_index = 1:num_animals
    pre_lesion_session_indices = (length(percent_correct_array_pre_lesion{animal_index})-num_sess_to_average+1):length(percent_correct_array_pre_lesion{animal_index}); %Pre-lesion session indices to average over
    post_lesion_session_indices = (length(percent_correct_array_post_lesion{animal_index})-num_sess_to_average+1):length(percent_correct_array_post_lesion{animal_index}); %Post-lesion session indices to average over
    mean_pre_lesion(animal_index) = mean(percent_correct_array_pre_lesion{animal_index}(pre_lesion_session_indices)); %Mean of 3 sessions pre-lesion
    mean_post_lesion(animal_index) = mean(percent_correct_array_post_lesion{animal_index}(post_lesion_session_indices)); %Mean of 3 session post-lesion
    plot(1, mean_pre_lesion(animal_index), '.', 'Color', 'k', 'MarkerSize', 20) %Plot pre-lesion average
    hold on
    plot(2, mean_post_lesion(animal_index), '.', 'Color', 'k', 'MarkerSize', 20) %Plot post-lesion average
    plot([1 2], [mean_pre_lesion(animal_index) mean_post_lesion(animal_index)], '-', 'Color', 'k', 'LineWidth', 2) %Plot line connecting means to show change  
end
all_pre_lesion = cellfun(@(x) x((end-num_sess_to_average+1):end), percent_correct_array_pre_lesion, 'UniformOutput', false); all_pre_lesion = vertcat(all_pre_lesion{:}); %Get pre-lesion percent correct for each animal
all_post_lesion = cellfun(@(x) x((end-num_sess_to_average+1):end), percent_correct_array_post_lesion, 'UniformOutput', false); all_post_lesion = vertcat(all_post_lesion{:}); %Get post-lesion percent correct for each animal
errorbar([1 2], [mean(all_pre_lesion), mean(all_post_lesion)], [std(all_pre_lesion)/sqrt(num_animals), std(all_post_lesion)/sqrt(num_animals)], '.', 'Color', 'r', 'MarkerSize', 25) %Plot mean +/- SE
plot([1 2], [mean(all_pre_lesion), mean(all_post_lesion)], '-', 'Color', 'r', 'MarkerSize', 20)
xlim([0 3]); ylim([0.5 1])
title('Change in percent correct pre-post lesion')
xticks([1 2]); xticklabels({'Pre-lesion', 'Post-lesion'})
ylabel('Percent correct')
set(gcf, 'Position', figure_position) %Set figure window size to be consistent
[h,p,ci,stats] = ttest(mean_pre_lesion, mean_post_lesion); %Perform paired-sample t-test

%% Helper functions

%Helper function to build percent correct array
function percent_correct_array = build_percent_correct_array(behav_data_structs)
    
    %Initialize
    session_mask_trials = 50; %Number of trials to remove at beginning and end of session when computing performance

    percent_correct_array = cell(1, length(behav_data_structs)); %Cell array to hold session percent correct vectors for each animal    
    for animal_index = 1:length(behav_data_structs) %Loop over animals

        %Preallocate
        percent_correct_vector = NaN(length(behav_data_structs{animal_index}.session_date), 1); %Preallocate vector to hold percent correct for each session       
        
        %Get logicals
        special_mode_logicals = (cellfun(@(x) ~strcmp(x, 'Normal'), {behav_data_structs{animal_index}.trial_settings(:).specialMode})' & ...
            cellfun(@(x) ~strcmp(x, 'Brutal'), {behav_data_structs{animal_index}.trial_settings(:).specialMode})') | ...
            cellfun(@(x) strcmp(x, 'Assist'), {behav_data_structs{animal_index}.trial_settings(:).lickportPrepositionState})'; %Get "special trials" (i.e. licking, timing, alternate) to be discounted from calculations
        correct_trials = [behav_data_structs{animal_index}.trial_outcomes(:).isCorrect]'; %Get correct trial logicals
        incorrect_trials = [behav_data_structs{animal_index}.trial_outcomes(:).isError]'; %Get incorrect trial logicals
        ignore_trials = [behav_data_structs{animal_index}.trial_outcomes(:).isIgnore]'; %Ignore trial logicals

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

            %Compute percent correct
            session_percent_correct = sum(correct_trials & this_session_logicals) / (sum(correct_trials & this_session_logicals) + sum(incorrect_trials & this_session_logicals)); %Compute percent correct (ignores ignore trials)
            percent_correct_vector(session_index) = session_percent_correct; %Add session percent correct to vector 
        end
        percent_correct_array{animal_index} = percent_correct_vector; %Add percent correct vector to cell array        
    end
end
