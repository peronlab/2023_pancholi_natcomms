%
% Fig 3 : choice is stimulus based
%       - overall response for a single mouse, ops+/ops-
%       - overall response across mice, ops+/ops-                                     
%
function figure_3_choice (force_redo)
    if (nargin < 1) ; force_redo = 0 ; end

    % --- load relevant data
    [all_dat settings] = get_all_data_pstim; 
    beh_dat = get_pstim_beh_data(); 

    ex_ani = 'an293010';
    ex_sess = '01-Oct'; % si=23

    ex_ai = find(strcmp(settings.anims, ex_ani));
    ex_si = find(strcmp(all_dat.quick_dat{ex_ai}.date_str,ex_sess));
    opwd = pwd;
    load([settings.base_dir filesep 'pertrial_summary_example_mouse.mat'])
    
    % --- setup plots
    fh = figure ('Position',[0 0 1200 800]);
    example_evoked_pop_single_mouse_ax(1) = subplot('Position', [0.1 0.55 .233 .35]); 
    example_evoked_pop_single_mouse_ax(2) = subplot('Position', [0.1 0.1 .233 .35]); 
    example_evoked_cross_mouse_avg_ax(1) = subplot('Position', [0.45 0.55 .233 .35]); 
    example_evoked_cross_mouse_avg_ax(2) = subplot('Position', [0.45 0.1 .233 .35]); 
    example_evoked_cross_mouse_avg_stat_ax(1) = subplot('Position', [0.8 0.55 .1 .35]); 
    example_evoked_cross_mouse_avg_stat_ax(2) = subplot('Position', [0.8 0.1 .1 .35]); 

    ax = [example_evoked_pop_single_mouse_ax example_evoked_cross_mouse_avg_ax example_evoked_cross_mouse_avg_stat_ax];

    % --- plot - core code 

    %% population evoked, single session
    pstimi = find(all_dat.quick_dat{ex_ai}.probability_response{end}(ex_si,:) > .25);
    
    plot_population_evoked(example_evoked_pop_single_mouse_ax(1), pt_dat, ex_si, intersect(pstimi,all_dat.quick_dat{ex_ai}.types.red), 'ops+, photoresponsive');
    plot_population_evoked(example_evoked_pop_single_mouse_ax(2), pt_dat, ex_si, intersect(pstimi,all_dat.quick_dat{ex_ai}.types.green), 'ops-, photoresponsive');
 
    %% all mice evoked (avg across *all* terminal behavior sessions) 
    choice_dat = get_choice_dat(force_redo, settings, all_dat);

    mu_dat_r = nan*zeros(6,5) ; 
    mu_dat_l = nan*zeros(6,5) ;     
    for a=1:6 ; 
        for p=1:5 ; 
            mu_dat_r(a,p) = nanmean(squeeze(choice_dat(a).r_median_response(p,:,1))); % 1 = ops+ 
            mu_dat_l(a,p) = nanmean(squeeze(choice_dat(a).l_median_response(p,:,1)));
        end ; 
    end
    Mr = max(mu_dat_r');
    Ml = max(mu_dat_l');
    M = max(Mr,Ml);
    hold(example_evoked_cross_mouse_avg_ax(1),'on');
    plot(example_evoked_cross_mouse_avg_ax(1),mu_dat_l'./M,'b')
    plot(example_evoked_cross_mouse_avg_ax(1),mu_dat_l'./M,'bo', 'MarkerSize', 7)
    mu = nanmean((mu_dat_l'./M)');
    sem = nanstd((mu_dat_l'./M)')/sqrt(6);
    plot(example_evoked_cross_mouse_avg_ax(1), mu,'o','Color', 'None','MarkerFaceColor', [1 0 0], 'MarkerSize', 12);
    for p=1:5
        plot(example_evoked_cross_mouse_avg_ax(1), p*[1 1], mu(p)+[1 -1]*sem(p), '-','Color', [0 0 0], 'LineWidth', 2);
    end
    plot(example_evoked_cross_mouse_avg_ax(1),mu_dat_r'./M,'r')    
    plot(example_evoked_cross_mouse_avg_ax(1),mu_dat_r'./M,'ro')    
    mu = nanmean((mu_dat_r'./M)');
    sem = nanstd((mu_dat_r'./M)')/sqrt(6);
    plot(example_evoked_cross_mouse_avg_ax(1), mu,'o','Color', 'None','MarkerFaceColor', [0 0 1], 'MarkerSize', 12);
    for p=1:5
        plot(example_evoked_cross_mouse_avg_ax(1), p*[1 1], mu(p)+[1 -1]*sem(p), '-','Color', [0 0 0], 'LineWidth', 2);
    end
    title(example_evoked_cross_mouse_avg_ax(1), 'Ops+');

    plot(example_evoked_cross_mouse_avg_stat_ax(1), 1+zeros(1,6), nanmean((mu_dat_l'./M)), 'o', 'Color', 'None' ,'MarkerFaceColor', [0 0 1], 'MarkerSize',7);
    hold(example_evoked_cross_mouse_avg_stat_ax(1), 'on');
    plot(example_evoked_cross_mouse_avg_stat_ax(1), 1, nanmean(nanmean((mu_dat_l'./M))), 'o', 'Color', 'None' ,'MarkerFaceColor', [0 0 1], 'MarkerSize',12);
    plot(example_evoked_cross_mouse_avg_stat_ax(1), 2+zeros(1,6), nanmean((mu_dat_r'./M)), 'o', 'Color', 'None' ,'MarkerFaceColor', [1 0 0], 'MarkerSize',7);
    plot(example_evoked_cross_mouse_avg_stat_ax(1), 2, nanmean(nanmean((mu_dat_r'./M))), 'o', 'Color', 'None' ,'MarkerFaceColor', [1 0 0], 'MarkerSize',12);
    axis(example_evoked_cross_mouse_avg_stat_ax(1), [0 3 0 1]);

    mu_l = nanmean((mu_dat_l'./M));
    mu_r = nanmean((mu_dat_r'./M));
    [h p ] = ttest(mu_l, mu_r);

    title(example_evoked_cross_mouse_avg_stat_ax(1), sprintf('%s pval: %0.3f', 'ops+', p));

    mu_dat_r = nan*zeros(6,5) ; 
    mu_dat_l = nan*zeros(6,5) ;     
    for a=1:6 ; 
        for p=1:5 ; 
            mu_dat_r(a,p) = nanmean(squeeze(choice_dat(a).r_median_response(p,:,2))); % 2 = ops- 
            mu_dat_l(a,p) = nanmean(squeeze(choice_dat(a).l_median_response(p,:,2)));
        end ; 
    end
    Mr = max(mu_dat_r');
    Ml = max(mu_dat_l');
    M = max(Mr,Ml);
    hold(example_evoked_cross_mouse_avg_ax(2),'on');
    plot(example_evoked_cross_mouse_avg_ax(2),mu_dat_l'./M,'bo');
    plot(example_evoked_cross_mouse_avg_ax(2),mu_dat_l'./M,'b');
        mu = nanmean((mu_dat_l'./M)');
    sem = nanstd((mu_dat_l'./M)')/sqrt(6);
    plot(example_evoked_cross_mouse_avg_ax(2), mu,'o','Color', 'None','MarkerFaceColor', [1 0 0], 'MarkerSize', 12);
    for p=1:5
        plot(example_evoked_cross_mouse_avg_ax(2), p*[1 1], mu(p)+[1 -1]*sem(p), '-','Color', [0 0 0], 'LineWidth', 2);
    end
    plot(example_evoked_cross_mouse_avg_ax(2),mu_dat_r'./M,'r');    
    plot(example_evoked_cross_mouse_avg_ax(2),mu_dat_r'./M,'ro');    
    mu = nanmean((mu_dat_r'./M)');
    sem = nanstd((mu_dat_r'./M)')/sqrt(6);
    plot(example_evoked_cross_mouse_avg_ax(2), mu,'o','Color', 'None','MarkerFaceColor', [0 0 1], 'MarkerSize', 12);
    for p=1:5
        plot(example_evoked_cross_mouse_avg_ax(2), p*[1 1], mu(p)+[1 -1]*sem(p), '-','Color', [0 0 0], 'LineWidth', 2);
    end
    title(example_evoked_cross_mouse_avg_ax(2), 'Ops-');

    plot(example_evoked_cross_mouse_avg_stat_ax(2), 1+zeros(1,6), nanmean((mu_dat_l'./M)), 'o', 'Color', 'None' ,'MarkerFaceColor', [0 0 1], 'MarkerSize',7);
    hold(example_evoked_cross_mouse_avg_stat_ax(2), 'on');
    plot(example_evoked_cross_mouse_avg_stat_ax(2), 1, nanmean(nanmean((mu_dat_l'./M))), 'o', 'Color', 'None' ,'MarkerFaceColor', [0 0 1], 'MarkerSize',12);
    plot(example_evoked_cross_mouse_avg_stat_ax(2), 2+zeros(1,6), nanmean((mu_dat_r'./M)), 'o', 'Color', 'None' ,'MarkerFaceColor', [1 0 0], 'MarkerSize',7);
    plot(example_evoked_cross_mouse_avg_stat_ax(2), 2, nanmean(nanmean((mu_dat_r'./M))), 'o', 'Color', 'None' ,'MarkerFaceColor', [1 0 0], 'MarkerSize',12);
    axis(example_evoked_cross_mouse_avg_stat_ax(2), [0 3 0 1]);

    mu_l = nanmean((mu_dat_l'./M));
    mu_r = nanmean((mu_dat_r'./M));
    [h p ] = ttest(mu_l, mu_r);

    title(example_evoked_cross_mouse_avg_stat_ax(2), sprintf('%s pval: %0.3f', 'ops-', p));

    % --- finalize axes
    for a=1:length(ax)
        set(ax(a), 'TickDir', 'out', 'FontSize', 15);
    end

function choice_dat = get_choice_dat(force_redo, settings, all_dat)
    choice_dat_fname = [settings.base_dir filesep 'choice_dat.mat'];
    pulse_counts = 1:2:9;
    resp_right_trial_type = [2 3];
    resp_left_trial_type = [1 4];

    global choice_dat;

    if (isempty(choice_dat))
        ld = load(choice_dat_fname);
        choice_dat = ld.choice_dat;
    end

function plot_population_evoked(ax, pt_dat, ex_si, c_vali, tstr)
    %  trial_type_str: {'HitL'  'HitR'  'ErrL'  'ErrR'  'IgnoreL'  'IgnoreR'} - indices below correspond to this
    resp_right_trial_type = [2 3];
    resp_left_trial_type = [1 4];

    ddff = pt_dat.session(ex_si).dff_stim_epoch - pt_dat.session(ex_si).dff_pre_stim_epoch;

    hold(ax,'on');
    pulse_counts = 1:2:9;
    for p=1:length(pulse_counts)
        rpi = find(pt_dat.session(ex_si).num_pulses == pulse_counts(p) & ismember(pt_dat.session(ex_si).trial_type,resp_right_trial_type));
        lpi = find(pt_dat.session(ex_si).num_pulses == pulse_counts(p) & ismember(pt_dat.session(ex_si).trial_type,resp_left_trial_type));
        r_qs = quantile(nanmean(ddff(c_vali,rpi)), [.01 .25 .5 .75 .99]);
        plot(ax,pulse_counts(p)*[1 1]-0.25, r_qs([1 5]), 'k-', 'LineWidth',1);
        rectangle('Parent', ax, 'Position', [pulse_counts(p)-0.25-0.125 r_qs(2) .25 r_qs(4)-r_qs(2)], 'EdgeColor','None','FaceColor', [1 0 0]);
        plot(ax, pulse_counts(p)-0.25+[-0.125 0.125], r_qs(3)*[1 1],'k-', 'LineWidth',2);

        l_qs = quantile(nanmean(ddff(c_vali,lpi)), [.01 .25 .5 .75 .99]);
        plot(ax,pulse_counts(p)*[1 1]+0.25, l_qs([1 5]), 'k-', 'LineWidth',1);
        rectangle('Parent', ax, 'Position', [pulse_counts(p)+0.25-0.125 l_qs(2) .25 l_qs(4)-l_qs(2)], 'EdgeColor','None','FaceColor', [0 0 1]);
        plot(ax, pulse_counts(p)+0.25+[-0.125 0.125], l_qs(3)*[1 1],'k-', 'LineWidth',2);

        [hv pv] = ttest2(nanmean(ddff(c_vali,lpi)), nanmean(ddff(c_vali,rpi)));
        tstr = sprintf('%s %0.3f', tstr, pv);
    end
    aa=axis(ax);
    axis(ax, [0 10 aa(3) aa(4)]);
    axis(ax, [0 10 -0.2 1.2]);
    xlabel(ax,'Pulse count');
    ylabel(ax,'Mean evoked \DeltaF/F across neurons');
    title(ax, tstr);
    set(ax,'XTick',pulse_counts);
    
