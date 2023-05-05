function figure_4gh_activity_and_behavior
    top_frac = 1; % take ALL neurons
    single_plot (top_frac, 1:3); % early
    single_plot (top_frac, 6:8); % late


function single_plot (top_frac, stat_di)
    % --- pull data
    [all_dat settings] = get_all_data_pstim;

    disp(sprintf('Top frac: %0.2f day %d to %d', top_frac, stat_di(1), stat_di(end)));
    
    % --- setup plots
    fh = figure ('Position',[0 0 1200 900], 'Name' , sprintf('Top frac: %0.2f day %d to %d', top_frac, stat_di(1), stat_di(end)));

    for c=1:4
        corr_ax(c) = subplot('Position', [0.05+0.25*(c-1) 0.75 .15 .2]);  % {presp , zsc dff} X {ops- , ops+}
        day_by_day_ax(c) = subplot('Position', [0.05+0.25*(c-1) 0.4 .15 .2]); % {presp , zsc dff} X {ops- , ops+}
        stats_ax(c) = subplot('Position', [0.05+0.25*(c-1) 0.05 .15 .2]);  % {presp , zsc dff} X {ops- , ops+}
    end

    ax = [day_by_day_ax stats_ax corr_ax];

    % --- do the dew
    learn_col = [1 .25 0];
    nlearn_col = [0 0 0 ];

    max_d_learni = 22;
    summ_dat_learni = nan*zeros(length(all_dat.settings.learni), max_d_learni, 2, 2); % anim , day , presp/zsc dff, ops+/ops-
    for a=1:length(all_dat.settings.learni)
        ai = all_dat.settings.learni(a);
        qd = all_dat.quick_dat{ai};

        zsc_dff = qd.dff_sd_normed_response{end};
        presp = qd.probability_response{end};
    
        vali = qd.types.green;
        mu_resp =  nanmean(presp(1:max_d_learni, vali));
        [irr sorti] = sort(mu_resp, 'descend');
        vali = vali(sorti(1:round(length(vali)*top_frac)));
        summ_dat_learni(a,:,1,1) = nanmean(presp(1:max_d_learni, vali)');
        summ_dat_learni(a,:,2,1) = nanmean(zsc_dff(1:max_d_learni, vali)');
 
        vali = qd.types.red;
        mu_resp =  nanmean(presp(1:max_d_learni, vali));
        [irr sorti] = sort(mu_resp, 'descend');
        vali = vali(sorti(1:round(length(vali)*top_frac)));
        summ_dat_learni(a,:,1,2) = nanmean(presp(1:max_d_learni, vali)');
        summ_dat_learni(a,:,2,2) = nanmean(zsc_dff(1:max_d_learni, vali)');
    end    

    max_d_nlearni = 8;
    summ_dat_nlearni = nan*zeros(length(all_dat.settings.nlearni), max_d_nlearni, 2, 2); % anim , day , presp/zsc dff, ops+/ops-
    for a=1:length(all_dat.settings.nlearni)
        ai = all_dat.settings.nlearni(a);
        qd = all_dat.quick_dat{ai};

        zsc_dff = qd.dff_sd_normed_response{end};
        presp = qd.probability_response{end};
    
        vali = qd.types.green;
        mu_resp =  nanmean(presp(1:max_d_nlearni, vali));
        [irr sorti] = sort(mu_resp, 'descend');
        vali = vali(sorti(1:round(length(vali)*top_frac)));
        summ_dat_nlearni(a,:,1,1) = nanmean(presp(1:max_d_nlearni, vali)');
        summ_dat_nlearni(a,:,2,1) = nanmean(zsc_dff(1:max_d_nlearni, vali)');
 
        vali = qd.types.red;
        mu_resp =  nanmean(presp(1:max_d_nlearni, vali));
        [irr sorti] = sort(mu_resp, 'descend');
        vali = vali(sorti(1:round(length(vali)*top_frac)));        
        summ_dat_nlearni(a,:,1,2) = nanmean(presp(1:max_d_nlearni, vali)');
        summ_dat_nlearni(a,:,2,2) = nanmean(zsc_dff(1:max_d_nlearni, vali)');
    end    

    maxval_stats_ax = [.15 .4 .25 1.2];
    for x=1:4
        hold(day_by_day_ax(x), 'on');
        aa = axis(stats_ax(x));
        axis(stats_ax(x), [aa(1) aa(2) aa(3) maxval_stats_ax(x)]);
    end

    summ_dat_learni(find(summ_dat_learni == 0)) = nan;
    summ_dat_nlearni(find(summ_dat_nlearni == 0)) = nan;

    if (0) % PRE POST PLOTS
        pre_di = 1:3;
        post_di = 6:8;

        lw_t = 3;

        delta_learni = nanmean(squeeze(summ_dat_learni(:,post_di,1,1))')- nanmean(squeeze(summ_dat_learni(:,pre_di,1,1))');
        delta_nlearni = nanmean(squeeze(summ_dat_blearni(:,post_di,1,1))')- nanmean(squeeze(summ_dat_nlearni(:,pre_di,1,1))');
        


    elseif (0) % individual plots?
        lw_t = 3;

        for a=1:length(all_dat.settings.learni) ; plot (day_by_day_ax(1), 1:max_d_learni, squeeze(summ_dat_learni(:,:,1,1)), 'Color', min([1 1 1], [.5 .5 .5]+learn_col));end
        plot(day_by_day_ax(1), 1:max_d_learni, nanmean(squeeze(summ_dat_learni(:,:,1,1))), 'Color', learn_col, 'LineWidth', lw_t);
        for a=1:length(all_dat.settings.learni) ; plot (day_by_day_ax(1), 1:max_d_nlearni, squeeze(summ_dat_nlearni(:,:,1,1)), 'Color', min([1 1 1], [.75 .75 .75]+nlearn_col));end
        plot(day_by_day_ax(1), 1:max_d_nlearni, nanmean(squeeze(summ_dat_nlearni(:,:,1,1))), 'Color', nlearn_col, 'LineWidth', lw_t);

        for a=1:length(all_dat.settings.learni) ; plot (day_by_day_ax(2), 1:max_d_learni, squeeze(summ_dat_learni(:,:,2,1)), 'Color', min([1 1 1], [.5 .5 .5]+learn_col));end
        plot(day_by_day_ax(2), 1:max_d_learni, nanmean(squeeze(summ_dat_learni(:,:,2,1))), 'Color', learn_col, 'LineWidth', lw_t);
        for a=1:length(all_dat.settings.learni) ; plot (day_by_day_ax(2), 1:max_d_nlearni, squeeze(summ_dat_nlearni(:,:,2,1)), 'Color', min([1 1 1], [.75 .75 .75]+nlearn_col));end
        plot(day_by_day_ax(2), 1:max_d_nlearni, nanmean(squeeze(summ_dat_nlearni(:,:,2,1))), 'Color', nlearn_col, 'LineWidth', lw_t);

        for a=1:length(all_dat.settings.learni) ; plot (day_by_day_ax(3), 1:max_d_learni, squeeze(summ_dat_learni(:,:,1,2)), 'Color', min([1 1 1], [.5 .5 .5]+learn_col));end
        plot(day_by_day_ax(3), 1:max_d_learni, nanmean(squeeze(summ_dat_learni(:,:,1,2))), 'Color', learn_col, 'LineWidth', lw_t);
        for a=1:length(all_dat.settings.learni) ; plot (day_by_day_ax(3), 1:max_d_nlearni, squeeze(summ_dat_nlearni(:,:,1,2)), 'Color', min([1 1 1], [.75 .75 .75]+nlearn_col));end
        plot(day_by_day_ax(3), 1:max_d_nlearni, nanmean(squeeze(summ_dat_nlearni(:,:,1,2))), 'Color', nlearn_col, 'LineWidth', lw_t);

        for a=1:length(all_dat.settings.learni) ; plot (day_by_day_ax(4), 1:max_d_learni, squeeze(summ_dat_learni(:,:,2,2)), 'Color', min([1 1 1], [.5 .5 .5]+learn_col));end
        plot(day_by_day_ax(4), 1:max_d_learni, nanmean(squeeze(summ_dat_learni(:,:,2,2))), 'Color', learn_col, 'LineWidth', lw_t);
        for a=1:length(all_dat.settings.learni) ; plot (day_by_day_ax(4), 1:max_d_nlearni, squeeze(summ_dat_nlearni(:,:,2,2)), 'Color', min([1 1 1], [.75 .75 .75]+nlearn_col));end
        plot(day_by_day_ax(4), 1:max_d_nlearni, nanmean(squeeze(summ_dat_nlearni(:,:,2,2))), 'Color', nlearn_col, 'LineWidth', lw_t);

    else % aggregate
        axes(day_by_day_ax(1));
        plot_error_poly (1:max_d_learni, nanmean(squeeze(summ_dat_learni(:,:,1,1))), nanstd(squeeze(summ_dat_learni(:,:,1,1)))/sqrt(length(all_dat.settings.learni)), learn_col, [1 0.75 0.5]);    
        plot_error_poly (1:max_d_nlearni, nanmean(squeeze(summ_dat_nlearni(:,:,1,1))), nanstd(squeeze(summ_dat_nlearni(:,:,1,1)))/sqrt(length(all_dat.settings.nlearni)), nlearn_col, [1 1 1]*0.75);


        axes(day_by_day_ax(2));
        plot_error_poly (1:max_d_learni, nanmean(squeeze(summ_dat_learni(:,:,2,1))), nanstd(squeeze(summ_dat_learni(:,:,2,1)))/sqrt(length(all_dat.settings.learni)), learn_col, [1 0.75 0.5]);    
        plot_error_poly (1:max_d_nlearni, nanmean(squeeze(summ_dat_nlearni(:,:,2,1))), nanstd(squeeze(summ_dat_nlearni(:,:,2,1)))/sqrt(length(all_dat.settings.nlearni)), nlearn_col, [1 1 1]*0.75);

        axes(day_by_day_ax(3));
        plot_error_poly (1:max_d_learni, nanmean(squeeze(summ_dat_learni(:,:,1,2))), nanstd(squeeze(summ_dat_learni(:,:,1,2)))/sqrt(length(all_dat.settings.learni)), learn_col, [1 0.75 0.5]);    
        plot_error_poly (1:max_d_nlearni, nanmean(squeeze(summ_dat_nlearni(:,:,1,2))), nanstd(squeeze(summ_dat_nlearni(:,:,1,2)))/sqrt(length(all_dat.settings.nlearni)), nlearn_col, [1 1 1]*0.75);

        axes(day_by_day_ax(4));
        plot_error_poly (1:max_d_learni, nanmean(squeeze(summ_dat_learni(:,:,2,2))), nanstd(squeeze(summ_dat_learni(:,:,2,2)))/sqrt(length(all_dat.settings.learni)), learn_col, [1 0.75 0.5]);    
        plot_error_poly (1:max_d_nlearni, nanmean(squeeze(summ_dat_nlearni(:,:,2,2))), nanstd(squeeze(summ_dat_nlearni(:,:,2,2)))/sqrt(length(all_dat.settings.nlearni)), nlearn_col, [1 1 1]*0.75);

        for x=1:4
            aa = axis(day_by_day_ax(x));
            axis(day_by_day_ax(x), [0 aa(2)+1 0 aa(4)]); 
        end
    end

    beh_corr_plot (corr_ax(1), all_dat, [nanmean(squeeze(summ_dat_learni(:,stat_di,1,1)')) nanmean(squeeze(summ_dat_nlearni(:,stat_di,1,1)'))])
    beh_corr_plot (corr_ax(2), all_dat ,[nanmean(squeeze(summ_dat_learni(:,stat_di,2,1)')) nanmean(squeeze(summ_dat_nlearni(:,stat_di,2,1)'))])
    beh_corr_plot (corr_ax(3), all_dat, [nanmean(squeeze(summ_dat_learni(:,stat_di,1,2)')) nanmean(squeeze(summ_dat_nlearni(:,stat_di,1,2)'))])
    beh_corr_plot (corr_ax(4), all_dat, [nanmean(squeeze(summ_dat_learni(:,stat_di,2,2)')) nanmean(squeeze(summ_dat_nlearni(:,stat_di,2,2)'))])

    axes(day_by_day_ax(1));
    title(day_by_day_ax(1), 'ops-');
    ylabel(day_by_day_ax(1), 'Presp');

    hold(day_by_day_ax(2),'on');
    title(day_by_day_ax(2), 'ops-');
    ylabel(day_by_day_ax(2), 'zsc dff');

    hold(day_by_day_ax(3),'on');
    title(day_by_day_ax(3), 'ops+');
    ylabel(day_by_day_ax(3), 'Presp');

    hold(day_by_day_ax(4),'on');
    title(day_by_day_ax(4), 'ops+');
    ylabel(day_by_day_ax(4), 'zsc dff');

    ecol = [1 1 1]*0.75;
    hold(stats_ax(1),'on');
    plot(stats_ax(1), ones(1,length(all_dat.settings.learni)), nanmean(squeeze(summ_dat_learni(:,stat_di,1,1)')), 'o', 'Color', ecol, 'MarkerFaceColor', learn_col);
    plot(stats_ax(1), 1, nanmean(nanmean(squeeze(summ_dat_learni(:,stat_di,1,1)'))), 'o', 'Color', learn_col, 'MarkerSize', 12);
    plot(stats_ax(1), 1+ones(1,length(all_dat.settings.nlearni)), nanmean(squeeze(summ_dat_nlearni(:,stat_di,1,1)')), 'o', 'Color', ecol, 'MarkerFaceColor', nlearn_col);
    plot(stats_ax(1), 2, nanmean(nanmean(squeeze(summ_dat_nlearni(:,stat_di,1,1)'))), 'o', 'Color', nlearn_col, 'MarkerSize', 12);
    aa = axis(stats_ax(1));
    axis(stats_ax(1), [0 3 0 aa(4)]);
    [h pv] = ttest2(nanmean(squeeze(summ_dat_learni(:,stat_di,1,1)')), nanmean(squeeze(summ_dat_nlearni(:,stat_di,1,1)')));
    title(stats_ax(1), sprintf('p=%0.3f', pv));
    disp(sprintf('presp ops- L: %0.3f/%0.3f NL: %0.3f/%0.3f p: %0.3f', nanmean(nanmean(squeeze(summ_dat_learni(:,stat_di,1,1)'))), nanstd(nanmean(squeeze(summ_dat_learni(:,stat_di,1,1)'))), ...
                                                                       nanmean(nanmean(squeeze(summ_dat_nlearni(:,stat_di,1,1)'))), nanstd(nanmean(squeeze(summ_dat_nlearni(:,stat_di,1,1)'))), pv));

    hold(stats_ax(2),'on');
    plot(stats_ax(2), ones(1,length(all_dat.settings.learni)), nanmean(squeeze(summ_dat_learni(:,stat_di,2,1)')), 'o', 'Color', ecol, 'MarkerFaceColor', learn_col);
    plot(stats_ax(2), 1, nanmean(nanmean(squeeze(summ_dat_learni(:,stat_di,2,1)'))), 'o', 'Color', learn_col, 'MarkerSize', 12);
    plot(stats_ax(2), 1+ones(1,length(all_dat.settings.nlearni)), nanmean(squeeze(summ_dat_nlearni(:,stat_di,2,1)')), 'o', 'Color', ecol, 'MarkerFaceColor', nlearn_col);
    plot(stats_ax(2), 2, nanmean(nanmean(squeeze(summ_dat_nlearni(:,stat_di,2,1)'))), 'o', 'Color', nlearn_col, 'MarkerSize', 12);
    aa = axis(stats_ax(2));
    axis(stats_ax(2), [0 3 0 aa(4)]);
    [h pv] = ttest2(nanmean(squeeze(summ_dat_learni(:,stat_di,2,1)')), nanmean(squeeze(summ_dat_nlearni(:,stat_di,2,1)')));
    title(stats_ax(2), sprintf('p=%0.3f', pv));
    disp(sprintf('zsc dff ops- L: %0.3f/%0.3f NL: %0.3f/%0.3f p: %0.3f', nanmean(nanmean(squeeze(summ_dat_learni(:,stat_di,2,1)'))), nanstd(nanmean(squeeze(summ_dat_learni(:,stat_di,2,1)'))), ...
                                                                       nanmean(nanmean(squeeze(summ_dat_nlearni(:,stat_di,2,1)'))), nanstd(nanmean(squeeze(summ_dat_nlearni(:,stat_di,2,1)'))), pv));


    hold(stats_ax(3),'on');
    plot(stats_ax(3), ones(1,length(all_dat.settings.learni)), nanmean(squeeze(summ_dat_learni(:,stat_di,1,2)')), 'o', 'Color', ecol, 'MarkerFaceColor', learn_col);
    plot(stats_ax(3), 1, nanmean(nanmean(squeeze(summ_dat_learni(:,stat_di,1,2)'))), 'o', 'Color', learn_col, 'MarkerSize', 12);
    plot(stats_ax(3), 1+ones(1,length(all_dat.settings.nlearni)), nanmean(squeeze(summ_dat_nlearni(:,stat_di,1,2)')), 'o', 'Color', ecol, 'MarkerFaceColor', nlearn_col);
    plot(stats_ax(3), 2, nanmean(nanmean(squeeze(summ_dat_nlearni(:,stat_di,1,2)'))), 'o', 'Color', nlearn_col, 'MarkerSize', 12);
    aa = axis(stats_ax(3));
    axis(stats_ax(3), [0 3 0 aa(4)]);
    [h pv] = ttest2(nanmean(squeeze(summ_dat_learni(:,stat_di,1,2)')), nanmean(squeeze(summ_dat_nlearni(:,stat_di,1,2)')));
    title(stats_ax(3), sprintf('p=%0.3f', pv));
    disp(sprintf('presp ops+ L: %0.3f/%0.3f NL: %0.3f/%0.3f p: %0.3f', nanmean(nanmean(squeeze(summ_dat_learni(:,stat_di,1,2)'))), nanstd(nanmean(squeeze(summ_dat_learni(:,stat_di,1,2)'))), ...
                                                                       nanmean(nanmean(squeeze(summ_dat_nlearni(:,stat_di,1,2)'))), nanstd(nanmean(squeeze(summ_dat_nlearni(:,stat_di,1,2)'))), pv));


    hold(stats_ax(4),'on');
    plot(stats_ax(4), ones(1,length(all_dat.settings.learni)), nanmean(squeeze(summ_dat_learni(:,stat_di,2,2)')), 'o', 'Color', ecol, 'MarkerFaceColor', learn_col);
    plot(stats_ax(4), 1, nanmean(nanmean(squeeze(summ_dat_learni(:,stat_di,2,2)'))), 'o', 'Color', learn_col, 'MarkerSize', 12);
    plot(stats_ax(4), 1+ones(1,length(all_dat.settings.nlearni)), nanmean(squeeze(summ_dat_nlearni(:,stat_di,2,2)')), 'o', 'Color', ecol, 'MarkerFaceColor', nlearn_col);
    plot(stats_ax(4), 2, nanmean(nanmean(squeeze(summ_dat_nlearni(:,stat_di,2,2)'))), 'o', 'Color', nlearn_col, 'MarkerSize', 12);
    aa = axis(stats_ax(4));
    axis(stats_ax(4), [0 3 0 aa(4)]);
    [h pv] = ttest2(nanmean(squeeze(summ_dat_learni(:,stat_di,2,2)')), nanmean(squeeze(summ_dat_nlearni(:,stat_di,2,2)')));
    title(stats_ax(4), sprintf('p=%0.3f', pv));
    disp(sprintf('zscdff ops- L: %0.3f/%0.3f NL: %0.3f/%0.3f p: %0.3f', nanmean(nanmean(squeeze(summ_dat_learni(:,stat_di,2,2)'))), nanstd(nanmean(squeeze(summ_dat_learni(:,stat_di,2,2)'))), ...
                                                                       nanmean(nanmean(squeeze(summ_dat_nlearni(:,stat_di,2,2)'))), nanstd(nanmean(squeeze(summ_dat_nlearni(:,stat_di,2,2)'))), pv));


    % --- finalize axes
    for a=1:length(ax)
        set(ax(a), 'TickDir', 'out', 'FontSize', 15);
    end


