%
% For all the decoding goodies ...
%
function figure_7_decoding

    [all_dat settings] = get_all_data_pstim; 
    
    % --- setup plots

    fh = figure ('Position',[0 0 1200 600]);

    example_mouse_auc_red_ax = subplot('Position', [0.05 0.5 .15 .4]); 
    example_mouse_auc_red_zoom_ax = subplot('Position', [0.75 0.5 .15 .2]); 
    example_mouse_auc_green_ax = subplot('Position', [0.05 0.05 .15 .4]);
    example_mouse_auc_green_zoom_ax = subplot('Position', [0.75 0.05 .15 .2]);
    top_cells_auc_red_ax = subplot('Position', [0.25 0.5 .2 .2]); 
    top_cells_auc_green_ax = subplot('Position', [0.25 0.05 .2 .2]); 
    n_good_cells_auc_red_ax = subplot('Position', [0.5 0.5 .2 .2]); 
    n_good_cells_auc_green_ax = subplot('Position', [0.5 0.05 .2 .2]); 
    corr_drift_beh_red_ax = subplot('Position', [0.75 .75 .1 .2]);
    corr_drift_beh_green_ax = subplot('Position', [0.875 .75 .1 .2  ]);

    ax = [example_mouse_auc_red_ax example_mouse_auc_green_ax top_cells_auc_red_ax top_cells_auc_green_ax ...
          n_good_cells_auc_red_ax n_good_cells_auc_green_ax example_mouse_auc_red_zoom_ax example_mouse_auc_green_zoom_ax];

    fh = figure ('Position',[0 0 1200 600]);

    example_mouse_auc_response_distro_ax = subplot('Position', [0.05 0.5 .15 .3]); 
    example_mouse_auc_respose_stim_hist_ax = subplot('Position', [0.3 0.5 .15 .3]); 
    example_mouse_auc_ax = subplot('Position', [0.6 0.5 .15 .3]);

    ax = [ax example_mouse_auc_response_distro_ax  example_mouse_auc_respose_stim_hist_ax example_mouse_auc_ax];

    % --- for example mouse, show ops+ and ops- trend for AUC stim decode
    ex_ani = 3;
    vali = find(~isnan(sum(all_dat.quick_dat{ex_ani}.discrim_specific_mat{1}')));
    dM = all_dat.quick_dat{ex_ani}.discrim_specific_mat{1}(vali,:);
    dM(find(dM < 0.5)) = 1-dM(find(dM<0.5));

    dMr = dM(:,all_dat.quick_dat{ex_ani}.types.red);
    [irr sorti] = sort(nanmean(dMr),'descend');
    imagesc(example_mouse_auc_red_ax, dMr(:,sorti)', [0.5 1]);
    title(example_mouse_auc_red_ax, 'ops+, range 0.5 1, AUC');
    colormap(example_mouse_auc_red_ax, colormap_human);
    imagesc(example_mouse_auc_red_zoom_ax, dMr(:,sorti(1:45))', [0.5 1]);
    colormap(example_mouse_auc_red_zoom_ax, colormap_human);

    dMg = dM(:,all_dat.quick_dat{ex_ani}.types.green);
    [irr sorti] = sort(nanmean(dMg),'descend');
    imagesc(example_mouse_auc_green_ax, dMg(:,sorti)', [0.5 1]);
    title(example_mouse_auc_green_ax, 'ops-, range 0.5 1, AUC');
    colormap(example_mouse_auc_green_ax, colormap_human);
    imagesc(example_mouse_auc_green_zoom_ax, dMg(:,sorti(1:118))', [0.5 1]);
    colormap(example_mouse_auc_green_zoom_ax, colormap_human);

    % --- trend in top 5/1% of AUC for ops+ and ops- ; learning and nonlearning on each
    learn_col = [1 0.25 1];
    learni = all_dat.settings.learni;
    nlearn_col = [1 1 1]*0.25;
    nlearni = all_dat.settings.nlearni;

    top_frac = .05; % .05 = top 5% ; change accordingly
    top_frac = .25; % .05 = top 5% ; change accordingly
    disp(sprintf('DECODING ANALYSIS USING TOP %d PERCENT OF CELLS', top_frac*100))

    ldi=1:22;
    nldi=1:8;
    n_ani = length(all_dat.quick_dat);
    top_auc_red = nan*zeros(n_ani, length(ldi));
    n_high_red = nan*zeros(n_ani, length(ldi));
    frac_high_red = nan*zeros(n_ani, length(ldi));
    top_auc_green = nan*zeros(n_ani, length(ldi));
    n_high_green = nan*zeros(n_ani, length(ldi));
    frac_high_green = nan*zeros(n_ani, length(ldi));
    for a=1:n_ani
        dM = all_dat.quick_dat{a}.discrim_specific_mat{1};
        dM(find(dM < 0.5)) = 1-dM(find(dM<0.5));

        dMr = dM(:,all_dat.quick_dat{a}.types.red);
        dMg = dM(:,all_dat.quick_dat{a}.types.green);

        if (ismember(a, learni)) ; di =ldi ; else ; di = nldi ; end

        for dd=1:length(di)
            d = di(dd);
            n_high_red(a,d) = length(find(dMr(d,:) > 0.75));
            frac_high_red(a,d) = length(find(dMr(d,:) > 0.75))/length(all_dat.quick_dat{a}.types.red);
            n_high_green(a,d) = length(find(dMg(d,:) > 0.75));
            frac_high_green(a,d) = length(find(dMg(d,:) > 0.75))/length(all_dat.quick_dat{a}.types.green);
        end

        [irr sorti] = sort(nanmean(dMr),'descend');
        nr = round(top_frac*size(dMr,2));
        top_auc_red(a,di) = nanmean(dMr(di,sorti(1:nr))');

        [irr sorti] = sort(nanmean(dMg),'descend');
        ng = round(top_frac*size(dMg,2));
        top_auc_green(a,di) = nanmean(dMg(di,sorti(1:ng))');
        
        if (a == ex_ani) ; disp(sprintf('For example animal, n_r=%d and n_g=%d', nr, ng)); end
    end
    nldi = di;

    di_compare = 1:3;
    di_compare = 7:8;

    axes(top_cells_auc_red_ax);
    plot_error_poly (ldi, nanmean(top_auc_red(learni,ldi)), nanstd(top_auc_red(learni,ldi))/sqrt(length(learni)), learn_col, [1 0.75 1]);    
    plot_error_poly (nldi, nanmean(top_auc_red(nlearni,nldi)), nanstd(top_auc_red(nlearni,nldi))/sqrt(length(nlearni)), nlearn_col, [1 1 1]*0.75);  
    xlabel(sprintf('Day; stats %d to %d', di_compare(1), di_compare(end)));
    ylabel('AUC, stimulus decoding');
    axis(top_cells_auc_red_ax, [0 23 0.5 1]);
    plot_comparison(top_cells_auc_red_ax, top_auc_red, learni, nlearni, di_compare, 'Ops+ AUC');

    beh_corr_plot (corr_drift_beh_red_ax, all_dat, [nanmean(top_auc_red(learni,di_compare)') nanmean(top_auc_red(nlearni,di_compare)') ], 'o+');

    axes(top_cells_auc_green_ax);
    plot_error_poly (ldi, nanmean(top_auc_green(learni,ldi)), nanstd(top_auc_green(learni,ldi))/sqrt(length(learni)), learn_col, [1 0.75 1]);    
    plot_error_poly (nldi, nanmean(top_auc_green(nlearni, nldi)), nanstd(top_auc_green(nlearni,nldi))/sqrt(length(nlearni)), nlearn_col, [1 1 1]*0.75);    
    axis(top_cells_auc_green_ax, [0 23 0.5 0.75]);
    plot_comparison(top_cells_auc_green_ax, top_auc_green, learni, nlearni, di_compare, 'Ops- AUC');

    beh_corr_plot (corr_drift_beh_green_ax, all_dat, [nanmean(top_auc_green(learni,di_compare)') nanmean(top_auc_green(nlearni,di_compare)') ], 'o-');

    % --- trend in number of cells with AUC > .75 for ops+ and ops- ; learning and nonlearning on each
    if (0) % n good cells
        axes(n_good_cells_auc_red_ax);
        plot_error_poly (ldi, nanmean(n_high_red(learni,ldi)), nanstd(n_high_red(learni,ldi))/sqrt(length(learni)), learn_col, [1 0.75 1]);    
        plot_error_poly (nldi, nanmean(n_high_red(nlearni,nldi)), nanstd(n_high_red(nlearni,nldi))/sqrt(length(nlearni)), nlearn_col, [1 1 1]*0.75);  
        xlabel('Day');
        ylabel('AUC, stimulus decoding');
        plot_comparison(n_good_cells_auc_red_ax, n_high_red, learni, nlearni, di_compare, 'Ops+ ncelss');
        axis(n_good_cells_auc_red_ax, [0 23 0 150]);

        axes(n_good_cells_auc_green_ax);
        plot_error_poly (ldi, nanmean(n_high_green(learni,ldi)), nanstd(n_high_green(learni,ldi))/sqrt(length(learni)), learn_col, [1 0.75 1]);    
        plot_error_poly (nldi, nanmean(n_high_green(nlearni,nldi)), nanstd(n_high_green(nlearni,nldi))/sqrt(length(nlearni)), nlearn_col, [1 1 1]*0.75);    
        plot_comparison(n_good_cells_auc_green_ax, n_high_green, learni, nlearni, di_compare, 'Ops- ncells');
        axis(n_good_cells_auc_green_ax, [0 23 0 65]);
    else % frac good cells (AUC > .75)
        axes(n_good_cells_auc_red_ax);
        plot_error_poly (ldi, nanmean(frac_high_red(learni,ldi)), nanstd(frac_high_red(learni,ldi))/sqrt(length(learni)), learn_col, [1 0.75 1]);    
        plot_error_poly (nldi, nanmean(frac_high_red(nlearni,nldi)), nanstd(frac_high_red(nlearni,nldi))/sqrt(length(nlearni)), nlearn_col, [1 1 1]*0.75);  
        xlabel('Day');
        ylabel('AUC, stimulus decoding');
        plot_comparison(n_good_cells_auc_red_ax, frac_high_red, learni, nlearni, di_compare, 'Ops+ frac');
        axis(n_good_cells_auc_red_ax, [0 23 0 .25]);

        axes(n_good_cells_auc_green_ax);
        plot_error_poly (ldi, nanmean(frac_high_green(learni,ldi)), nanstd(frac_high_green(learni,ldi))/sqrt(length(learni)), learn_col, [1 0.75 1]);    
        plot_error_poly (nldi, nanmean(frac_high_green(nlearni,nldi)), nanstd(frac_high_green(nlearni,nldi))/sqrt(length(nlearni)), nlearn_col, [1 1 1]*0.75);    
        plot_comparison(n_good_cells_auc_green_ax, frac_high_green, learni, nlearni, di_compare, 'Ops- frac');
        axis(n_good_cells_auc_green_ax, [0 23 0 .04]);    
    end

    % --- example walking reader through how AUC is calculated for stim responses
    ex_sess = 22;
    ex_cell_id = 10040535; % good but not best - so nice example
    ex_ci = find(all_dat.quick_dat{ex_ani}.ids == ex_cell_id);

    load([settings.base_dir filesep 'pertrial_summary_example_mouse.mat'])
    ddff = pt_dat.session(ex_sess).dff_stim_epoch - pt_dat.session(ex_sess).dff_pre_stim_epoch;

    resp_vec = ddff(ex_ci,:);
    stim_vec = pt_dat.session(ex_sess).num_pulses;

    hold(example_mouse_auc_response_distro_ax, 'on');

    ms = 5;
    ii = find(stim_vec == 1 | stim_vec == 3);
    rnd = rand(1,length(ii))-0.5;
    plot(example_mouse_auc_response_distro_ax, stim_vec(ii)+rnd, resp_vec(ii), 'o', 'MarkerEdge', 'None', 'MarkerFaceColor', [0 0 1], 'MarkerSize', ms);
    ii = find(stim_vec == 5);
    rnd = rand(1,length(ii))-0.5;
    plot(example_mouse_auc_response_distro_ax, stim_vec(ii)+rnd, resp_vec(ii), 'o', 'MarkerEdge', 'None', 'MarkerFaceColor', [0 0 0], 'MarkerSize', ms);
    ii = find(stim_vec == 7 | stim_vec == 9);
    rnd = rand(1,length(ii))-0.5;
    plot(example_mouse_auc_response_distro_ax, stim_vec(ii)+rnd, resp_vec(ii), 'o', 'MarkerEdge', 'None', 'MarkerFaceColor', [1 0 0], 'MarkerSize', ms);
    axis(example_mouse_auc_response_distro_ax, [0 10 -1 2]);
    xlabel(example_mouse_auc_response_distro_ax, 'Pulse count');
    ylabel(example_mouse_auc_response_distro_ax, 'Evoked dff');
    
    hold(example_mouse_auc_respose_stim_hist_ax, 'on');
    ii = find(stim_vec == 1 | stim_vec == 3);
    resp_low_vec = resp_vec(ii);
    histogram(example_mouse_auc_respose_stim_hist_ax, resp_vec(ii), -1:.1:2, 'EdgeColor', 'none', 'FaceColor',[0 0 1], 'FaceAlpha', 0.5);
    ii = find(stim_vec == 7 | stim_vec == 9);
    resp_hi_vec = resp_vec(ii);
    histogram(example_mouse_auc_respose_stim_hist_ax, resp_vec(ii), -1:.1:2, 'EdgeColor', 'none', 'FaceColor',[1 0 0], 'FaceAlpha', 0.5);
    axis(example_mouse_auc_respose_stim_hist_ax, [-1 2 0 35]);
    xlabel(example_mouse_auc_response_distro_ax, 'Evoked dff');
    ylabel(example_mouse_auc_response_distro_ax, 'Count');

    thresh_vec = -1:.1:2;
    hr = nan*thresh_vec;
    far = nan*thresh_vec;
    for t=1:length(thresh_vec)
        n_hit_t = length(find(resp_hi_vec > thresh_vec(t)));
        n_miss_t = length(find(resp_hi_vec < thresh_vec(t)));
        hr(t) = n_hit_t/(n_hit_t+n_miss_t);

        n_fa_t = length(find(resp_low_vec > thresh_vec(t)));
        n_cr_t = length(find(resp_low_vec < thresh_vec(t)));
        far(t) = n_fa_t/(n_fa_t+n_cr_t);
    end
    hold(example_mouse_auc_ax, 'on');
    plot(example_mouse_auc_ax, far, hr, 'k-', 'LineWidth', 3);
    plot(example_mouse_auc_ax, far, hr, 'o', 'Color','none', 'MarkerSize', 5, 'MarkerFaceColor', [1 1 1]*0.75);
    plot(example_mouse_auc_ax, [0 1], [0 1], 'k:', 'LineWidth', 1);
    axis(example_mouse_auc_ax, [-0.1 1.1 -0.1 1.1]);
    xlabel(example_mouse_auc_ax, 'False alarm rate');
    ylabel(example_mouse_auc_ax, 'Hit rate');

    % --- finalize axes
    for a=1:length(ax)
        set(ax(a), 'TickDir', 'out', 'FontSize', 15);
    end                          

function plot_comparison (ax, dat_mat, learni, nlearni, di, tstr)
    learn_col = [1 0.25 1];
    nlearn_col = [1 1 1]*0.25;

    plot(ax, max(di)+1*ones(1,length(learni)),  nanmean(dat_mat(learni,di)'), 'o', 'Color', 'none', 'MarkerFaceColor', [1 .5 1]);
    plot(ax, max(di)+1,  nanmean(nanmean(dat_mat(learni,di)')), 'o', 'Color', learn_col, 'MarkerSize', 10, 'LineWidth', 2);
    plot(ax, max(di)+2*ones(1,length(nlearni)),  nanmean(dat_mat(nlearni,di)'), 'o', 'Color', 'none', 'MarkerFaceColor', [1 1 1]*0.5);
    plot(ax, max(di)+2,  nanmean(nanmean(dat_mat(nlearni,di)')), 'o', 'Color', nlearn_col, 'MarkerSize', 10, 'LineWidth', 2);

    [h p] = ttest2(nanmean(dat_mat(learni,di)'), nanmean(dat_mat(nlearni,di)'));
    title(ax, sprintf('%s L v nL: %0.3f', tstr, p));

    disp(sprintf('%s L v. nL p: %0.3f L mu/sd: %0.3f/%0.3f nL: %0.3f/%0.3f days: %d to %d', tstr, p,  ...
            nanmean(nanmean(dat_mat(learni,di)')),  nanstd(nanstd(dat_mat(learni,di)')), ...
            nanmean(nanmean(dat_mat(nlearni,di)')),  nanstd(nanstd(dat_mat(nlearni,di)')), di(1), di(end)));

    dis = [1:3];
    disp(sprintf('%s L mu/sd: %0.3f/%0.3f nL: %0.3f/%0.3f  days: %d to %d', tstr,  ...
            nanmean(nanmean(dat_mat(learni,dis)')),  nanstd(nanstd(dat_mat(learni,dis)')), ...
            nanmean(nanmean(dat_mat(nlearni,dis)')),  nanstd(nanstd(dat_mat(nlearni,dis)')), dis(1), dis(end)));
    dif = [20:22];
    [h p] = ttest2(nanmean(dat_mat(learni,dif)'), nanmean(dat_mat(learni,dis)'));
    disp(sprintf('  %s d1-3 v d20-22 p: %0.3f L mu/sd: %0.3f/%0.3f nL: %0.3f/%0.3f days: %d to %d', tstr, p,  ...
            nanmean(nanmean(dat_mat(learni,dif)')),  nanstd(nanstd(dat_mat(learni,dif)')), ...
            nanmean(nanmean(dat_mat(nlearni,dif)')),  nanstd(nanstd(dat_mat(nlearni,dif)')), dif(1), dif(end)));

