function figure_5_intraday_stability (force_redo)
    if (nargin < 1) ; force_redo = 0 ; end
    
    % --- get data
    [all_dat settings] = get_all_data_pstim;
    rep_drift_dat = get_representational_drift_data(all_dat, settings, force_redo);

    % --- setup figures
    ax = [];  

    fh = figure ('Position',[0 0 1200 600]);

    for i=1:3
        example_session_response_matrix_ax(i) = subplot('Position', [0.05 0.1+.33*(i-1) .1 .2]);
        ax(end+1) = example_session_response_matrix_ax(i);
    end
    example_intrasession_t_to_t_corr_mat_r_ax = subplot('Position', [0.2 0.5 .2 .4]); % t_to_t = trial to trial
    example_intrasession_t_to_t_corr_mat_r_avg_ax = subplot('Position', [0.2 0.35 .2 .1]);
    example_intrasession_t_to_t_corr_mat_g_ax = subplot('Position', [0.45 0.5 .2 .4]);
    example_intrasession_t_to_t_corr_mat_g_avg_ax = subplot('Position', [0.45 0.35 .2 .1]);
    example_intrasession_t_to_t_corr_cb_ax = subplot('Position', [0.675 0.5 .025 .4]);

    intrasession_t_to_t_progression_by_mouse_r_ax = subplot('Position', [0.725 0.7 .25 .2]);
    intrasession_t_to_t_progression_by_mouse_g_ax = subplot('Position', [0.725 0.45 .25 .2]);

    corr_intrasession_var_beh_g_ax = subplot('Position', [0.7 0.05 .125 .25]);
    corr_intrasession_var_beh_r_ax = subplot('Position', [0.85 0.05 .125 .25]);

    ax = [ ax example_intrasession_t_to_t_corr_mat_r_ax  example_intrasession_t_to_t_corr_mat_r_avg_ax ...
           example_intrasession_t_to_t_corr_mat_g_ax example_intrasession_t_to_t_corr_mat_g_avg_ax ...
           corr_intrasession_var_beh_g_ax corr_intrasession_var_beh_r_ax ...
           example_intrasession_t_to_t_corr_cb_ax intrasession_t_to_t_progression_by_mouse_r_ax intrasession_t_to_t_progression_by_mouse_g_ax];
     
    fh = figure ('Position',[0 0 1200 600]);

    for i=1:3
        example_cell_ax(i) = subplot('Position', [0.05 0.1+.33*(i-1) .3 .2]);
        ax(end+1) = example_cell_ax(i);
    end
    example_session_r_ax = subplot('Position', [0.4 0.1 .15 .4]); % t_to_t = trial to trial
    example_session_r_good_ax = subplot('Position', [0.4 0.55 .15 .4]); % t_to_t = trial to trial
    example_session_g_ax = subplot('Position', [0.65 0.1 .15 .4]); % t_to_t = trial to trial
    example_session_g_good_ax = subplot('Position', [0.65 0.55 .15 .4]); % t_to_t = trial to trial

    ax = [ ax example_session_g_ax example_session_r_good_ax example_session_r_ax example_session_g_good_ax];
      
 
    %% 1) example response over the course of a (late) session
    ex_ani = 3;
    ex_sess = 12;

    load([settings.base_dir filesep 'pertrial_summary_example_mouse.mat'])
    
    ninei = find(pt_dat.session(ex_sess).num_pulses == 9);
    greeni = all_dat.quick_dat{ex_ani}.types.green;
    redi = all_dat.quick_dat{ex_ani}.types.red;
    ddff = pt_dat.session(ex_sess).dff_stim_epoch - pt_dat.session(ex_sess).dff_pre_stim_epoch;
    vali = find(nanmean( pt_dat.session(ex_sess).dff_sig_response' ) > 0.1)';
    ddff = ddff(vali,:);
    greeni = find(ismember(vali,greeni));
    redi = find(ismember(vali,redi));
    
    [irr gssi] = sort(nanmean(ddff(greeni,ninei)'), 'descend');
    [irr rssi] = sort(nanmean(ddff(redi,ninei)'), 'descend');

    ii = [5 50 83 90];
    for i=1:3
        v1 =  ddff(greeni,ninei(ii(1)));
        v1 = v1/quantile(v1,.95);
        v2 =  ddff(greeni,ninei(ii(i+1)));
        v2 = v2/quantile(v2,.95);
        plot(example_session_response_matrix_ax(i), v1,v2, 'o', 'Color', [1 1 1]*0.5, 'MarkerFaceColor', settings.colors.ops_negative, 'MarkerSize',7); 
        hold(example_session_response_matrix_ax(i), 'on');
        plot(example_session_response_matrix_ax(i), [-0.1 1], [-0.1 1], 'k:');
        xlabel(example_session_response_matrix_ax(i), sprintf ('Nrmd trial %d response', ii(1)));
        ylabel(example_session_response_matrix_ax(i), sprintf ('Nrmd trial %d response', ii(i+1)));
        [cv pv] = nancorr(v1,v2);
        title(example_session_response_matrix_ax(i), sprintf('Corr: %0.2f p: %0.3f', cv, pv));
        axis(example_session_response_matrix_ax(i), [-0.2 1.5 -0.2 1.5]);
    end

    %% 2) example ops+/ops- trial pairwise corr mat 
    cmr = rep_drift_dat.ani(ex_ani).cmr{ex_sess};
    cmg = rep_drift_dat.ani(ex_ani).cmg{ex_sess};
    goodi = find(nansum(cmr) > 0 & nansum(cmg) > 0); % eliminate weird trials

    M = .6;
    imagesc(example_intrasession_t_to_t_corr_mat_r_ax, cmr(goodi,goodi), [0 M]);
    colormap(example_intrasession_t_to_t_corr_mat_r_ax, colormap_human);
    imagesc(example_intrasession_t_to_t_corr_mat_g_ax, cmg(goodi,goodi), [0 M]);
    colormap(example_intrasession_t_to_t_corr_mat_g_ax, colormap_human);

    title(example_intrasession_t_to_t_corr_mat_r_ax, 'Example corr, ops+, 9p only');

    title(example_intrasession_t_to_t_corr_mat_g_ax, ['Example corr, ops-, 9p nt: ' num2str(length(goodi))]);
 
    cb = colorbar(example_intrasession_t_to_t_corr_cb_ax);
    colormap(example_intrasession_t_to_t_corr_cb_ax, colormap_human);
    cb.Position = [0.675 0.5 0.025 0.4];
    cb.Ticks = [0 1];
    cb.TickLabels = {'0', num2str(M)};
    set(example_intrasession_t_to_t_corr_cb_ax, 'XTick',[], 'YTick', []);

    %% 3) unidim version to compare both - e.g., mean within 10 trials of corr
    muval_r = nan*zeros(1,length(goodi));
    muval_g = nan*zeros(1,length(goodi));
    for t=5:length(goodi)-5
        muval_r(t) = nanmean(reshape(cmr(goodi(t-4:t+4), goodi(t-4:t+4)),[],1));
        muval_g(t) = nanmean(reshape(cmg(goodi(t-4:t+4), goodi(t-4:t+4)),[],1));
    end

    plot(example_intrasession_t_to_t_corr_mat_r_avg_ax, muval_r, '-', 'Color', settings.colors.ops_positive, 'LineWidth', 3);
    xlabel(example_intrasession_t_to_t_corr_mat_r_avg_ax, 'Trial');
    ylabel(example_intrasession_t_to_t_corr_mat_r_avg_ax, 'Mean correlation');
    aa=axis(example_intrasession_t_to_t_corr_mat_r_avg_ax) ; 
    axis(example_intrasession_t_to_t_corr_mat_r_avg_ax, [aa(1) aa(2) 0.2 0.7]);

    plot(example_intrasession_t_to_t_corr_mat_g_avg_ax, muval_r, '-', 'Color', [1 1 1]*0.8, 'LineWidth', 3);
    hold(example_intrasession_t_to_t_corr_mat_g_avg_ax, 'on');
    plot(example_intrasession_t_to_t_corr_mat_g_avg_ax, muval_g, '-', 'Color', settings.colors.ops_negative, 'LineWidth', 3);
    set(example_intrasession_t_to_t_corr_mat_g_avg_ax, 'TickDir', 'out', 'FontSize', 15);
    aa=axis(example_intrasession_t_to_t_corr_mat_g_avg_ax) ; 
    axis(example_intrasession_t_to_t_corr_mat_g_avg_ax, [aa(1) aa(2) 0.2 0.7]);

    %% 4) within-day similarity is higher among learners in ops- but not ops+
    learn_col = [1 0.25 1];
    learni = settings.learni;
    nlearn_col = [1 1 1]*0.25;
    nlearni = settings.nlearni;
    di = 1:8;

    rep_drift_dat.all_mur(find(rep_drift_dat.all_mur == 0)) = nan;
    rep_drift_dat.all_mug(find(rep_drift_dat.all_mug == 0)) = nan;

    axes(intrasession_t_to_t_progression_by_mouse_r_ax);
    plot_error_poly (di, nanmean(rep_drift_dat.all_mur(learni,di)), nanstd(rep_drift_dat.all_mur(learni,di))/sqrt(length(learni)), learn_col, [1 0.75 1]);    
    plot_error_poly (di, nanmean(rep_drift_dat.all_mur(nlearni,di)), nanstd(rep_drift_dat.all_mur(nlearni,di))/sqrt(length(nlearni)), nlearn_col, [1 1 1]*.75);
    axis(intrasession_t_to_t_progression_by_mouse_r_ax, [0 max(di)+3 -0.01 0.8]);
    set(intrasession_t_to_t_progression_by_mouse_r_ax, 'XTick',[]);


    plot(intrasession_t_to_t_progression_by_mouse_r_ax, max(di)+1*ones(1,length(learni)),  nanmean(rep_drift_dat.all_mur(learni,di)'), 'o', 'Color', 'none', 'MarkerFaceColor', [1 .5 1]);
    plot(intrasession_t_to_t_progression_by_mouse_r_ax, max(di)+1,  nanmean(nanmean(rep_drift_dat.all_mur(learni,di)')), 'o', 'Color', learn_col, 'MarkerSize', 10, 'LineWidth', 2);
    plot(intrasession_t_to_t_progression_by_mouse_r_ax, max(di)+2*ones(1,length(nlearni)),  nanmean(rep_drift_dat.all_mur(nlearni,di)'), 'o', 'Color', 'none', 'MarkerFaceColor', [1 1 1]*0.5);
    plot(intrasession_t_to_t_progression_by_mouse_r_ax, max(di)+2,  nanmean(nanmean(rep_drift_dat.all_mur(nlearni,di)')), 'o', 'Color', nlearn_col, 'MarkerSize', 10, 'LineWidth', 2);

    beh_corr_plot (corr_intrasession_var_beh_r_ax, all_dat, [nanmean(rep_drift_dat.all_mur(learni,di)') nanmean(rep_drift_dat.all_mur(nlearni,di)')], 'o+');

    [h p ] = ttest2(nanmean(rep_drift_dat.all_mur(learni,di)'),nanmean(rep_drift_dat.all_mur(nlearni,di)'));
    title(intrasession_t_to_t_progression_by_mouse_r_ax, sprintf('o+ learn v. nonlearn: %0.3f', p));
    ml = nanmean(nanmean(rep_drift_dat.all_mur(learni,di)));
    mn = nanmean(nanmean(rep_drift_dat.all_mur(nlearni,di)));
    sdl = nanstd(nanmean(rep_drift_dat.all_mur(learni,di)));
    sdn = nanstd(nanmean(rep_drift_dat.all_mur(nlearni,di))); 
    disp(sprintf('o+ learn v. nonlearn: %0.3f learn mu/sd: %0.3f / %0.3f ; nonlearn: %0.3f / %0.3f', p, ml, sdl, mn, sdn));


    hold(intrasession_t_to_t_progression_by_mouse_g_ax, 'on');
    axes(intrasession_t_to_t_progression_by_mouse_g_ax);
    plot_error_poly (di, nanmean(rep_drift_dat.all_mug(learni,di)), nanstd(rep_drift_dat.all_mug(learni,di))/sqrt(length(learni)), learn_col, [1 0.75 1]);    
    plot_error_poly (di, nanmean(rep_drift_dat.all_mug(nlearni,di)), nanstd(rep_drift_dat.all_mug(nlearni,di))/sqrt(length(nlearni)), nlearn_col, [1 1 1]*.75);
    axis(intrasession_t_to_t_progression_by_mouse_g_ax, [0 max(di)+3 -0.01 0.7]);
    
    plot(intrasession_t_to_t_progression_by_mouse_g_ax, max(di)+1*ones(1,length(learni)),  nanmean(rep_drift_dat.all_mug(learni,di)'), 'o', 'Color', 'none', 'MarkerFaceColor', [1 .5 1]);
    plot(intrasession_t_to_t_progression_by_mouse_g_ax, max(di)+1,  nanmean(nanmean(rep_drift_dat.all_mug(learni,di)')), 'o', 'Color', learn_col, 'MarkerSize', 10, 'LineWidth', 2);
    plot(intrasession_t_to_t_progression_by_mouse_g_ax, max(di)+2*ones(1,length(nlearni)),  nanmean(rep_drift_dat.all_mug(nlearni,di)'), 'o', 'Color', 'none', 'MarkerFaceColor', [1 1 1]*0.5);
    plot(intrasession_t_to_t_progression_by_mouse_g_ax, max(di)+2,  nanmean(nanmean(rep_drift_dat.all_mug(nlearni,di)')), 'o', 'Color', nlearn_col, 'MarkerSize', 10, 'LineWidth', 2);

    beh_corr_plot (corr_intrasession_var_beh_g_ax, all_dat, [nanmean(rep_drift_dat.all_mug(learni,di)') nanmean(rep_drift_dat.all_mug(nlearni,di)')], 'o-');

    [h p ] = ttest2(nanmean(rep_drift_dat.all_mug(learni,di)'),nanmean(rep_drift_dat.all_mug(nlearni,di)'));
    title(intrasession_t_to_t_progression_by_mouse_g_ax, sprintf('o- learn v. nonlearn: %0.3f', p));
    ml = nanmean(nanmean(rep_drift_dat.all_mug(learni,di)));
    mn = nanmean(nanmean(rep_drift_dat.all_mug(nlearni,di)));
    sdl = nanstd(nanmean(rep_drift_dat.all_mug(learni,di)));
    sdn = nanstd(nanmean(rep_drift_dat.all_mug(nlearni,di)));
    disp(sprintf('o- learn v. nonlearn: %0.3f learn mu/sd: %0.3f / %0.3f ; nonlearn: %0.3f / %0.3f', p, ml, sdl, mn, sdn));
    xlabel(intrasession_t_to_t_progression_by_mouse_g_ax, 'Training day');
   
    %% 6) example session
    ninei = find(pt_dat.session(ex_sess).num_pulses == 9);
    greeni = all_dat.quick_dat{ex_ani}.types.green;
    redi = all_dat.quick_dat{ex_ani}.types.red;
    ddff = pt_dat.session(ex_sess).dff_stim_epoch - pt_dat.session(ex_sess).dff_pre_stim_epoch;
    [irr gssi] = sort(nanmean(ddff(greeni,ninei)'), 'descend');
    [irr rssi] = sort(nanmean(ddff(redi,ninei)'), 'descend');

    imagesc(example_session_r_ax, ddff(redi(rssi),ninei),[0 .5]);
    colormap(example_session_r_ax, colormap_human);
    imagesc(example_session_g_ax, ddff(greeni(gssi),ninei),[0 .5]);
    colormap(example_session_g_ax, colormap_human);
    title(example_session_r_ax, ['ops+, max: 0.5 len: ' num2str(length(ninei))]);

    % best cells
    vali = find(nanmean( pt_dat.session(ex_sess).dff_sig_response' ) > 0.05)';
    ddff = ddff(vali,:);
    greeni = find(ismember(vali,greeni));
    redi = find(ismember(vali,redi));
    [irr gssi] = sort(nanmean(ddff(greeni,ninei)'), 'descend');
    [irr rssi] = sort(nanmean(ddff(redi,ninei)'), 'descend');

    imagesc(example_session_r_good_ax, ddff(redi(rssi),ninei),[0 .5]);
    colormap(example_session_r_good_ax, colormap_human);
    imagesc(example_session_g_good_ax, ddff(greeni(gssi),ninei),[0 .5]);
    colormap(example_session_g_good_ax, colormap_human);
    title(example_session_r_good_ax, 'ops+, max: 0.5');
    title(example_session_g_good_ax, 'ops-');
    
    % -- finalize axes
    for a=1:length(ax)
        set(ax(a), 'TickDir', 'out', 'FontSize', 15);
    end
    
    
function rep_drift_dat = get_representational_drift_data(all_dat, settings, force_redo);
    rep_drift_dat_fname = [settings.base_dir filesep 'rep_drift_dat.mat'];

    global rep_drift_dat;
    ld = load(rep_drift_dat_fname);
    rep_drift_dat = ld.rep_drift_dat;

    % quickly-computed things that you may wanat to alter:
    for a=1:length(rep_drift_dat.ani)
        d_dat = rep_drift_dat.ani(a);
        mug = [] ; 
        mur = [] ; 
        cmg = d_dat.cmg;
        cmr = d_dat.cmr;

        for i=1:length(cmg); 
            if(isempty(cmg{i})); continue ; end; 
            tmg = cmg{i}; 
            tmr = cmr{i}; 
    
            goodi = find(nansum(tmr) > 0 & nansum(tmg) > 0); % eliminate weird trials

            muval_r = nan*zeros(1,length(goodi));
            muval_g = nan*zeros(1,length(goodi));
            for t=5:length(goodi)-5
                muval_r(t) = nanmean(reshape(tmr(goodi(t-4:t+4), goodi(t-4:t+4)),[],1));
                muval_g(t) = nanmean(reshape(tmg(goodi(t-4:t+4), goodi(t-4:t+4)),[],1));
            end

            mug(i) = nanmean(muval_g);
            mur(i)=nanmean(muval_r);
        end

        rep_drift_dat.all_mug(a,1:length(mug)) = mug;
        rep_drift_dat.all_mur(a,1:length(mur)) = mur;
    end

