function figure_6_interday_stability
    
    % --- get data
    [all_dat settings] = get_all_data_pstim;
    to_dat = get_pstim_turnover_data;
    to_dat_pstim = to_dat.to_dat_pstim_presp;

    % --- setup figure
    fsize = 10;

    % example data get
    ex_ani = find(strcmp(all_dat.settings.anims, 'an293157'));
    quick_dat = all_dat.quick_dat{ex_ani};
    p_response_thresh = 0.25;
    p_response_streak_thresh = 0.1;
    [Mr dMr] = plot_example_turnover(nan, quick_dat, p_response_thresh, p_response_streak_thresh, quick_dat.types.red, 'ops+, P-resp',fsize );
    [Mg dMg]= plot_example_turnover(nan, quick_dat, p_response_thresh, p_response_streak_thresh, quick_dat.types.green, 'ops-', fsize);


    % --- newer plots -- figure 3, o+/o- pairwise correlation across days
    fh3 = figure ('Position',[0 0 1200 600]);
    example_mouse_corr(1) = subplot('Position',[.05 .05 .15 .3]); 
    example_mouse_corr(2) = subplot('Position',[.05 .55 .15 .3]);

    cross_mouse_corr(1) = subplot('Position',[.25 .05 .15 .3]); 
    cross_mouse_corr(2) = subplot('Position',[.25 .55 .15 .3]);
    cross_mouse_corr(3) = subplot('Position',[.45 .05 .15 .3]); 
    cross_mouse_corr(4) = subplot('Position',[.45 .55 .15 .3]);
    cross_mouse_corr(5) = subplot('Position',[.65 .65 .15 .3]);
    cross_mouse_corr(6) = subplot('Position',[.825 .65 .15 .3]);

    cross_mouse_stats(1) = subplot('Position',[.65 .1 .3 .2]); 
    cross_mouse_stats(2) = subplot('Position',[.65 .4 .3 .2]);

    for i=1:2 
        set(example_mouse_corr(i), 'TickDir', 'out','FontSize', 15);
        set(cross_mouse_stats(i), 'TickDir', 'out','FontSize', 15);
    end
    for i=1:4
        set(cross_mouse_corr(i), 'TickDir', 'out','FontSize', 15);
    end

    ex_ani = 1;
    nn = isnan(to_dat.to_dat_pstim_presp(ex_ani).day_corr_good_cells_green); % eliminate days that are nan - bad
    vali = find(sum(nn) < 5);
    imagesc(example_mouse_corr(1), to_dat.to_dat_pstim_presp(ex_ani).day_corr_good_cells_red(vali,vali),[ 0 1]);
    colormap(example_mouse_corr(1), colormap_human);
    title(example_mouse_corr(1), 'Example mouse, ops+');
    xlabel(example_mouse_corr(1), 'Day');
    ylabel(example_mouse_corr(1), 'Day');

    imagesc(example_mouse_corr(2), to_dat.to_dat_pstim_presp(ex_ani).day_corr_good_cells_green(vali,vali),[ 0 1]);
    colormap(example_mouse_corr(2), colormap_human);
    title(example_mouse_corr(2), 'Example mouse, ops-');
    xlabel(example_mouse_corr(2), 'Day');
    ylabel(example_mouse_corr(2), 'Day');
        
    learn_ops_pos = plot_corr_across_days (cross_mouse_corr(1), to_dat_pstim(settings.learni), all_dat, settings, dMr, 'day_corr_zsc_dff_red', fsize);
%    learn_ops_pos = plot_corr_across_days (cross_mouse_corr(5), to_dat_pstim(settings.learni), all_dat, settings, dMr, 'day_corr_zsc_dff_red', fsize, 8);
    learn_ops_neg = plot_corr_across_days (cross_mouse_corr(2), to_dat_pstim(settings.learni), all_dat, settings, dMr, 'day_corr_zsc_dff_green', fsize);
%    learn_ops_neg = plot_corr_across_days (cross_mouse_corr(6), to_dat_pstim(settings.learni), all_dat, settings, dMr, 'day_corr_zsc_dff_green', fsize, 8);
    nlearn_ops_pos = plot_corr_across_days (cross_mouse_corr(3), to_dat_pstim(settings.nlearni), all_dat, settings, dMr, 'day_corr_zsc_dff_red', fsize);
    nlearn_ops_neg = plot_corr_across_days (cross_mouse_corr(4), to_dat_pstim(settings.nlearni), all_dat, settings, dMr, 'day_corr_zsc_dff_green', fsize);

    for x=1:4 ; axis(cross_mouse_corr(x), [0.5 22.5 0.5 22.5]); end

    learn_col = [1 0.25 1];
    nlearn_col = [1 1 1]*0.25;

    hold(cross_mouse_stats(1), 'on');
    di =  1:3; % we have 8 days for nonlearners but this will underestimate 
%    di =  1:8; % we have 8 days for nonlearners but this will underestimate 
    axes(cross_mouse_stats(1));
    plot_error_poly (di, nanmean(learn_ops_pos(:,di)), nanstd(learn_ops_pos(:,di))/sqrt(length(settings.learni)), learn_col, [1 0.75 1]);    
    plot_error_poly (di, nanmean(nlearn_ops_pos(:,di)), nanstd(nlearn_ops_pos(:,di))/sqrt(length(settings.nlearni)), nlearn_col, [1 1 1]*0.75);    
    axis(cross_mouse_stats(1), [0 di(end)+2 0 8]);
    xlabel(cross_mouse_stats(1), 'Training day');
    ylabel(cross_mouse_stats(1), 'Days until R2 < 0.5');
    rnd = (rand(1,length(settings.learni))-0.5)*.1;
    plot(4+0*(1:length(settings.learni))+rnd,  nanmean(learn_ops_pos(:,di)'), 'o', 'Color', 'None', 'MarkerFaceColor', learn_col, 'MarkerSize', 5);
    plot(4,  nanmean(nanmean(learn_ops_pos(:,di)')), 'o', 'Color', learn_col ,'MarkerFaceColor', 'None', 'MarkerSize', 10, 'LineWidth', 2);
    rnd = (rand(1,length(settings.nlearni))-0.5)*.1;
    plot(4.5+0*(1:length(settings.nlearni))+rnd,  nanmean(nlearn_ops_pos(:,di)'), 'o', 'Color', 'None', 'MarkerFaceColor', nlearn_col, 'MarkerSize', 5);
    plot(4.5,  nanmean(nanmean(nlearn_ops_pos(:,di)')), 'o', 'Color', nlearn_col ,'MarkerFaceColor', 'None', 'MarkerSize', 10, 'LineWidth', 2);
    [h p] = ttest2(nanmean(nlearn_ops_pos(:,di)'), nanmean(learn_ops_pos(:,di)'));
    title(cross_mouse_stats(1), sprintf('Ops+ p=%0.3f', p));
    disp(sprintf('Ops+ p=%0.3f; lrn mu/sd: %0.3f/%0.3f nlrn mu/sd: %0.3f/%0.3f', p, nanmean(nanmean(learn_ops_pos(:,di))), nanstd(nanmean(learn_ops_pos(:,di))), ...
                                                                                  nanmean(nanmean(nlearn_ops_pos(:,di))), nanstd(nanmean(nlearn_ops_pos(:,di)))));
    beh_corr_plot (cross_mouse_corr(5), all_dat, [nanmean(learn_ops_pos(:,di)') nanmean(nlearn_ops_pos(:,di)')], 'o+');


    hold(cross_mouse_stats(2), 'on');
    axes(cross_mouse_stats(2));
    plot_error_poly (di, nanmean(learn_ops_neg(:,di)), nanstd(learn_ops_neg(:,di))/sqrt(length(settings.learni)), learn_col, [1 0.75 1]);    
    plot_error_poly (di, nanmean(nlearn_ops_neg(:,di)), nanstd(nlearn_ops_neg(:,di))/sqrt(length(settings.nlearni)), nlearn_col, [1 1 1]*0.75);    
    axis(cross_mouse_stats(2), [0 di(end)+2 0 8]);
    rnd = (rand(1,length(settings.learni))-0.5)*.1;
    plot(4+0*(1:length(settings.learni))+rnd,  nanmean(learn_ops_neg(:,di)'), 'o', 'Color', 'None', 'MarkerFaceColor', learn_col, 'MarkerSize', 5);
    plot(4,  nanmean(nanmean(learn_ops_neg(:,di)')), 'o', 'Color', learn_col ,'MarkerFaceColor', 'None', 'MarkerSize', 10, 'LineWidth', 2);
    rnd = (rand(1,length(settings.nlearni))-0.5)*.1;
    plot(4.5+0*(1:length(settings.nlearni))+rnd,  nanmean(nlearn_ops_neg(:,di)'), 'o', 'Color', 'None', 'MarkerFaceColor', nlearn_col, 'MarkerSize', 5);
    plot(4.5,  nanmean(nanmean(nlearn_ops_neg(:,di)')), 'o', 'Color', nlearn_col ,'MarkerFaceColor', 'None', 'MarkerSize', 10, 'LineWidth', 2);
    [h p] = ttest2(nanmean(nlearn_ops_neg(:,di)'), nanmean(learn_ops_neg(:,di)'));
    disp(sprintf('Ops- p=%0.3f; lrn mu/sd: %0.3f/%0.3f nlrn mu/sd: %0.3f/%0.3f', p, nanmean(nanmean(learn_ops_neg(:,di))), nanstd(nanmean(learn_ops_neg(:,di))), ...
                                                                                  nanmean(nanmean(nlearn_ops_neg(:,di))), nanstd(nanmean(nlearn_ops_neg(:,di)))));
    beh_corr_plot (cross_mouse_corr(6), all_dat, [nanmean(learn_ops_neg(:,di)') nanmean(nlearn_ops_neg(:,di)')], 'o-');

    di1=1:3;
    di2=8:10;
    [h p] = ttest2(nanmean(learn_ops_pos(:,di1)), nanmean(learn_ops_pos(:,di2)));
    disp(sprintf('Ops+ d1 v 10 p=%0.3f; d1 lrn mu/sd: %0.3f/%0.3f d10 lrn mu/sd: %0.3f/%0.3f', p, nanmean(nanmean(learn_ops_pos(:,di1))), nanstd(nanmean(learn_ops_pos(:,di1))), ...
                                                                                  nanmean(nanmean(learn_ops_pos(:,di2))), nanstd(nanmean(learn_ops_pos(:,di2)))));
    [h p] = ttest2(nanmean(learn_ops_neg(:,di1)), nanmean(learn_ops_neg(:,di2)));
    disp(sprintf('Ops- d1 v 10 p=%0.3f; d1 lrn mu/sd: %0.3f/%0.3f d10 lrn mu/sd: %0.3f/%0.3f', p, nanmean(nanmean(learn_ops_neg(:,di1))), nanstd(nanmean(learn_ops_neg(:,di1))), ...
                                                                                  nanmean(nanmean(learn_ops_neg(:,di2))), nanstd(nanmean(learn_ops_neg(:,di2)))));


function plot_corr_trend (ax, to_dat_pstim, settings, trend_shown)
    day_centers = [2 8  10 ];
    hold(ax, 'on');
    ms = 10;

    col_ops = settings.colors.ops_positive;
    col_nops =  settings.colors.ops_negative;

    vmat_r = nan*zeros(length(day_centers), length(to_dat_pstim));
    vmat_g = nan*zeros(length(day_centers), length(to_dat_pstim));

    for c=1:length(day_centers)
        vals_g = [];
        vals_r = [];
        for a=1:length(to_dat_pstim)
            ci = day_centers(c)+[-1 0 1];
            if (strcmp(trend_shown, 'zsc_dff'))
                vals_r(a) = nanmean(to_dat_pstim(a).day_corr_zsc_dff_min_days_to_thresh_red(ci));
                vals_g(a) = nanmean(to_dat_pstim(a).day_corr_zsc_dff_min_days_to_thresh_green(ci));
            elseif (strcmp(trend_shown, 'p_resp'))
                vals_g(a) = nanmean(to_dat_pstim(a).day_corr_p_resp_min_days_to_thresh_green(ci));
                vals_r(a) = nanmean(to_dat_pstim(a).day_corr_p_resp_min_days_to_thresh_red(ci));
            end
        end
        plot(ax, -0.1+c*[1 1], nanmean(vals_g)+nanstd(vals_g)/sqrt(length(to_dat_pstim))*[-1 1], '-', 'Color', col_nops, 'LineWidth', 2);
        plot(ax, 0.1+c*[1 1], nanmean(vals_r)+nanstd(vals_r)/sqrt(length(to_dat_pstim))*[-1 1], '-', 'Color', col_ops, 'LineWidth', 2);
 
        plot(ax, c-0.1, nanmean(vals_g), 'o', 'Color', [1 1 1]*0.8, 'MarkerFaceColor', col_nops, 'MarkerSize', ms);
        plot(ax, c+0.1, nanmean(vals_r), 'o', 'Color', [1 1 1]*0.8, 'MarkerFaceColor', col_ops, 'MarkerSize', ms);

        vmat_r(c,:) = vals_r;
        vmat_g(c,:) = vals_g;
    end
    axis(ax,[0 length(day_centers)+1 0 10]);
    set(ax, 'TickDir', 'out', 'XTick', [1 2 3], 'XTickLabels', {'d1-3', '7-9', '14-16'});
    ylabel(ax, 'Days until R^2 < 0.5');

    [h p] = ttest(vmat_r(1,:), vmat_r(3,:));
    disp(sprintf('red early v. late days to R2 0.5 %0.1f +/- %0.1f ; %0.1f +/- %0.1f pval: %0.3f', nanmean(vmat_r(1,:)),nanstd(vmat_r(1,:)),  nanmean(vmat_r(3,:)), nanstd(vmat_r(3,:)), p));

    [h p] = ttest(vmat_g(1,:), vmat_g(3,:));
    disp(sprintf('green early v. late days to R2 0.5 %0.1f +/- %0.1f ; %0.1f +/- %0.1f pval: %0.3f', nanmean(vmat_g(1,:)),nanstd(vmat_g(1,:)),  nanmean(vmat_g(3,:)), nanstd(vmat_g(3,:)), p));

    [h p] = ttest(vmat_g(1,:), vmat_r(1,:));
    disp(sprintf('early g v r: %0.3f', p));
    [h p] = ttest(vmat_g(3,:), vmat_r(3,:));
    disp(sprintf('late g v r: %0.3f', p));


function mu_per_ani = plot_corr_across_days (ax, to_dat_pstim, all_dat, settings, M_ex_ani, mat_fld_name, fsize, m_size)
    thresh_corr = 0.5;

    n_days = [to_dat_pstim.num_days];
    if (nargin < 8 || isempty(m_size))
        m_size = min(n_days);
    end
m_size = 10;
    col_ops = settings.colors.ops_positive;
    col_nops =  settings.colors.ops_negative;

    % correlation matrix generate (Avg across days)
    n_mu_days = 10; 
    day_corr_all = zeros(m_size, m_size, length(to_dat_pstim));
    mu_per_ani = nan*zeros(length(to_dat_pstim),n_mu_days);
    for a=1:length(to_dat_pstim)
        n_days_ani = min(m_size, to_dat_pstim(a).num_days);
        day_corr_all(1:n_days_ani,1:n_days_ani,a) = to_dat_pstim(a).(mat_fld_name)(1:n_days_ani, 1:n_days_ani);
        for d=1:n_days_ani
            i1 = min(find( to_dat_pstim(a).(mat_fld_name)(d:end, d) < sqrt(thresh_corr)));
            if (isempty(i1)) ; i1 = nan ; end
            mu_per_ani(a,d) = i1-1;
%        mu_per_ani(a,:) = nanmean( to_dat_pstim(a).(mat_fld_name)(1:n_mu_days, 1:n_mu_days));
        end
    end
    day_corr = nanmean(day_corr_all,3);

    % track down the line of 0.5
    min_thresh_d1 = nan*zeros(1,size(day_corr,1));
    for d=1:size(day_corr,1)
        candi = find(day_corr(:,d) > thresh_corr);
        mi = candi(min(find(candi > d)));
        if (length(mi) == 1)
            min_thresh_d1(d) = mi+d;
        end
    end

    axes(ax(1));
    imagesc(day_corr, [0 1]);
    colormap(colormap_human);
    hold(ax(1),'on');
    plot(1:size(day_corr,1), min_thresh_d1, '-','Color',[1 1 1]);
    title(ax(1), strrep(mat_fld_name,'_','-'));
    set(ax(1), 'TickDir','out','FontSize',fsize, 'XTick', [1 size(day_corr,1)], 'YTick', [1 size(day_corr,2)]);

    % example pair days on axes 2 and 3
    if (length(ax) == 4) 
        ex_ani = 1;
        
        day_pairs = {[1 2], [1 5], [1 10]};
        if (strfind(mat_fld_name, 'red'))
            colr = col_ops;
        else
            colr = col_nops;
        end

        for d=1:length(day_pairs)
            hold(ax(d+1), 'on');
            plot(ax(d+1), [0 1], [0  1], 'k:');
            plot(ax(d+1), M_ex_ani(day_pairs{d}(1),:)/max(M_ex_ani(day_pairs{d}(1),:)), M_ex_ani(day_pairs{d}(2),:)/max(M_ex_ani(day_pairs{d}(2),:)), ...
                 'o', 'MarkerFaceColor', colr, 'MarkerEdgeColor', [1 1 1]*0.2, 'MarkerSize', 8);
            axis(ax(d+1),[-0.1 1.1 -0.1 1.1]);
            set(ax(d+1), 'TickDir','out','FontSize',fsize, 'XTick',[0 1], 'YTick', [0 1]);
            xlabel(ax(d+1), sprintf('Day %d, normd. zdff', day_pairs{d}(1)));
            ylabel(ax(d+1), sprintf('Day %d, normd. zdff', day_pairs{d}(2)));
            title(ax(d+1), [ 'R: ' num2str(nancorr( M_ex_ani(day_pairs{d}(1),:), M_ex_ani(day_pairs{d}(2),:)))]);
        end
    end


function [retM retMdff] = plot_example_turnover(ax, quick_dat, p_response_thresh, p_response_streak_thresh, restricti, tstr, fsize)

    % Turnover ... last days are sometimes weird due to lesion, so drop
    vali = find(nanmax(quick_dat.probability_response{6}(1:end-1,:)) > p_response_thresh);
    vali = intersect(vali,restricti);
    M = quick_dat.probability_response{6}(1:end-1,vali);
    Mdff = quick_dat.dff_sd_normed_response{6}(1:end-1,vali);
    [maxval peaki] = nanmax(M);

    % find first and last ...
    firsti = [];
    lasti = [];
    for i=1:size(M,2); 
        mf = min(find(M(:,i) > p_response_streak_thresh)) ;
        ml = max(find(M(:,i) > p_response_streak_thresh)) ;
        if (length(mf) == 1)
            firsti(i) = mf;
        end
        if (length(ml) == 1)
            lasti(i) = ml;
        end
    end
    peaki=firsti;
    [irr sorti] = sort(peaki);

    % subsort by LAST day
    R = quick_dat.probability_response{6}(1:end-1,vali);
    R(find(R > 0.1)) = 1;
    R(find(R <= 0.1)) = 0;

    up = unique(peaki);
    nsorti = [];
    for p=1:length(up)
        ii = find(peaki == up(p));
        [irr subsorti] = sort(lasti(ii));
        nsorti = [nsorti ii(subsorti)];
    end
    sorti = nsorti;
    
    % plot
    if (~isnan(ax))
        axes(ax);
        imagesc(M(:,sorti)', [0 1]); 
        colormap (ax,colormap_human);
        title(ax, tstr);
        set(ax, 'TickDir','out','FontSize',fsize, 'XTick', [1 size(M,1)], 'YTick', [1 size(M,2)]);
    end

    % return 
    retM = M(:,sorti);
    retMdff = Mdff(:,sorti);
