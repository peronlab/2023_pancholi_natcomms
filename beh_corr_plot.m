
function beh_corr_plot (ax, all_dat, val_all_vec, tstr)
    if (nargin < 4) ; tstr = ''; end

    learni = all_dat.settings.learni;
    nlearni = all_dat.settings.nlearni;
    
    o_di_beh = 6:8; % who will we use for behavior (will use LATEST day available of this range)
   % o_di_beh = 9:12; % who will we use for behavior (will use LATEST day available of this range)
o_di_beh =1:8;
    bp = get_pstim_beh_data;
    for a=1:length(bp)
        bp(a).frac_correct = bp(a).frac_correct(find(~isnan(bp(a).frac_correct)));
        bp(a).frac_correct_9pulse = bp(a).frac_correct_9pulse(find(~isnan(bp(a).frac_correct_9pulse)));
    %    disp(sprintf('%s : %d', all_dat.settings.anims{a}, length(bp(a).frac_correct)));


    end

    beh_all_vec = nan*zeros(1,length(learni)+length(nlearni));
    n_days_thresh = nan*zeros(1,length(learni)+length(nlearni));
    alli = [learni nlearni];
    for a=1:length(alli)
        ai = alli(a);
        di_beh = o_di_beh(find(o_di_beh <= length((bp(ai).frac_correct_9pulse))));
    
        if (length(di_beh) == 1)
            di_beh = [1 1]*di_beh;
        elseif (length(di_beh) > 1)
            di_beh = [1 1]*di_beh(end);
        end

        if (length(di_beh) > 0)
            beh_all_vec(ai) = nanmean(bp(ai).frac_correct_9pulse(di_beh));
            beh_all_vec(ai) = nanmax(bp(ai).peak_dprime(di_beh)/5);
            %beh_all_vec(ai) = nanmean(bp(ai).frac_correct_9pulse);

            if (length(find(bp(ai).frac_correct_9pulse(3:end) > .75)) > 0)
                n_days_thresh(ai) = 2+min(find(bp(ai).frac_correct_9pulse(3:end) > .75));
            end
        end
    end
    n_days_thresh(nlearni) = max(n_days_thresh(learni)+1);
%    beh_all_vec = n_days_thresh/max(n_days_thresh)

    [r p] = nancorr(val_all_vec, beh_all_vec);
    hold (ax, 'on');
    plot(ax, val_all_vec(learni),  beh_all_vec(learni),'o', 'MarkerEdgeColor', [1 1 1]*.75, 'MarkerFaceColor', [1 0.25 0], 'MarkerSize', 8);
    plot(ax,val_all_vec(nlearni), beh_all_vec(nlearni),  'o', 'MarkerEdgeColor', [1 1 1]*.75, 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    title(ax, sprintf('%s r=%0.3f p=%0.3f', tstr, r, p));
    aa = axis(ax);
    axis(ax, [aa(1) aa(2) 0.35 1.05 ]);

