function beh_dat = get_pstim_beh_data
    [all_dat settings] = get_all_data_pstim;

    fname = [settings.base_dir filesep 'pstim_behavior_summary.mat'];
    global beh_dat;
    if (isempty(beh_dat))
        load(fname);
    end
