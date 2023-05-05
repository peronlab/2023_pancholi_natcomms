function [all_dat settings] = get_all_data_pstim 
    % CHANGE THE LINE BELOW TO WHERE THE DATA LIVES
    settings.base_dir = '~/src/plab/publishedCode/2023_pancholi_pstim/data';
    if (~exist(settings.base_dir, 'dir'))
        disp(['You should change settings.base_dir in get_all_data_pstim - specified directory does not exist: [' settings.base_dir ']' ]);
        return;
    end

    settings.anims = {'an293157','an293167', 'an293010', 'an295578', 'an295678' 'an010852', 'an295137','an295331','an293156','an294755', 'an010580'};
    settings.whisker_idx = [2 2 2 2 2 2 nan nan nan nan nan]; % which whisker to use for touch analyses (C1-C3) ; nan = did not do bumpin

    settings.colors.ops_positive = [1 0 0];
    settings.colors.ops_negative = [0 0.5 1];

    settings.colors.w2_t = [0.5 0.25 1];
    settings.colors.w1_t = [1 0.25 1];
    settings.colors.w3_t = [1 0.5 .75];

    settings.colors.wh = [0.36 0.78 0.53];

    settings.learni = 1:6; % learned
    settings.nlearni = 7:11; % animals that did not lern

    all_dat_fname = [settings.base_dir filesep 'all_quick_dat.mat'];

    global all_dat;
    if (isempty(all_dat))
        disp(['Attempting to load ' all_dat_fname]);
        ld = load(all_dat_fname);
        all_dat = ld.all_dat;
    end
    
    all_dat.settings = settings;
