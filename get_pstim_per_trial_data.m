% This will gather per trial data, meaning
%   - trial type
%   - for each cell, if cell responded
%   - amplitude of said response (pre-post delta)
% Run this within the session_neuropilone directory
function pt_dat = get_pstim_per_trial_data(force_redo)
    if (nargin < 1) ; force_redo = 0 ;end

    fname = [pwd filesep 'pstim_per_trial_summary.mat'];

    global pt_dat;
    if (isempty(pt_dat))
        load(fname);
    end

