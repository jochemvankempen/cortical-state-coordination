function cfg = set_cfg(recInfo)
% cfg = set_cfg(recInfo)
% 
% Parameters
% ----------
% recInfo : table
%     table with single row (one recording)
%
% Returns
% -------
% cfg : struct
%     struct with info on analysis settings, channel exclusion etc.
% 
% 
% Notes
% -----
% cfg.snr : signal-noise ratio threshold (channels below this are excluded)
% cfg.channelSelectionCriteria : criteria: criteria for channel in/exclusion
% 
%     - 'snr'       : select channels based on SNR value, the threshold for inclusion is set in variable cfg.snr
%     - 'depth'     : select channels based on the recording depth
%     - 'stability  : Based on recording stability
% 
% 

% task specific settings
switch recInfo.Task
    case 'gratc'
        cfg.snr = 3; % signal-noise ratio threshold (channels below this are excluded)
        cfg.channelSelectionCriteria = {'snr','depth','stability'}; 
        
        cfg.ignoreDirection = true;% boolean, collapse across stimulus direction

        switch recInfo.Subject{1}
            %%% set different time windows for different subjects when aligning to
            %%% fixation. This takes into account the visual transient after
            %%% fixation onset. Because the time between fixation and stimulus
            %%% onset is shorter for W compared to T, we take a slightly
            %%% shorter time. Luckily, W also has a shorter visual
            %%% transient after fix onset. 
            %%% Analyses will only run until next event
            case 'J'
                cfg.fixAlignTimeWindow = [0.150 0.600];
            case 'T'
                cfg.fixAlignTimeWindow = [0.200 0.600]; %
            case 'W'
                cfg.fixAlignTimeWindow = [0.100 0.400];
            otherwise
                error('List path not defined for %s', Subject)
                cfg.fixAlignTimeWindow = [NaN];

        end
        
end

