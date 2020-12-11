function [rate, sRate] = spike_rate(  mu, pickChannel, time, window )
% compute mean firing rate across trials
%
% arguments:
% ----------
%   mu              -   cell array of size (numChannel, numTrial) which contains spike
%                       times
%   pickChannel     -   channel number for which to perform analysis
%   time            -   time bins for rate computation
%   window          -   window for counting spikes in each time bin and on
%                       each trial (numTrial or 1, numTimeBin, 2), on the last
%                       dimension [t_beg t_end]
%
% outputs:
% --------
%   rate            -   firing rate
%   sRate           -   firing rate jackknife errorbars

    numTrial = size( mu,2 );
    if (size(window,1)==1)
        window = repmat(window, numTrial, 1);
    end;
    
    count = zeros(numTrial,length(time));
    
    for iTime = 1:length(time)
        for iTrial = 1:numTrial
            spikes = mu{pickChannel,iTrial, : };
            count(iTrial,iTime) = sum( spikes>=window(iTrial,iTime,1) & spikes<window(iTrial,iTime,2) );
        end;
    end;
    
    count = count./squeeze((window(:,:,2)-window(:,:,1)));
    
    jackRate = jackknife(@mean, count);
    rate = mean(jackRate,1);
    sRate = std(jackRate,1)*sqrt(numTrial-1);

end

