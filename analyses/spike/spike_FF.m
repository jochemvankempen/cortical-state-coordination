function [FF, sFF, meanCount, varCount] = spike_FF(  msu, pickChannel, time, window, numSubWindow)
% compute Fano Factor for specified channel
%
% arguments:
% ----------
%   msu             -   cell array of size (numChannel, numTrial) which contains spike
%                       times
%   pickChannel     -   channel number for which to perform analysis
%   time            -   time bins for the Fano Factor computation
%   window          -   window size bins for the Fano Factor computation
%   numSubWindow    -   array that contains number of subwindows for the FF
%                       computation, for each time bin or for each window size
%
% outputs:
% --------
%   FF              -   Fano Factor
%   sFF             -   Fano Factor jackknife errorbars
%   meanCount       -   mean spike count in each time bin for each window
%                       size (numTimeBin, numWinSize)
%   varCount        -   variance of spike count in each time bin for each
%                       window size (numTimeBin, numWinSize)


    numTrial = size( msu,2 );

    FF = zeros(length(time),length(window));
    sFF = zeros(length(time),length(window));
    meanCount = zeros(length(time),length(window));
    varCount = zeros(length(time),length(window));
    
    spikes = msu(pickChannel,:);
    for iTime = 1:length(time)
       for iWindow = 1:length(window)

            numSW = numSubWindow(iWindow,:);
            win = window(iWindow);
            
            if (length(numSW)==1)
                allCount = zeros(numSW*numTrial, 1);
                for iWin = 1:numSW
                    allCount(1+(iWin-1)*numTrial:iWin*numTrial) = cellfun(@(x) sum( x>=time(iTime)+(iWin-1)*win & x<time(iTime)+iWin*win ), spikes);
                end;
            else
                assert(length(numSW)==numTrial,'Length of array numSW should be equal to the number of trials');
                allCount = zeros(sum(numSW), 1);
                ind = [0 cumsum(numSW)];
                for iTrial = 1:numTrial
                    winArray = time(iTime):win:time(iTime)+(numSW(iTrial)-1)*win;
                    allCount(1+ind(iTrial):ind(iTrial+1)) = arrayfun(@(x) sum(spikes{iTrial}>=x & spikes{iTrial}<x+win), winArray);
                end;
            end;      
            
            jackVar = jackknife(@var, allCount, 1);
            jackMean = jackknife(@mean, allCount);
            jackFF = jackVar./jackMean;
            FF(iTime, iWindow) = nanmean(jackFF,1);
            sFF(iTime, iWindow) = nanstd(jackFF,1).*sqrt( max(sum(~isnan(jackFF),1)-1,0) );
            
            meanCount(iTime, iWindow) = mean(allCount, 1);
            varCount(iTime, iWindow) = var(allCount, 1);
       end;
    end;


end

