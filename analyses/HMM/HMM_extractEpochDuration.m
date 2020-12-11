function [epochDuration, fracTimeState, residTimePDF] = HMM_extractEpochDuration(timeEpochState, numState)


% probability density function
binsize_PDF = 0.05;
nbins_PDF = 25;
% bincentres_PDF = (0+binsize_PDF/2):binsize_PDF:(binsize_PDF*nbins_PDF - binsize_PDF/2);
edges = 0:binsize_PDF:binsize_PDF*nbins_PDF;


fracTimeState = NaN(1,numState);
epochDuration = NaN(1,numState);
residTimePDF = NaN(nbins_PDF+1,2,numState);
for istate = 1:numState
    %%% get the duration of each epoch for each state
    timeEpochThisState = vertcat( timeEpochState{:,istate} );
    residTime = timeEpochThisState(:,2) - timeEpochThisState(:,1);
    
    fracTimeState(istate) = sum(residTime);
    
    epochDuration(istate) = nanmean(residTime);
    
    
    %[f, x] = ecdf( residTime );
    %residTimeCDF{istate} = [x, f];
    %[n, x] = hist(residTime, bincentres_PDF);
    [n] = histc(residTime(:)', edges);
    n = n/sum(n);
    residTimePDF(:,:,istate) = [edges', n'];
    
end
fracTimeState = fracTimeState/sum(fracTimeState); % /totalTime
end

