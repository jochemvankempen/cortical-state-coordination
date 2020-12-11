function [ emissionSeq, timeBins ] = getEmissionSeq( mu, binSize, timeEpoch )
% [ emissionSeq, timeBins ] = getEmissionSeq( mu, binSize, timeEpoch )
%
% generates emission sequence from the spike-times data
% emission sequence contains the spike count for each
% neuron in the population at each time bin
%
% arguments:
% ---------
%   mu          -   cell array of size (numChannel, numTrial) that contains
%                   spike times
%   binSize     -   bin size used to bin spiking data for generating
%                   emission sequence
%   timeEpoch   -   time epoch that should be used for the analysis, it can
%                   be either an interval [tStart, tEnd], or matrix of size
%                   (numTrial, 2) that contains tStaet and tEnd for each
%                   trial
%
% outputs:
% -------
%   emissionSeq     -   cell array of length numTrial, each element
%                       contains an array of size (numChannel, numBin),
%                       where numChannel is the number of neurons and
%                       numBin is the number of time bins on that trial
%   timeBins        -   cell array of length numTrial, that contains time
%                       bins edges at which spike counts were obtained on
%                       each trial to generate the emission sequence

    numChannel = size(mu,1);
    numTrial = size(mu,2);
    
    if size(timeEpoch,1) == 1
        timeEpoch = repmat(timeEpoch,numTrial,1);
    end

    % generate emission sequence from spikes
    emissionSeq = cell(numTrial,1);
    timeBins = cell(numTrial,1);
    for iTrial = 1:numTrial
        timeBins{iTrial} = timeEpoch(iTrial,1):binSize:timeEpoch(iTrial,2);
        numBin = length(timeBins{iTrial})-1;
        emissionSeq{iTrial} = zeros(numChannel, numBin);
    end;

    for iTrial = 1:numTrial
        for iChann = 1:numChannel
            spikes = mu{iChann,iTrial, : };
            spikes = spikes( spikes>timeEpoch(iTrial,1) & spikes<=timeEpoch(iTrial,2) );
            if (~isempty(spikes))
                n = histc( spikes, timeBins{iTrial}  );
                emissionSeq{iTrial}(iChann,:) = n(1:end-1);
            end;
        end;
        timeBins{iTrial} = timeBins{iTrial}(1:end-1) + 0.5*binSize;
    end;
    
end

