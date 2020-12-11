function estEmis = estimateRateState(seq, states, timeBins, numState, binsize, rateUnits)
% estEmis = estimateRateState(trainSeq, trainStates, numState)
%
%   estimate firing rate of each unit in each of the states
%
% arguments:
% ----------
%   seq         -   emission sequence, cell array of size numTrial, each cell contains the spike
%                   count for each neuron in each time bin on that trial
%   states      -   sequence of states on each of the trials decoded from the HMM fit
%   timeBins    -   sequence of corresponding time bins on each trial, cell array
%   numState    -   number of states in the HMM
%   rateUnits   -   string 'Hz' or 'prob', specifying which units should be used for estEmis
%
% outputs:
% --------
%   estEmis     -   estimated emission matrix, array of size (numState, numChannel)


    numTrial = length(states);
    numChannel = size(seq{1},1);
    
    estEmis = zeros(numState, numChannel);
    
    [ timeEpochState, durationState, ~ ] = extractUpDownTimeEpoch( states, timeBins, numState, binsize, true);
    

    for iState = 1:numState
        norm = 0;
        for iTrial = 1:numTrial
            thisTimeEpochState = timeEpochState{iTrial,iState};
            numEpoch = size(thisTimeEpochState,1);
            ind = false(1, length(timeBins{iTrial}));
            for iEpoch = 1:numEpoch
                ind = ind | ( timeBins{iTrial}>=thisTimeEpochState(iEpoch,1) & timeBins{iTrial}<thisTimeEpochState(iEpoch,2) );
            end;
            estEmis(iState,:) = estEmis(iState,:) + sum(seq{iTrial}(:,ind),2)';
            norm = norm + sum(ind);
        end;
        switch rateUnits
            case 'Hz'
                estEmis(iState,:) = estEmis(iState,:)/sum( durationState(:,iState) );
            case 'prob'
                estEmis(iState,:) = estEmis(iState,:)/norm;
            otherwise
                error(['Unknown specification of emision rate units:' rateUnits ', use Hz or prob']);
        end;
    end;
    
   
end

