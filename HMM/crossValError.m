function [ cvError, looError, veError, totalCount ] = crossValError(testSeq, estTrans, estEmis, estPi0, chann, timeWin)
% for each neuron computes the cross-validation error using estimated
% parameters of the HMM, and the activities of numChannel-1 neurons to
% decode the sequence of states on the test trials
%
% arguments:
% ---------
%   testSeq                     -   emission sequence on the test trials
%   estTrans, estEmis, estPi0   -   parameters of the HMM fitted on the
%                                   train trials
%   timeWin                     -   array that contains windows in which to perform
%                                   cross-validation, in units of the number of bins
%
% outputs:
% -------
%   cvError     - array of size (numChannel,numWin) with the average
%                 cross-validation error for each neuron and for each
%                 window size

    
    numWin = length(timeWin);

    numTrial = length(testSeq);
    numChannel = size(estEmis,2);
    cvError = zeros(numChannel,numWin);
    looError = zeros(numChannel,numWin);
    veError = zeros(numChannel,numWin);
    totalCount = zeros(numChannel,1);

    % determine indeces of window sizes for each there is not enough long
    % data for cross-validation, to set them to NaN
    maxNumBin = max( cellfun(@(x) size(x,2), testSeq) );
    indNan = false(size(timeWin));
    indNan(timeWin>maxNumBin) = true;
   
    % leave out neurons one by one, also do not include nearby channels
    for iChann = 1:numChannel
        cve = zeros(numWin,1); % cve for each window size
        cve(indNan) = NaN;
        loo = zeros(numWin,1); % cve for each window size
        loo(indNan) = NaN;
        vee = zeros(numWin,1); % vee for each window size
        vee(indNan) = NaN;
        for iTrial = 1:numTrial
            ind = true(numChannel,1); % exclude units from the same and nearby electrodes as the current unit
            ind(chann==chann(iChann) | chann==chann(iChann)+1 | chann==chann(iChann)-1) = false;
            
            indLoo = true(numChannel,1); % exclude units from the same electrode as the current unit
            indLoo(chann==chann(iChann)) = false;
            
            emis = estEmis(:, ind);
            seq = testSeq{iTrial}(ind, :);
            % decode states based on the activity of remaining neurons
            state = hmmviterbiPoisson(seq,estTrans,emis, estPi0);
            % decode states based on the activity of remaining neurons
            stateLoo = hmmviterbiPoisson(testSeq{iTrial}(indLoo, :),estTrans,estEmis(:, indLoo), estPi0);
            % decode states based on the activity of all neurons
            stateFull = hmmviterbiPoisson(testSeq{iTrial},estTrans,estEmis, estPi0);
            % predict firing rate
            predCount = arrayfun(@(x) estEmis(x, iChann), state);
            predCountLoo = arrayfun(@(x) estEmis(x, iChann), stateLoo);
            predCountFull = arrayfun(@(x) estEmis(x, iChann), stateFull);
            diffCount = predCount - testSeq{iTrial}(iChann,:); % difference between predicted and actual spike count
            diffCountLoo = predCountLoo - testSeq{iTrial}(iChann,:); % difference between predicted and actual spike count
            diffCountFull = predCountFull - testSeq{iTrial}(iChann,:); % difference between predicted and actual spike count
            totalNumBin  = length(diffCount);
            for iWin = 1:numWin
                numBin = timeWin(iWin); % number of bins over which to compute error
                for iBin = 1:floor(totalNumBin/numBin)
                    cve(iWin) = cve(iWin) + (sum(diffCount(1+(iBin-1)*numBin:iBin*numBin) ))^2;
                    loo(iWin) = loo(iWin) + (sum(diffCountLoo(1+(iBin-1)*numBin:iBin*numBin) ))^2;
                    vee(iWin) = vee(iWin) + (sum(diffCountFull(1+(iBin-1)*numBin:iBin*numBin) ))^2;
                end;
            end;
        end;
        totalCount(iChann) = sum( cellfun(@(x) sum(x(iChann,:)), testSeq) ); % total spike count
        cvError(iChann,:) = cve;    %/totCount;
        looError(iChann,:) = loo;   %/totCount;
        veError(iChann,:) = vee; %/totCount;
    end;
    
    

end