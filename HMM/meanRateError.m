function mrError = meanRateError(trainSeq, testSeq, timeWin)
% mrError = meanRateError(trainSeq, testSeq)
%
% for each neuron computes the cross-validation error using prediction
% based on its mean firing rate estimated from the training sequence
%
% arguments:
% ---------
%   trainSeq    -   training sequence from which the mean firing rate is
%                   estimated
%   testSeq     -   test sequence
%   timeWin                     -   array that contains windows in which to perform
%                                   cross-validation, in units of the number of bins
%
% outputs:
% -------
%   mrError     -   array of size numChannel with the average
%                   cross-validation error for each neuron
    
    numWin = length(timeWin);
    numChannel = size(trainSeq{1},1);
    numTrial = length(testSeq);
    
    % determine indeces of window sizes for each there is not enough long
    % data for cross-validation, to set them to NaN
    maxNumBin = max( cellfun(@(x) size(x,2), testSeq) );
    indNan = false(size(timeWin));
    indNan(timeWin>maxNumBin) = true;
    
    % estimate mean firing rate from training sequence
    meanRate = zeros(numChannel, 1);
    for iChann = 1:numChannel
        meanRate(iChann) = sum( cellfun(@(x) sum( x(iChann,:)), trainSeq ) )/sum( cellfun(@(x) size(x,2), trainSeq ) );
    end;
    
    % compute cross-validation error
    mrError = zeros(numChannel,numWin);
    for iChann = 1:numChannel      
        cve = zeros(numWin,1); % cve for each window size
        cve(indNan) = NaN;
        for iTrial = 1:numTrial
            diffCount = meanRate(iChann) - testSeq{iTrial}(iChann,:); % difference between predicted and actual spike count
            totalNumBin  = length(diffCount);
            for iWin = 1:numWin
                numBin = timeWin(iWin); % number of bins over which to compute error
                for iBin = 1:floor(totalNumBin/numBin)
                    cve(iWin) = cve(iWin) + (sum(diffCount(1+(iBin-1)*numBin:iBin*numBin) ))^2;
                end;
            end;
        end;
        totCount = sum( cellfun(@(x) sum(x(iChann,:)), testSeq) ); % total spike count
        mrError(iChann,:) = cve; %/totCount;
    end;

end
