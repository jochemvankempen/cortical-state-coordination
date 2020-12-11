function [states, seq]= hmmGeneratePoisson(tr, e, pi0, timeBins)
% [states, seq]= hmmGeneratePoisson(tr, e, pi0, timeBins)
% generates a realization of HMM model with Poisson emissions
%
% arguments:
% ---------
%   tr          -   transition matrix for the HMM, size (numState, numState)
%   e           -   emission matrix, size (numState, numChannel), firing
%                   probability for the Poisson emission of each channel lambda*deltaT
%   pi0         -   stationary distribution for the Markov chain
%   timeBins    -   cell array of length numTrial, that contains centers of time
%                   bins at which spike counts will be emitted on
%                   each trial
%
% outputs:
% --------
%   states      -   cell array of length numTrial, contains the generated
%                   sequence of states for each trial
%   seq         -   emission sequence, cell array of length numTrial, which
%                   contains spike counts for each channel in each time bin
%                   on each trial, seq{numTrial}(numChannel, numBin)


    numTrial = length(timeBins);

    seq = cell(numTrial,1);
    states = cell(numTrial,1);

    % tr must be square
    numStates = size(tr,1);
    checkTr = size(tr,2);
    if checkTr ~= numStates
        error(message('stats:hmmgenerate:BadTransitions'));
    end

    % number of rows of e must be same as number of states
    checkE = size(e,1);
    if checkE ~= numStates
        error(message('stats:hmmgenerate:InputSizeMismatch'));
    end
    numChannel = size(e,2);

    % calculate cumulative probabilities
    trc = cumsum(tr,2);
    pi0c = cumsum(pi0);
   
    % normalize these just in case they don't sum to 1.
    trc = trc./repmat(trc(:,end),1,numStates);
    pi0c = pi0c/pi0c(end);
   
    for iTrial = 1:numTrial

        numBin = size(timeBins{iTrial},2);
        states{iTrial} = zeros(1,numBin);
        seq{iTrial} = zeros(numChannel,numBin);

        % choose initial state according to pi0
        currentstate = binarySearch(rand(),pi0c);
 
        % create random numbers for state changes
        statechange = rand(1,numBin);
        
        % calculate states and emissions
        for iBin = 1:numBin
            
            state = binarySearch(statechange(iBin),trc(currentstate,:));
            states{iTrial}(iBin) = state;
            seq{iTrial}(:,iBin) = poissrnd(e(state,:));
            currentstate = state;
            
        end;
    end;

end








