function hmmModelSelection( mu, binSize, timeEpoch, maxNumState, numFold, recording, resultSave )
% hmmAnalysisCrossVal( mu, binSize, timeEpoch, maxNumState, numFold, layerData, resultSave )
%
% fits the Hidden Markov Model to the data and performs
% Leave-One-Neuron-Out cross-valadation analysis
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
%   numFold     -   number of cross-validation folds
%   layerData   -   layer-related information about each channel included
%                   in the analysis
%   resultSave  -   data structure that contains field dataDirName, where
%                   the resultd will be saved
%
% saved data:
% ----------
%   file:  resultHMMfit.mat
%   -----------------------
%      estTrans    -   estimated transition matrix, size (numState, numState)
%      estEmis     -   estimated emission matrix, size (numState, numChannel)
%      estPi0      -   estimated initial probability distribution, size
%                      numChannel
%      states      -   cell array of size numTrial, each element contains
%                      the decoded most likely sequence of states on that
%                      trial
%      Pstates     -   cell array of size numTrial, each element is an
%                      array of size (numState, numBin) that contains the
%                      probability of each state throughout each trial
%
%   file:   crossValHMMfit.mat
%   ----------------------------------
%       cvError     -   array of size (numChannel, numFold) with the
%                       leave-one-neuron out cross-validation prediction
%                       error based on HMM fit




    numTrial = size(mu,2);
    numChannel = size(mu,1);
    
    numIter = 5;   %5
    
    % on which channel was the unit recorded
    chann = cellfun(@(x) x(1), recording.msUnitList);
    
    if (size(timeEpoch,1)==1)
        timeEpoch = repmat( timeEpoch, numTrial,1);
    end;
    
    % generate emission sequence from the spike data
    [ emissionSeq, ~ ] = getEmissionSeq( mu, binSize, timeEpoch );


    timeWin = [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]; % time windows to evaluate cv Error, in units of number of bins
    numWin = length(timeWin);
    
    % =============================
    % perform cross-validation
    % =============================

    randTrial = randperm(numTrial); % randomized trials
    chunkEdge = floor(linspace(1, numTrial+1, numFold+1)); % edges of chunkes used for cross-validation
    
    
    cvError = zeros(numChannel, maxNumState, numWin, numFold);   % cross-validation error of the HMM fit
    looError = zeros(numChannel, maxNumState, numWin, numFold);   % cross-validation error of the HMM fit
    veError = zeros(numChannel, maxNumState, numWin, numFold);   % variance-explained error of the HMM fit
    totalCount = zeros(numChannel, numFold);   % total spike count for each channel on each cross-validation fold
        
    dirName = resultSave.dataDirName;
    fileName = [dirName 'selectHMM_maxNumState_' num2str(maxNumState) '.mat'];
    if exist(fileName, 'file')
        fprintf('%s already exists, skipping...\n',fileName)
        return
    end
    for numState = 1:maxNumState
        fprintf('Running state %d/%d\n', numState, maxNumState)
        parfor iFold = 1:numFold
            % Set cross-validation masks
            testMask = false(1, numTrial);        
            testMask(chunkEdge(iFold):chunkEdge(iFold+1)-1) = true;
            trainMask = ~testMask;

            indTrain = randTrial(trainMask);
            indTest = randTrial(testMask);    
            trainSeq = emissionSeq(indTrain);
            testSeq = emissionSeq(indTest);

            % train on train trials
            [estTrans, estEmis, estPi0] = fitHMM( trainSeq, numState, binSize, numIter );
            % compute error on test trials
            [ cvError(:,numState,:,iFold), looError(:,numState,:,iFold), veError(:,numState,:,iFold), totalCount(:,iFold) ] = crossValError(testSeq, estTrans, estEmis, estPi0, chann, timeWin);
            
        end;
    end;
    
    % =============================
    % save data
    % =============================
    % estimate mean firing rate from all data
    meanRate = zeros(numChannel, 1);
    for iChann = 1:numChannel
        meanRate(iChann) = sum( cellfun(@(x) sum( x(iChann,:)), emissionSeq ) )/sum( cellfun(@(x) size(x,2), emissionSeq ) );
    end;
    meanRate = meanRate/binSize;
    timeWindow = timeWin*binSize;
    layerData = recording.layerData;
    %
    if (resultSave.data)
        if (~isdir(dirName))
            mkdir(dirName);
        end;
        save(fileName, 'cvError', 'looError', 'veError', 'totalCount', 'meanRate', 'layerData', 'maxNumState', 'timeWindow', 'binSize');
    end;
    data = load(fileName);
    if (resultSave.figures)
    	plotHMMselection( data, resultSave );
    end
    
    
end























