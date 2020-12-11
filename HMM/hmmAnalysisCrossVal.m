function hmmAnalysisCrossVal( mu, binSize, timeEpoch, minNumState, maxNumState, numFold, recording, resultSave )
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
%       mrError     -   array of size (numChannel, numFold) with the
%                       cross-validation error based on mean firing rate
%       R2          -   coefficient of determination R2, defined as 1-cvError/mrError
%       spError1    -   array of size (numChannel, numFold) with the
%                       cross-validation error based on predicting the
%                       hidden state from activity of a single neuron vs.
%                       from activity of all neurons
%       spError2    -   array of size (numChannel, numFold) with the
%                       cross-validation error based on predicting the
%                       hidden state from activity of all but one neurons
%                       vs. from all neurons

numTrial = size(mu,2);
numChannel = size(mu,1);

numIter = 10;   %5 %number of EM-optimization (Expectation Maximization) random initialisations of HMM parameters (10 in Engel paper).
nBoot = 10;  %20;

% on which channel was the unit recorded
chann = cellfun(@(x) x(1), recording.msUnitList);

if (size(timeEpoch,1)==1)
    timeEpoch = repmat( timeEpoch, numTrial,1);
end;

% generate emission sequence from the spike data
[ emissionSeq, timeBins ] = getEmissionSeq( mu, binSize, timeEpoch );


% =============================
% fit all trials
% =============================

N = length(emissionSeq);
bootsam = unidrnd(N,N,nBoot);

aic = zeros(maxNumState,1);
bic = zeros(maxNumState,1);
aicBoot = zeros(nBoot,maxNumState);
bicBoot = zeros(nBoot,maxNumState);

% fit the HMM model
for numState = minNumState:maxNumState
    dirName = resultSave.dataDirName;
    fileName = [dirName 'HMMfit_numState_' num2str(numState) '.mat'];
    
    if ~exist(fileName,'file') || resultSave.overWriteExistingFiles
        fprintf('testing numstate %1d\n', numState)
        tStart = tic;
        fprintf('estimating parameters\n')
        [estTrans, estEmis, estPi0, aic(numState), bic(numState)] = fitHMM( emissionSeq, numState, binSize, numIter );
        % determine states and state probabilities
        states = cell(numTrial,1);
        Pstates = cell(numTrial,1);
        fprintf('fitting model\n')
        for iTrial = 1:numTrial
            states{iTrial} = hmmviterbiPoisson(emissionSeq{iTrial},estTrans,estEmis, estPi0);
            Pstates{iTrial} = hmmDecodePoisson(emissionSeq{iTrial},estTrans,estEmis, estPi0);
        end;
        
        % analysis of residence time distribution
        [ timeEpochState, ~, ~ ] = extractUpDownTimeEpoch( states, timeBins, numState, binSize, true);
        residTimeCDF = cell(numState,1);
        residTimePDF = cell(numState,1);
        fracTimeState = zeros(numState,1);
        for iState = 1:numState
            timeEpochThisState = vertcat( timeEpochState{:,iState} );
            
            if (~isempty(timeEpochThisState))
                residTime = timeEpochThisState(:,2) - timeEpochThisState(:,1);
                [f, x] = ecdf( residTime );
                residTimeCDF{iState} = [x, f];
                [n, x] = hist(residTime, 15);
                n = n/sum(n);
                residTimePDF{iState} = [x', n'];
                fracTimeState(iState) = sum(residTime);
            else
                residTimeCDF{iState} = zeros(0,2);
                residTimePDF{iState} = zeros(0,2);
                fracTimeState(iState) = 0.0;
            end;
        end;
        fracTimeState = fracTimeState/sum(fracTimeState); % /totalTime
        
        % bootstrap parameters and model selection
        fprintf('bootstrapping\n')
        [ bootstat, bootTrans, bootEmis, bootPi0 ] = bootstrapLL(emissionSeq, bootsam, numState, binSize, numIter);
        aicBoot(:,numState) = bootstat(:,1);
        bicBoot(:,numState) = bootstat(:,2);
        
        estTransBcb = prctile( bootTrans, [5 95] );
        estEmisBcb = prctile( bootEmis, [5 95] );
        estPi0Bcb = prctile( bootPi0, [5 95] );
        
        
        % =============================
        % save data
        % =============================
        if (resultSave.data)
            if (~exist(dirName,'dir'))
                mkdir(dirName);
            end;
            save(fileName, 'estTrans', 'estEmis', 'estPi0', 'aic', 'bic', 'aicBoot', 'bicBoot', 'estTransBcb', 'estEmisBcb', 'estPi0Bcb', 'residTimeCDF', 'residTimePDF', 'fracTimeState',...
                'states', 'Pstates', 'timeBins', 'binSize', 'timeEpoch', 'numState');
        end;
        
        tStop = toc(tStart);
        fprintf(1,'Duration: %1.0fm %1.0fs \n',floor(tStop/60),rem(tStop,60));
        
        if resultSave.figures
            data = load(fileName);
            plotHMMfit( mu, data, resultSave );
        end
    else
        fprintf('loading existing file %s\n',fileName)
        tmp = load(fileName);
        
        if resultSave.figures
            tmp = load(fileName);
            plotHMMfit( mu, tmp, resultSave );
        end
        aicBoot(:,numState) = tmp.aicBoot(:,numState);
        bicBoot(:,numState) = tmp.bicBoot(:,numState);
        aic(numState) = tmp.aic(numState);
        bic(numState) = tmp.bic(numState);
    end
    
end;

% frequency of model selection in bootstrap samples
[~, selectKaic] = min(aicBoot,[],2);
[~, selectKbic] = min(bicBoot,[],2);
nAic = histc(selectKaic,minNumState:maxNumState)/nBoot;
nBic = histc(selectKbic,minNumState:maxNumState)/nBoot;

if resultSave.figures
    plotNaicNbic(aic, bic, nAic, nBic, resultSave);
end

timeWin = [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]; % time windows to evaluate cv Error, in units of number of bins
numWin = length(timeWin);

% =============================
% perform cross-validation
% =============================

randTrial = randperm(numTrial); % randomized trials
chunkEdge = floor(linspace(1, numTrial+1, numFold+1)); % edges of chunkes used for cross-validation

for numState = minNumState:maxNumState
    dirName = resultSave.dataDirName;
    fileName = [dirName 'crossValHMMfit_numState_' num2str(numState) '.mat'];
    
    if ~exist(fileName,'file') || resultSave.overWriteExistingFiles
        fprintf('cross validating numstate %1d\n', numState)

        cvError = zeros(numChannel, numWin, numFold);   % cross-validation error of the HMM fit
        looError = zeros(numChannel, numWin, numFold);   % cross-validation error of the HMM fit
        veError = zeros(numChannel, numWin, numFold);   % variance-explained error of the HMM fit
        mrError = zeros(numChannel, numWin, numFold);   % cross-validation error for prediction based on mean firing rate
        R2 = zeros(numChannel, numWin, numFold);    % R2 for cross-validation
        veR2 = zeros(numChannel, numWin, numFold);    % variance explained R2 for cross-validation
        looR2 = zeros(numChannel, numWin, numFold);    % variance explained R2 for cross-validation
        spError1 = zeros(numChannel, numFold);   % cross-validation error for predicting state based on activity of single neurons
        spError2 = zeros(numChannel, numFold);
        
        for iFold = 1:numFold
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
            % compute error of test trials
            [ cvError(:,:,iFold), looError(:,:,iFold), veError(:,:,iFold) ] = crossValError(testSeq, estTrans, estEmis, estPi0, chann, timeWin);
            mrError(:,:,iFold) = meanRateError(trainSeq, testSeq, timeWin);
            R2(:,:,iFold) = 1 - cvError(:,:,iFold)./mrError(:,:,iFold);
            veR2(:,:,iFold) = 1 - veError(:,:,iFold)./mrError(:,:,iFold);
            looR2(:,:,iFold) = 1 - looError(:,:,iFold)./mrError(:,:,iFold);
            
            [spError1(:,iFold), spError2(:,iFold)] = statePredError(testSeq, estTrans, estEmis, estPi0);
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
            %             if (~isdir(dirName))
            %                 mkdir(dirName);
            %             end;
            save(fileName, 'cvError', 'mrError', 'R2', 'veError', 'looError', 'veR2', 'looR2', 'spError1', 'spError2', 'meanRate', 'layerData', 'numState', 'timeWindow');
        end;
    else
        fprintf('loading existing cross validation %s\n', fileName)

    end
    
    if resultSave.figures
        
        data = load(fileName);
        
        plotHMMcrossVal( data, resultSave );
    end
end;



end





function [ bootstat, bootTrans, bootEmis, bootPi0 ] = bootstrapLL(emissionSeq, bootsam, numState, binSize, numIter)

nBoot = size(bootsam,2);
numChannel = size(emissionSeq{1},1);
bootstat = zeros(nBoot,2); % aic, bic
bootaic = zeros(nBoot,1);
bootbic = zeros(nBoot,1);
bootTrans = zeros(nBoot, numState, numState); % bootTrans
bootEmis = zeros(nBoot, numState, numChannel); % bootEmis
bootPi0 = zeros(nBoot, numState);

parfor iBoot = 1:nBoot
    [bootTrans(iBoot,:,:), bootEmis(iBoot,:,:), bootPi0(iBoot,:), bootaic(iBoot,1), bootbic(iBoot,1)] = ...
        fitHMM( emissionSeq(bootsam(:,iBoot)), numState, binSize, numIter );
end;

bootstat(:,1) = bootaic;
bootstat(:,2) = bootbic;

end





function [spError1, spError2 ]= statePredError(testSeq, estTrans, estEmis, estPi0)
% for each neuron computes the cross-validation error of predicting the
% hidden state, which is the accuracy with each state decoded from
% individual neurons overlaps with state decoded from all neurons
%
% arguments:
% ---------
%   testSeq                     -   emission sequence on the test trials
%   estTrans, estEmis, estPi0   -   parameters of the HMM fitted on the
%                                   train trials
%
% outputs:
% -------
%   spError1    -   array of size numChannel with the cross-vaidation error
%                   of type 1 normalized by bin count; type 1 is a discrepancy between
%                   the hidden state decoded from the activity of single neuron
%                   and decoded from all neurons
%   spError2    -   array of size numChannel with the cross-validation
%                   error of type 2 normalized by bin count; type 2 is a
%                   discrepance between the hidden state decoded with one
%                   neuron left out and decoded from all neurons

numTrial = length(testSeq);
numChannel = size(estEmis,2);
spError1 = zeros(numChannel,1);
spError2 = zeros(numChannel,1);

stateAll = cell(numTrial,1);
% decode state sequence based on all channels
for iTrial = 1:numTrial
    stateAll{iTrial} = hmmviterbiPoisson(testSeq{iTrial},estTrans,estEmis, estPi0);
end;

% decode state from single neurons one by one
for iChann = 1:numChannel
    cve1 = 0;
    cve2 = 0;
    for iTrial = 1:numTrial
        % decode states based on the activity of a single neuron
        emis = estEmis(:, iChann);
        seq = testSeq{iTrial}(iChann, :);
        state = hmmviterbiPoisson(seq,estTrans,emis, estPi0);
        % compare to the state decoded from all neurons
        cve1 = cve1 + sum(state~=stateAll{iTrial}); % number of bins with different decode
        
        % decode states based on the activity of all but one neuron (leave-one-out)
        ind = [1:iChann-1 iChann+1:numChannel];
        emis = estEmis(:, ind);
        seq = testSeq{iTrial}(ind, :);
        % decode states based on the activity of remaining neurons
        state = hmmviterbiPoisson(seq,estTrans,emis, estPi0);
        % compare to the state decoded from all neurons
        cve2 = cve2 + sum(state~=stateAll{iTrial}); % number of bins with different decode
    end;
    totCount = sum( cellfun(@(x) size(x,2), testSeq) ); % total bin count
    spError1(iChann) = cve1/totCount;
    spError2(iChann) = cve2/totCount;
end;

end











