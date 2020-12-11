function crossValBasedOnState( mu, binSize, timeEpoch, numFold, recording, resultSave )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    indHash = cellfun(@(x) strcmp(x,'hash'), recording.dataType);
    hash = mu(indHash,:);
    indMSU = cellfun(@(x) any(strcmpi(x,{'su', 'suinv', 'mu', 'unit'})), recording.dataType);
    msu = mu(indMSU,:);
    % on which channel was the unit recorded
    channHash = cellfun(@(x) x(1), recording.msUnitList(indHash));
    channMSU = cellfun(@(x) x(1), recording.msUnitList(indMSU));
    
    numChannel = size(hash,1);
    numUnit = size(msu,1);
    numTrial = size(hash,2);
    
    numIter = 10;   %5
    
    if (size(timeEpoch,1)==1)
        timeEpoch = repmat( timeEpoch, numTrial,1);
    end;
    
    % generate emission sequence from the spike data
    [ emissionHash, timeBins ] = getEmissionSeq( hash, binSize, timeEpoch );
    [ emissionMSU, timeBins ] = getEmissionSeq( msu, binSize, timeEpoch );
    
    
    timeWin = [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]; % time windows to evaluate cv Error, in units of number of bins
    numWin = length(timeWin);
    
    randTrial = randperm(numTrial); % randomized trials
    chunkEdge = floor(linspace(1, numTrial+1, numFold+1)); % edges of chunkes used for cross-validation
    
    for numState = 2;   %1:maxNumState
        dirName = resultSave.dataDirName;
        fileName = [dirName 'crossValStateBased_numState_' num2str(numState) '.mat'];

        if exist(fileName,'file') && ~resultSave.overWriteExistingFiles
            fprintf('%s already exists\n',fileName) 
            continue
        end
        
        cvError = zeros(numUnit, numWin, numFold);   % cross-validation error of the HMM fit
        veError = zeros(numUnit, numWin, numFold);   % variance-explained error of the HMM fit
        mrError = zeros(numUnit, numWin, numFold);   % cross-validation error for prediction based on mean firing rate
        R2 = zeros(numUnit, numWin, numFold);    % R2 for cross-validation
        veR2 = zeros(numUnit, numWin, numFold);    % variance explained R2 for cross-validation
 
        for iFold = 1:numFold
            % Set cross-validation masks
            testMask = false(1, numTrial);        
            testMask(chunkEdge(iFold):chunkEdge(iFold+1)-1) = true;
            trainMask = ~testMask;

            indTrain = randTrial(trainMask);
            indTest = randTrial(testMask);    
            trainSeqHash = emissionHash(indTrain);
            trainSeqMSU = emissionMSU(indTrain);
            testSeqHash = emissionHash(indTest);
            testSeqMSU = emissionMSU(indTest);

            % train using hash on train trials
            [estTrans, estEmisHash, estPi0] = fitHMM( trainSeqHash, numState, binSize, numIter );
            % decode states based on hash on train trials
            trainStates = cell(length(indTrain),1);
            for iTrial = 1:length(indTrain)
                trainStates{iTrial} = hmmviterbiPoisson(trainSeqHash{iTrial}, estTrans, estEmisHash, estPi0);
            end;
            % estimate firing rate of mu and su in each of the states
            estEmisMSU = estimateRateState(trainSeqMSU, trainStates, timeBins(indTrain), numState, binSize, 'prob');
            
            % compute error of test trials
            [ cvError(:,:,iFold), veError(:,:,iFold) ] = cvErrorStateBased(testSeqHash, testSeqMSU, estTrans, estEmisHash, estPi0, estEmisMSU, channHash, channMSU, timeWin);
            mrError(:,:,iFold) = meanRateError(trainSeqMSU, testSeqMSU, timeWin);
            R2(:,:,iFold) = 1 - cvError(:,:,iFold)./mrError(:,:,iFold);
            veR2(:,:,iFold) = 1 - veError(:,:,iFold)./mrError(:,:,iFold);

        end;
        
        % =============================
        % save data
        % =============================
        % estimate mean firing rate from all data
        meanRate = zeros(numUnit, 1);
        for iChann = 1:numUnit
            meanRate(iChann) = sum( cellfun(@(x) sum( x(iChann,:)), emissionMSU ) )/sum( cellfun(@(x) size(x,2), emissionMSU ) );
        end;
        meanRate = meanRate/binSize;
        timeWindow = timeWin*binSize;
        layerData = recording.layerData(indMSU);
        dataType = recording.dataType(indMSU);
        msUnitList = recording.msUnitList(indMSU);
        %
%         dirName = resultSave.dataDirName;
        if (resultSave.data)
            if (~isdir(dirName))
                mkdir(dirName);
            end;
%             fileName = [dirName 'crossValStateBased_numState_' num2str(numState) '.mat'];
            save(fileName, 'cvError', 'mrError', 'R2', 'veError', 'veR2', 'meanRate', 'layerData', 'dataType', 'msUnitList', 'numState', 'timeWindow');
        end;
        data = load(fileName);
        plotHMMcrossVal( data, resultSave );
        
    end;

end



function [ cvError, veError ] = cvErrorStateBased(testSeqHash, testSeqMSU, estTrans, estEmisHash, estPi0, estEmisMSU, channHash, channMSU, timeWin)
% for each neuron computes the cross-validation error using emission rate
% in each state estimated based on HMM fit
%
% arguments:
% ---------
%   trainSeq                    -   emisison sequence used to estimate firing rate in each state
%   testSeq                     -   emission sequence on the test trials
%   trainStates                 -   sequence of states on train trials decoded from the HMM fit
%   testStates                  -   sequence of states on test trials decoded from the HMM fit
%   numState                    -   number of states in the HMM
%
% outputs:
% -------
%   cvError     - array of size numChannel with the average cross-validation error for each neuron

    numWin = length(timeWin);

    numTrial = length(testSeqHash);
    numChannel = size(testSeqHash{1},1);
    numUnit = size(testSeqMSU{1},1);
    cvError = zeros(numUnit,numWin);
    veError = zeros(numUnit,numWin);
    
    % determine indeces of window sizes for each there is not enough long
    % data for cross-validation, to set them to NaN
    maxNumBin = max( cellfun(@(x) size(x,2), testSeqMSU) );
    indNan = false(size(timeWin));
    indNan(timeWin>maxNumBin) = true;
    
    
    for iUnit = 1:numUnit
        cve = zeros(numWin,1); % cve for each window size
        cve(indNan) = NaN;
        vee = zeros(numWin,1); % cve for each window size
        vee(indNan) = NaN;
        for iTrial = 1:numTrial
            ind = true(numChannel,1); % exclude hash channels from the same and nearby electrodes as the current MS-unit
            ind(channHash==channMSU(iUnit) | channHash==channMSU(iUnit)+1 | channHash==channMSU(iUnit)-1) = false;
            
            emis = estEmisHash(:, ind);
            seq = testSeqHash{iTrial}(ind, :);
            % decode states based on the activity of the remaining hash channels
            state = hmmviterbiPoisson(seq,estTrans,emis, estPi0);
            % decode states based on the activity of all hash channels
            stateFull = hmmviterbiPoisson(testSeqHash{iTrial},estTrans,estEmisHash, estPi0);
            
            % predict firing rate
            predCount = arrayfun(@(x) estEmisMSU(x, iUnit), state);
            predCountFull = arrayfun(@(x) estEmisMSU(x, iUnit), stateFull);
            diffCount = predCount - testSeqMSU{iTrial}(iUnit,:);
            diffCountFull = predCountFull - testSeqMSU{iTrial}(iUnit,:);
            totalNumBin  = length(diffCount);
            for iWin = 1:numWin
                numBin = timeWin(iWin); % number of bins over which to compute error
                for iBin = 1:floor(totalNumBin/numBin)
                    cve(iWin) = cve(iWin) + (sum(diffCount(1+(iBin-1)*numBin:iBin*numBin) ))^2;
                    vee(iWin) = vee(iWin) + (sum(diffCountFull(1+(iBin-1)*numBin:iBin*numBin) ))^2;
                end;
            end;
        end;
        totCount = sum( cellfun(@(x) sum(x(iUnit,:)), testSeqMSU) ); % total spike count
        cvError(iUnit,:) = cve;%/totCount;
        veError(iUnit,:) = vee;%/totCount;
    end;

end







