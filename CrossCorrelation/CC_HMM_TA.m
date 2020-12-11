function CC_HMM_TA(msu, timeEpoch, SF, analysisParam, dataType, msUnitList, areaList, resultSave)
% CC_HMM_TA(msu, timeEpoch, SF, analysisParam, dataType, msUnitList, areaList, resultSave)
%
% compute cross-correlations between HMM transitions and various types of
% data, such as spike trains/analog data for a collection of channels 
%
% Parameters
% ----------
% msu : cell 
%     cell array of size (numChannel+2, numTrial) which contains spike
%     times, the last two rows contain states and timeBins
% timeEpoch : array
%     array [t_begin t_end] time epoch 
% SF : float
%     float Sampling Frequency
% analysisParam : struct
%     structure, fields contain analysis parameters
% dataType : cell 
%     cell array that contains data type string for each channel, e.g 'su',
%     'mu' or 'hash' 
% msUnitList : cell 
%     cell array that contains labels of each channel: [chanNum] for hash,
%     and [chanNum unitNum] for su and mu 
% areaList : cell 
%     cell array that contains which area the unit/channel was recorded
%     from 
% resultSave : struct 
%     structure that contains information relevant for saving and
%     plotting the result 
%
% Returns
% -------
% crossCorr : array
%     array of size (numChan x length(timeLags)) with cross correlations
% meanRateDU : array 
%     array of size (numChannel, 1), mean transition rate 
% timeBinsLag : array
%     time lag axis
% labelCombinations : cell
%     cell with string combinations of labels for trigger and triggered channel
% label : cell
%     cell with label of state transitions (e.g. On-Off)
% msUnitList : cell 
%     cell array that contains labels of each channel: [chanNum] for hash,
%     and [chanNum unitNum] for su and mu 
% areaList : cell 
%     cell array that contains which area the channel was recorded from
%
%
% **HMM_TTA_[CC_type].mat**
%     file that includes all the variables listed under *Returns*. CC_type is set in
%     analysisParam and can be 'cpt' (coincidences per trigger) or
%     'jitter'. If 'jitter', jitterTimeWindow is included in filename 
% 

fprintf('Performing HMM transition triggered average\n')

if size(timeEpoch,1)>1
    %error('Only takes a fixed timeEpoch')
end

binSize = analysisParam.binSize;
maxTimeLag = analysisParam.maxTimeLag;

numState = analysisParam.numState;
numChannel = size(msu,1)-numState;
numTrial = size(msu,2);


timeBins = msu(numChannel+1,:)';
states = msu(numChannel+2,:)';

[ stateChangeTimes, label ] = HMM_extractUpDownTransitionTimes( states, timeBins, analysisParam.numState, binSize, [], analysisParam.minEpochDurationWin, 0);
% 
% [ timeEpochState, durationState, numTransition ] = extractUpDownTimeEpoch( states, timeBins, analysisParam.numState, binSize, 1);
% [ timeEpochState ] = HMM_transformUpDownTimeEpoch( stateChangeTimes, numTransition );
% 

% correct for numChannl, in analog signal the last channel (other than HMM) are the TimeStamps
switch dataType
    case 'HMM'
        
        if numChannel>1
            error('not implemented yet')
        end
        [ tmp, msUnitList ] = HMM_extractUpDownTransitionTimes( msu(1,:)', timeBins, analysisParam.numState, binSize, timeEpoch, analysisParam.minEpochDurationWin);
        
        msu = tmp';
        numChannel = numChannel+1;
        
        areaList = repmat(areaList,[2,1]);
        
        timeBinsLag = -maxTimeLag:binSize:maxTimeLag;
        numBin = length(timeBinsLag);
        
    case {'NCS','LFP','LFPb','MUAe','hash_conv','unit_conv'}
        numChannel = size(msu,1)-numState-1;
        TimeStamps = msu(numChannel+1,:)';
        
        for ichannel = 1:numChannel
            analogSignal.Samples(ichannel,:,:) = cell2mat(msu(ichannel,:)');
        end
        analogSignal.Fs = SF;
        analogSignal.TimeStamps = TimeStamps{1};

        timeBinsLag = [fliplr(1000/SF:1000/SF:(maxTimeLag*1000)).*(-1) 0 [1000/SF:1000/SF:(maxTimeLag*1000)]]/1000;
        numBin = length(timeBinsLag);
        
    case {'hash','unit'}
        numChannel = size(msu,1)-numState;
        TimeStamps = cell(numTrial,1);
        for itrial = 1:numTrial
            TimeStamps{itrial} = timeEpoch(1):(1/SF):timeEpoch(2);
        end
        
        timeBinsLag = -maxTimeLag:binSize:maxTimeLag;
        numBin = length(timeBinsLag);
end

if istable(stateChangeTimes)
    % consider individual transitions as separate trials to allow computing
    % cross correlations from/until the previous/next transition
    msu = msu(1:numChannel, stateChangeTimes.trial);
    
    for iChannel = 1:analysisParam.numState
        idx_state = stateChangeTimes.UD==iChannel;
        
        msuState = cell(1,height(stateChangeTimes));
        msuState(idx_state) = num2cell(stateChangeTimes.transition(idx_state))';
        msu = [msu ; msuState];
    end
    timeEpoch = [stateChangeTimes.previous stateChangeTimes.next];
else
    msu = [msu(1:numChannel,:) ; stateChangeTimes'];
end

crossCorr = zeros(numChannel*analysisParam.numState, numBin);
straightCrossCorr = zeros(numChannel*analysisParam.numState, numBin);
trueExpectPredictorCross = zeros(numChannel*analysisParam.numState, numBin);
trueExpectPredictorCrossVar = zeros(numChannel*analysisParam.numState, numBin);

labelCombinations = cell(numChannel*analysisParam.numState, 1);

meanRateDU = zeros(analysisParam.numState, 1);

% compute cross-correlations
k = 1;
for iChannel = 1:analysisParam.numState
    for jChannel = 1:numChannel
        pickChannel = [(numChannel+iChannel) jChannel];
        
        switch dataType
            case {'NCS','LFP','LFPb','MUAe','hash_conv','unit_conv'}
                
                switch analysisParam.CC_type
                    case 'cpt'
                        tmpAnalogSignal = analogSignal;
                        if istable(stateChangeTimes)
                            tmpAnalogSignal.Samples = tmpAnalogSignal.Samples(pickChannel(2),stateChangeTimes.trial,:);
                            WinExtendFlag=0;
                            [STA,BinCentres,SpikeNum,SP,CentreBinNr] = CC_SpikeTA(tmpAnalogSignal,msu(pickChannel(1),:)',maxTimeLag,timeEpoch,WinExtendFlag,stateChangeTimes.trial);
                        else
                            tmpAnalogSignal.Samples = tmpAnalogSignal.Samples(pickChannel(2),:,:);
                            WinExtendFlag=1;
                            [STA,BinCentres,SpikeNum,SP,CentreBinNr] = CC_SpikeTA(tmpAnalogSignal,msu(pickChannel(1),:)',maxTimeLag,timeEpoch,WinExtendFlag);
                        end
                        
                        meanRate_iChan = sum( cellfun(@(x) length(x(x>=timeEpoch(1,1) & x<=timeEpoch(1,2))), msu(pickChannel(1),:) ) ) ...
                            /(size(msu,2)*(timeEpoch(1,2)-timeEpoch(1,1)));
                        
                        if analysisParam.demean
                            % get rid of DC offset
                            STA = STA - repmat(nanmean(STA,2), [1 length(BinCentres)]);
                            SP = SP - repmat(nanmean(SP,2), [1 length(BinCentres)]);
                        end
                        
                        % normalize to Coincidences/spike
                        if istable(stateChangeTimes)
                            STA=STA .* (repmat(ones(1,length(BinCentres))./sum(~isnan(STA),1), [size(STA,1) 1]));
                            SP=SP .* (repmat(ones(1,length(BinCentres))./sum(~isnan(SP),1), [size(SP,1) 1]));
                        else
                            %STA=STA*(1/numTrial)*(1/meanRate_iChan);
                            %SP=SP*(1/numTrial)*(1/meanRate_iChan);
                            STA=STA*(1/sum(SpikeNum));
                            SP=SP*(1/sum(SpikeNum));
                        end
                        % correct shift predictor
                        crossCorr(k,:) = nansum(STA,1) - nansum(SP,1);
                        
                    case 'jitter'
                        
                        error('jitter for analog signal is not implemented yet')
                end

            case {'HMM','hash','unit'}
                
                switch analysisParam.CC_type
                    case 'cpt'
                
                        %%% cpt, coincidences per trigger
                        [STA,BinCentres,SpikeNum,SP,CentreBinNr] = CC_SpikeSpikeTA(msu(pickChannel(2),:)', msu(pickChannel(1),:)', binSize, maxTimeLag, timeEpoch);
                        
                        meanRate_iChan = sum( cellfun(@(x) length(x(x>=timeEpoch(1,1) & x<=timeEpoch(1,2))), msu(pickChannel(1),:) ) ) ...
                            /(size(msu,2)*(timeEpoch(1,2)-timeEpoch(1,1)));
                        meanRate_jChan = sum( cellfun(@(x) length(x(x>=timeEpoch(1,1) & x<=timeEpoch(1,2))), msu(pickChannel(2),:) ) ) ...
                            /(size(msu,2)*(timeEpoch(1,2)-timeEpoch(1,1)));
                        
                        switch analysisParam.normalisation
                            case 'cps'
                                % normalize to Coincidences/spike
                                STA=STA*(1/numTrial)*(1/sqrt(meanRate_iChan*meanRate_jChan));
                                SP=SP*(1/numTrial)*(1/sqrt(meanRate_iChan*meanRate_jChan));
                            case 'cpt'
                                % normalize to spike counts/trigger
                                STA=STA*(1/sum(SpikeNum));
                                SP=SP*(1/sum(SpikeNum));
                            case 'Hz'
                                % normalize to Firing rate [Hz]/trigger
                                STA=STA*(1/binSize)*(1/sum(SpikeNum));
                                SP=SP*(1/binSize)*(1/sum(SpikeNum));
                        end
                        % correct shift predictor
                        crossCorr(k,:) = sum(STA,1) - sum(SP,1);
                        
                    case 'jitter'
                        straightCrossCorr(k,:) = CC_computeStraightCrossCorrBiased( msu, timeEpoch, pickChannel, binSize, maxTimeLag );
                        [  trueExpectPredictorCross(k,:), trueExpectPredictorCrossVar(k,:) ] = ...
                            CC_trueExpectedPredStraightCrossCorrJitterPSTHbiased( msu, pickChannel, timeEpoch, maxTimeLag, binSize, analysisParam.jitterBinSize);
                        
                        crossCorr(k,:) = straightCrossCorr(k,:) - trueExpectPredictorCross(k,:);
                        ind = trueExpectPredictorCrossVar(k,:)~=0;
                        crossCorr(k,ind) = crossCorr(k,ind)./sqrt(trueExpectPredictorCrossVar(k,ind));
                        
                end
        end
                
        if iscell(msUnitList(jChannel))
            labelCombinations{k,1} = [areaList{jChannel} msUnitList{jChannel} '_' label{iChannel}];
        else
            labelCombinations{k,1} = [char(areaList(jChannel,:)) '_' num2str(msUnitList(jChannel)) '_' label{iChannel}];
        end
        k = k + 1;
        
    end;
    
    meanRateDU(iChannel) = sum( cellfun(@(x) length(x(x>=timeEpoch(1,1) & x<=timeEpoch(1,2))), msu((numChannel+iChannel),:) ) ) ...
        /(size(msu,2)*(timeEpoch(1,2)-timeEpoch(1,1)));
end;


dirName = resultSave.dataDirName;
% save data
if (resultSave.data)
    if (~isdir(dirName))
        mkdir(dirName);
    end;
    %fileName = [dirName 'spikeTransTA_sp_binsize(' num2str(analysisParam.binSize) ')_epoch(' num2str(analysisParam.minEpochDurationWin) ')_win(' num2str(analysisParam.timechop) ').mat'];
    switch analysisParam.CC_type
        case 'cpt'
            fileName = [dirName sprintf('HMM_TTA_cpt.mat')];
        case 'jitter'
            fileName = [dirName sprintf('HMM_TTA_jitter%1.3f.mat',analysisParam.jitterBinSize)];

    end
    % save auto- and cross-correlations
    %save(fileName, 'allTA', 'stateswitchLabels', 'timeWin', 'dataType', 'msUnitList','areaList');
    
    save(fileName, 'crossCorr*', ...
        'meanRateDU', 'timeBinsLag', 'labelCombinations', 'label', 'msUnitList', 'areaList');

end;


