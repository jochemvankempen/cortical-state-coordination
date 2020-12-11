function HMM_correlateEpochDuration(msu, analysisParam, dataType, resultSave)
% HMM_correlateEpochDuration(msu, analysisParam, dataType, resultSave)
%
% compute correlations with epoch durations and fraction of time spent in state
% 
% Parameters
% ----------
% msu : cell 
%     cell array of size (numChannel+2, numTrial) which contains spike
%     times, the last row contains microsaccade times 
% analysisParam : struct
%     structure, fields contain analysis parameters
% dataType : cell 
%     cell array that contains data type string for each channel, 'su',
%     'mu', 'hash', microsaccades, pupil
% resultSave : struct
%     structure that contains information relevant for saving and plotting
%     the result 
%
% Returns
% -------
% HMM_corr.mat : file
%     file that contains all the variables below
% R_fractimeState : array
%     array of size (numChannel x numState) with correlations for with
%     fractimestate
% R_fractimeStateShuffle : array
%     array of size (numChannel x numState) with correlations for with
%     fractimestate when trials are shuffled
% R_epochDuration : array
%     array of size (numChannel x numState) with correlations for with
%     epochDuration
% R_epochDurationShuffle : array
%     array of size (numChannel x numState) with correlations for with
%     epochDuration when trials are shuffled
% dataType : cell
%     cell of size (numChannel,1) with strings indicating type of data
%     correlated to HMM stats 
% epochDuration : cell
%     cell of size (numState,1) with epoch durations per state.
% fracTimeState : cell
%     cell of size (numState,1) with fracTimeState per state.
% msu : cell
%     cell with raw data
%
%
% **HMM_corr.mat : file**
%     file that includes all the variables listed under *Returns*
% 
% 

fprintf('\t\tComputing correlation HMM epoch duration\n')

binSize = analysisParam.binSize;

numChannel = size(msu{1},1);
numTrial = size(msu{1},2);

timeBins = msu{2}(1,:)';
timeBins = cellfun(@(x)(x-binSize/2),timeBins,'UniformOutput',false);
timeBins = cellfun(@(x)([x x(end)+binSize]),timeBins,'UniformOutput',false);

states = msu{2}(2,:)';

[ timeEpochState, durationState, numTransition ] = extractUpDownTimeEpoch( states, timeBins, analysisParam.numState, binSize, analysisParam.inclBorder);
% [ timeEpochState ] = transformUpDownTimeEpoch( timeEpochState, numTransition );

R_fractimeState = NaN(numChannel, analysisParam.numState);
R_fractimeStateShuffle = NaN(numChannel, analysisParam.numState);
R_epochDuration = NaN(numChannel, analysisParam.numState);
R_epochDurationShuffle = NaN(numChannel, analysisParam.numState);

epochDuration = cell(analysisParam.numState,1);
fracTimeState = cell(analysisParam.numState,1);

for istate = 1:analysisParam.numState
    

    fracTimeState{istate} = durationState(:,istate) ./ sum(durationState,2); % fraction of time spent in On state per trial

    epochDuration{istate} = zeros(numTrial,1);
    trIdx = ~cellfun(@(X) isempty(X), timeEpochState(:,istate));
    epochDuration{istate}(trIdx) = cellfun(@(X) (mean( X(:,2) - X(:,1)) ), timeEpochState(trIdx,istate), 'UniformOutput', true);

    
    for iChannel = 1:numChannel
        
        % calculate correlation        
        jackCorr = jackknife(@corr, msu{1}(iChannel,:), fracTimeState{istate});
        jackCorrShuffle = jackknife(@corr, msu{1}(iChannel,:), fracTimeState{istate}(randperm(numTrial)));
                
        R_fractimeState(iChannel,istate) = nanmean(jackCorr,1);
        R_fractimeStateShuffle(iChannel,istate) = nanmean(jackCorrShuffle,1);
        
        jackCorr = jackknife(@corr, msu{1}(iChannel,:), epochDuration{istate});
        jackCorrShuffle = jackknife(@corr, msu{1}(iChannel,:), epochDuration{istate}(randperm(numTrial)));

        R_epochDuration(iChannel,istate) = nanmean(jackCorr,1);
        R_epochDurationShuffle(iChannel,istate) = nanmean(jackCorrShuffle,1);
        
    end;
end;


% save data
dirName = resultSave.dataDirName ;
if (resultSave.data)
    if (~isfolder(dirName))
        mkdir(dirName);
    end;
    fileName = [dirName 'HMM_corr.mat'];
    save(fileName, 'R_fractimeState', 'R_fractimeStateShuffle', 'R_epochDuration', 'R_epochDurationShuffle', 'dataType', 'epochDuration', 'fracTimeState', 'msu');
end;







