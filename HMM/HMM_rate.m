function HMM_rate(msu, analysisParam, dataType, msUnitList, areaList, resultSave)
% HMM_rate(msu, analysisParam, dataType, msUnitList, areaList, resultSave)
%
% compute firing rate across HMM states
%
% Parameters
% ----------
% msu : cell 
%     cell array of size (numChannel+2, numTrial) which contains spike
%     times, two rows contain HMM timebins and states
% analysisParam : struct
%     structure, fields contain analysis parameters
% dataType : cell
%     cell with string indicating type of data
% msUnitList : cell 
%     cell array that contains labels of each channel: [chanNum] for hash,
%     and [chanNum unitNum] for su and mu 
% areaList : cell 
%     cell array that contains which area the channel was recorded from
% resultSave : struct
%     structure that contains information relevant for saving and plotting
%     the result 
%
% Returns
% -------
% rate : array
%     array of size (numChannel x numState x time) with firing rate 
% sRate : array
%     array of size (numChannel x numState x time) with firing rate
%     jackknife error bars
% time : array
%     time points to compute rate. Here set to 1 as we apply this to HMM
%     epochs
% dataType : cell
%     cell with string indicating type of data
% msUnitList : cell 
%     cell array that contains labels of each channel: [chanNum] for hash,
%     and [chanNum unitNum] for su and mu 
% areaList : cell 
%     cell array that contains which area the channel was recorded from
%
% 
% **HMM_rate.mat : file**
%     file that includes all the variables listed under *Returns*
% 
%
%

fprintf('Computing rate per HMM state\n')

numstate = 2;

binSize = analysisParam.binSize;

numChannel = size(msu,1)-numstate;
numTrial = size(msu,2);

timeBins = msu(numChannel+1,:)';
timeBins = cellfun(@(x)(x-0.005),timeBins,'UniformOutput',false);
timeBins = cellfun(@(x)([x x(end)+0.01]),timeBins,'UniformOutput',false);

states = msu(numChannel+2,:)';

[ timeEpochState, durationState, numTransition ] = extractUpDownTimeEpoch( states, timeBins, analysisParam.numState, binSize, analysisParam.inclBorder);
[ timeEpochState ] = HMM_transformUpDownTimeEpoch( timeEpochState, numTransition );

rate = zeros(numChannel, analysisParam.numState, 1);
sRate = zeros(numChannel, analysisParam.numState, 1);

for istate = 1:analysisParam.numState

    % only use epochs with minimum duration
    epochIdx = ( (timeEpochState{istate}(:,3)-timeEpochState{istate}(:,2)) >= analysisParam.minEpochDurationWin );
    timeEpochState{istate} = timeEpochState{istate}(epochIdx,:);
    
    nEpoch = length(find(epochIdx));
    
    % create msu matrix for each epoch, rather than trial
    msu_epoch = cell(numChannel, nEpoch);
    for iepoch = 1:nEpoch
        trialIdx = timeEpochState{istate}(iepoch,1);
        msu_epoch(:,iepoch) = msu(1:numChannel,trialIdx);
    end;
    
    clear window
    window(:,1,:) = timeEpochState{istate}(:,2:end);
    time = 1;
    
    for iChannel = 1:numChannel
        % calculate rate
        [rate(iChannel,istate,:), sRate(iChannel,istate,:)] = spike_rate(  msu_epoch, iChannel, time, window );
        
    end;
end;


% save data
dirName = resultSave.dataDirName ;
if (resultSave.data)
    if (~isdir(dirName))
        mkdir(dirName);
    end;
    fileName = [dirName 'HMM_rate.mat'];
    save(fileName, 'rate', 'sRate', 'time', 'dataType', 'msUnitList', 'areaList');
end;







