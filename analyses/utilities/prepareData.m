function [outData, timeEpoch, varargout] = prepareData(recInfo, areaInfo, analysisParam, varargin)
% [outData, timeEpoch, varargout] = prepareData(recInfo, areaInfo, analysisParam, varargin)
% 
% Prepare data for analyis
% 
% Example
% -------
% ::
% 
%     [outData, timeEpoch, varargout] = prepareData(recInfo, analysisParam, 'Event', 'Fix', 'Condition', 2, 'Channels', 3:14)
%     [outData, timeEpoch, varargout] = prepareData(recInfo, analysisParam, 'Event', 'Fix', 'Condition', 2, 'Area', 'V1')
%
% Parameters
% ----------
% recInfo : table
%     table with single row (one recording)
% areaInfo : struct
%     struct with info about each area, obtained by 'getAreaInfo'
% analysisParam : struct
%     struct with analysis parameters, obtained by 'specifyAnalysisParam
% varargin : optional
% 
%     - 'Event' : Event name, e.g. 'Fix', 'Cue'
%     - 'Condition' : condition index
%     - 'Area' : area index
%     - 'Channels' : array with channel index, logical/double array
%
% Returns
% -------
% outData : struct
%     struct with fields :
%     
%     - *Event*Align : cell array of size (numChannel x numTrial) with data
%     aligned to specified event
%     - numUnit : number of channels/units
%     - unitList : channel that the unit was recorded on
%     - area : list with areas per channel/unit
%
% timeEpoch : array
%     array of size (numTrial, 2) with time epochs for each trial defined
%     in 'analysisParam'
% varargout : cell, optional
% 
%     - {'HMM', HMM states or transitions}
% 



%% check input
assert(height(recInfo)==1, 'height recInfo ~= 1')

%% inspect varargin
% check varargin for selection criteria
flag.event = find(strcmpi(varargin, 'Event'))+1;
flag.condition = find(strcmpi(varargin, 'Condition'))+1;
flag.area = find(strcmpi(varargin, 'Area'))+1;
flag.channel = find(strcmpi(varargin, 'Channels'))+1;

% unpack varargin, get event
if ~isempty(flag.event)
    idx.event = varargin{flag.event};
else
    error('Need event selection')
end

% unpack varargin, get condition
if ~isempty(flag.condition)
    idx.condition = varargin{flag.condition};
else
    % for these analysisTypes, condition selection is necessary
    switch analysisParam.analysisType
        case {'HMM','HMMselect','HMM_2area',...
                'HMM_crossCorrelation','HMM_crossCorrelationPupil',...
                'HMM_cvStateBased',...
                'HMM_transitionTime',...
                'HMM_transitionTriggeredAverage'}
            error('No condition selected')
    end
end

% unpack varargin, get area
if ~isempty(flag.area)
    idx.area = varargin{flag.area};
else
    % for these analysisTypes, area selection is necessary
    switch analysisParam.analysisType
        case {'HMM_crossCorrelation'}
            error('No area selected')
        otherwise
            % warning('No area selected')
    end
end

% unpack varargin, get channels
if ~isempty(flag.channel)
    idx.channels = varargin{flag.channel};
end

varargout = {};
%% load trialdata, event/condition/microsaccade indices and time epoch

% first get trialdata and trial indices
load( fullfile(recInfo.path_read,'trialdata.mat') );
trialdata = selectStabilityWindow(recInfo, trialdata, 'trialData');% only select stable periods

switch analysisParam.analysisType
    
    case {'HMM_RT'}
        %%% if getting RT, get it before eliminating/selecting trials. If
        %%% performing z-score, you want to do this before subselection.
        varargout{1} = ([trialdata.targetDim])';
        
        RT = [trialdata.(analysisParam.RT2use) ];
        trIdx.eliminate = (RT > analysisParam.RTlim);
        
        RT(trIdx.eliminate) = NaN;

        if ~analysisParam.concatenateAcrossEvents
            switch analysisParam.RT2transform
                case 'zscore_log'
                    RT(~trIdx.eliminate) = zscore(log(RT(~trIdx.eliminate)));
            end
        end
end

% get trial indices for specific events
trIdx.event = getEventSpecificTrialIdx(recInfo, trialdata, idx.event);

% get trials that belong to condition
ievnt = find(strcmpi(analysisParam.event.Epoch, idx.event));
[condNames, ~, trIdx.cond] = getConditions(recInfo, analysisParam, analysisParam.event.ConcatenateCond(ievnt), trialdata);
trIdx.cond_true = trIdx.cond;
timeEpoch_allTrials = selectTimeEpoch(trialdata, analysisParam.event(ievnt,:), analysisParam);

switch analysisParam.analysisType
    case {'HMM_RT'}
        varargout{2} = trIdx.cond(:);
end

% exclude trials with microsaccades
[trIdx.ms, dataWarning] = getMicroSaccadeTrialIdx(recInfo, analysisParam, trialdata, ievnt, idx.event);

% select trials
if ~isempty(flag.condition)
    trIdx.cond = (trIdx.cond==idx.condition);
else
    trIdx.cond = true(length(trialdata),1);
end
trialdata = selectTrials(trialdata, 'trialData', (trIdx.event & trIdx.cond & ~trIdx.ms) );

% get times epoch
timeEpoch = selectTimeEpoch(trialdata, analysisParam.event(ievnt,:), analysisParam);

if dataWarning
    outData=[];
    return
end
       

%% now load eye/neural data

switch analysisParam.analysisType
    case {'rate','FF','noiseCorr',... % spike
            'HMM','HMMselect','HMM_2area',... HMM fit
            'HMM_microsaccades', 'HMM_pupil', 'HMM_rate', 'HMM_transitionTriggeredAverage', 'HMM_spectrogram', ... % HMM - data
            'microsaccade_analyse_parameters', ...
            }
        
        outData = getData(recInfo, analysisParam, trIdx, flag, idx, timeEpoch_allTrials);
        if isempty(outData)
            varargout = {[]};
            return 
        end
        
    case 'HMM_cvStateBased' % HMM CV

        assert( strcmpi(analysisParam.dataType,'unit'), 'wrong datatype for HMM_cvStateBased, needs to be "unit", "mu" or "su"')
        
        tmpAnalysisParam = analysisParam;

        % get hash for HMM training
        tmpAnalysisParam.dataType = 'hash';
        outData = getData(recInfo, tmpAnalysisParam, trIdx, flag, idx, timeEpoch_allTrials);
        outData.dataType = repmat(tmpAnalysisParam.dataType, [outData.numUnit 1]);
        
        % get unit/mu/su for HMM R2
        outData_unit = getData(recInfo, analysisParam, trIdx, flag, idx, timeEpoch_allTrials);
        %outData_unit.dataType = repmat(analysisParam.dataType, [outData_unit.numUnit 1]);
        outData_unit.dataType = outData_unit.unitClassification.Unit;

        % concat
        outData.([idx.event 'Align']) = vertcat(outData.([idx.event 'Align']), outData_unit.([idx.event 'Align']));
        outData.numUnit = outData.numUnit + outData_unit.numUnit;

        outData.channel = vertcat(outData.unitList, outData_unit.channel);
        outData.unitList = vertcat(outData.unitList, outData_unit.unitList);
        outData.area = vertcat(outData.area, outData_unit.area);
        outData.dataType = vertcat(outData.dataType, outData_unit.dataType);
                
    case {'HMM_crossCorrelation','HMM_crossCorrelationPupil', 'HMM_transitionTime'}

        HMM = getHMM(recInfo, analysisParam, idx.event, condNames{idx.condition}, idx.area);

        %         trIdx.ms = selectTrials(trIdx.ms, 'trialData', (trIdx.event & trIdx.cond) );
        
        %         outData.states = HMM.states(~trIdx.ms);
        %         outData.timeBins = HMM.timeBins(~trIdx.ms);
        %
        outData.states = HMM.states;
        outData.timeBins = HMM.timeBins;
        
        switch analysisParam.analysisType
            case {'HMM_crossCorrelationPupil'}
                tmpAnalysisParam = analysisParam;
                tmpAnalysisParam.event = tmpAnalysisParam.pupilEvent;
                tmpAnalysisParam.dataType = 'pupil';
                
                pupilDat = getData(recInfo, tmpAnalysisParam, trIdx, flag, idx);
                
                varargout{1} = {'pupil',pupilDat};
        end
        return

    case 'HMM_RT'
        
        %%% make subselection of RT, same as trialdata
        RT = selectTrials(RT, 'trialData', (trIdx.event & trIdx.cond & ~trIdx.ms) );
        
        outData = RT;
        outData = outData(:);
        
    otherwise
        fprintf('Unknown analysis option\n')
        outData = [];
end


% also load HMM states for specific analyses, output via varargout
switch analysisParam.analysisType
    case {'HMM_microsaccades', 'HMM_pupil', 'HMM_rate', 'HMM_spectrogram', 'HMM_transitionTriggeredAverage'}
        
        %trIdx.ms = selectTrials(trIdx.ms, 'trialData', (trIdx.event & trIdx.cond) );
        
        switch analysisParam.HMMtype
            case 'HMM'
                
                for iarea = 1:length(areaInfo)
                    HMM = getHMM(recInfo, analysisParam, idx.event, condNames{idx.condition}, areaInfo(iarea).name);
                    
                    %outData2.area{iarea} = areaInfo(iarea).name;
                    %outData2.states(:,iarea) = HMM.states(~trIdx.ms);
                    outData2.area{iarea} = areaInfo(iarea).name;
                    outData2.states(:,iarea) = HMM.states;

                end
                %outData2.timeBins = HMM.timeBins(~trIdx.ms);
                outData2.timeBins = HMM.timeBins;
                varargout{1} = {'HMM', outData2};
                
            case 'HMM_2area'
                HMM = getHMM(recInfo, analysisParam, idx.event, condNames{idx.condition}, sprintf('%s_%s',areaInfo(1).name,areaInfo(2).name));
                
                %outData2.states = HMM.states(~trIdx.ms);
                %outData2.timeBins = HMM.timeBins(~trIdx.ms);
                outData2.states = HMM.states;
                outData2.timeBins = HMM.timeBins;
                varargout{1} = {'HMM', outData2};

        end
end

% convert mat 2 cell, used for analog signal
switch analysisParam.analysisType
    case {'HMM_transitionTriggeredAverage'}
        switch analysisParam.dataType
            case {'NCS','LFP','LFPb','MUAe'}
                outData = convertMat2Cell(outData);
        end
end

% convolve spikes with Gaussian
if analysisParam.convoluteSpikes
    switch analysisParam.dataType
        case {'hash', 'unit','MUA_spont20'}
            minmaxTimeEpoch = repmat( [min(timeEpoch(:,1)) max(timeEpoch(:,2))], [size(timeEpoch,1), 1]);
            outData = convoluteSpikes(outData, minmaxTimeEpoch + [-analysisParam.maxTimeLag analysisParam.maxTimeLag], analysisParam.convoluteSpikesSigma);
    end
end

switch analysisParam.analysisType
    case 'HMM_microsaccades'
        
        if isfield(analysisParam, 'HMMbinSize')
            binSize = analysisParam.HMMbinSize;
        else
            binSize = analysisParam.binSize;
        end

        [outData2.msu, outData2.msulabel] = getHMM_transitionTimes(outData2.states, outData2.timeBins, analysisParam.numState, binSize, analysisParam.HMMtype);
        idx = find(strcmpi(varargout{1},'HMM'));
        varargout{1}{idx+1} = outData2;% overwrite HMM states with HMM_transitions
end

%%% select trials for varargout
if exist('varargout','var')
    for ivar = 1:length(varargout)
        try
            varargout{ivar} = selectTrials(varargout{ivar}, 'trialData', (trIdx.event & trIdx.cond & ~trIdx.ms) );
        catch
%             tmpfields = fieldnames(varargout{ivar});
%             for ifields = 1:length(tmpfields)
%                 varargout{ivar}.(tmpfields{ifields}) = selectTrials(varargout{ivar}.(tmpfields{ifields}), 'HMM', (trIdx.event & trIdx.cond) );
%             end
        end
    end
end

% set varargout that doesn't require trial selection
positionStimuli = table(trialdata(1).positionRF, trialdata(1).positionOut1, trialdata(1).positionOut2, 'VariableNames', {'RF','Out1','Out2'});

nvar = length(varargout);
varargout{nvar+1} = {'positionStimuli',positionStimuli};

end
%% extra functions

%%% ------------------------------------------------------------------- %%%
function outData = getData(recInfo, analysisParam, trIdx, flag, idx, timeEpoch)

    readName = fullfile(recInfo.path_read, sprintf('%s.mat', analysisParam.dataType));
    switch analysisParam.dataType
        case {'NCS','LFP','LFPb','MUAe'}
            readData = load(readName, [idx.event 'Align'],'numUnit','unitList','area','SF');
        case {'hash','MUA_spont20'}
            readData = load(readName, [idx.event 'Align'],'numUnit','unitList','area');
        case {'unit'}
            readData = load(readName, [idx.event 'Align'],'numUnit','unitList','area','channel','unitClassification');
        case {'microsaccades','pupil'}
            readData = load(readName, [idx.event 'Align']);
    end
    
    switch analysisParam.dataType
        case {'pupil'}
            emptystruct = ~isfield( readData, [idx.event 'Align'] );
            
            if emptystruct
                fprintf('No pupil data for event: %s, skipping...\n', idx.event)
                outData=[];
                return
            end
                
            readData = getPupil(readData, analysisParam, idx.event);
    end
    readData = selectStabilityWindow(recInfo, readData, analysisParam.dataType);% only select stable periods
    
    % discard units that show drift
    switch analysisParam.dataType
        case {'unit'}
            
            if analysisParam.excludeUnits
                stimData = load(readName, ['StimAlign'],'numUnit','unitList','area');
                stimData = selectStabilityWindow(recInfo, stimData, analysisParam.dataType);% only select stable periods

                readData = excludeUnits(stimData, readData, analysisParam, timeEpoch);
            end
    end
    
    % discard spikes to match rate across conditions
    switch analysisParam.dataType
        case {'hash','unit','MUA_spont20'}
            if analysisParam.matchRateAcrossCond
                
                trials2equate = trIdx.cond_true;
                trials2equate(~trIdx.event | trIdx.ms) = NaN;
                
                readData = spike_matchRateAcrossCond(readData, trials2equate, timeEpoch, idx.condition);
            end
    end
    
    % select trials
    readData = selectTrials(readData, analysisParam.dataType, (trIdx.event & trIdx.cond & ~trIdx.ms) );
    
    % select channels/area
    switch analysisParam.dataType
        case 'pupil'
            % no channel selection
            outData = readData;            
        case {'microsaccades'}
            % no channel selection
            outData = readData;
        otherwise
            %%% only select included channels
            if ~isempty(flag.channel)
                outData = selectChannels(readData, analysisParam.dataType, 'channels', idx.channels);
            elseif ~isempty(flag.area)
                outData = selectChannels(readData, analysisParam.dataType, 'area', idx.area);
            end
    end
end

%%% ------------------------------------------------------------------- %%%
function trIdxEvent = getEventSpecificTrialIdx(recInfo, inData, event)
% find trials in which certain events occurred 

    switch recInfo.Task
        case 'gratc'
            %%% If event is dim2 or 3, select trials where dim2 occurred
            switch event
                case {'Dim2','Dim3'}
                    switch event
                        case 'Dim2'
                            fprintf('\tFocusing on time period after first dimming, eliminating trials where stimulus in RF or target stimulus dimmed at first dimming\n')
                            trIdxEvent = ([inData.targetDim] > 1) & ([inData.rfDim] ~= 1); %%% Target dimming cannot be 1. After first dimming, RF location can dim at 2nd or 3rd (it can be 0, this means it didn't and isn't going to dim this trial)
                            fprintf('\teliminate %d out of %d trials\n', length(find(~trIdxEvent)), length(inData))
                        case 'Dim3'
                            fprintf('\tFocusing on time period after second dimming, eliminating trials where stimulus in RF or target stimulus dimmed at first or second dimming\n')
                            trIdxEvent = ([inData.targetDim] > 2) & ([inData.rfDim] ~= 1 | [inData.rfDim] ~= 2); %%% Target dimming cannot be 1 or 2. After first dimming, RF location can dim at 3rd
                            fprintf('\teliminate %d out of %d trials\n', length(find(~trIdxEvent)), length(inData))
                            
                    end
                    trIdxEvent = trIdxEvent(:);
                otherwise
                    trIdxEvent = true(length(inData),1);
            end
    end
end

%%% ------------------------------------------------------------------- %%%
function [trIdxMS, dataWarning] = getMicroSaccadeTrialIdx(recInfo, analysisParam, inData, ievent, event)
% find trials in which microsaccades occurred (in between specific times)

    trIdxMS = false(length(inData),1);
    trIdxMS([inData.fixbreak]) = true;% remove trials with fixation break 
    
    dataWarning = false;
    if isfield(analysisParam,'excludeMicrosaccades')
        if analysisParam.excludeMicrosaccades

            tmpTimeEpoch = selectTimeEpoch(inData, analysisParam.event(ievent,:), analysisParam);
            
            readName = fullfile(recInfo.path_read, 'microsaccades.mat');
            readData = load(readName, [event 'Align']);
            
            emptystruct = ~isfield( readData, [event 'Align'] );
            
            if emptystruct
                fprintf('No microsaccade data for event: %s, skipping...\n', event)
                trIdxMS = true(length(inData),1);
                dataWarning = true;
                return
            end
            
            for itrial = 1:length(inData)
                
                if any(isnan(tmpTimeEpoch(itrial,:)))
                    continue
                end
                
                ms = readData.([event 'Align']).ms(itrial).StartTime;
                if any( (ms>=tmpTimeEpoch(itrial,1)) & (ms<=tmpTimeEpoch(itrial,2)) )
                    trIdxMS(itrial) = true;
                end
            end
        end
    end
end

%%% ------------------------------------------------------------------- %%%
function outData = getPupil(pupilData, analysisParam, event)
    % outData = getPupil(pupilData, analysisParam, event)
    % 
    % Get average pupil data within a specified time window
    

    ievent = find(strcmpi(analysisParam.event.Epoch, event));

    switch analysisParam.pupilType
        case 'baseline'

            timeIdx(1) = find(pupilData.([event 'Align']).TimeStamps > analysisParam.event.TimeWindow(ievent,1), 1, 'first');
            timeIdx(2) = find(pupilData.([event 'Align']).TimeStamps > analysisParam.event.TimeWindow(ievent,2), 1, 'first');

            pupilDat = squeeze(nanmean(pupilData.([event 'Align']).Samples(1,:,timeIdx(1):timeIdx(2)), 3));

    end
    
    outData = pupilDat;
    outData = outData(:);
end
    
%%% ------------------------------------------------------------------- %%%
function outData = selectChannels(inData, inDataType, varargin)
    % select trials for specific condition

    allFields = fieldnames(inData);
    allEventIdx = find(contains(allFields, 'Align'));

    % NLX data
    outData = inData;
    
    switch varargin{1}
        case 'area'
            keyboard
            %         chanIdx =
        case 'channels'
            switch islogical(varargin{2})
                case 1
                    chanIdx = varargin{2};
                case 0
                    switch inDataType
                        case 'unit'
                            chanIdx = ismember(inData.channel, varargin{2});
                            
                        otherwise
                            %chanIdx = logical(sum(inData.unitList == varargin{2}, 2));
                            chanIdx = ismember(inData.unitList, varargin{2});
                    end
            end
    end
    
    outData.numUnit = length(find(chanIdx));
    outData.area = outData.area(chanIdx,:);
    outData.unitList = outData.unitList(chanIdx);

    switch inDataType
        case 'unit'
            outData.channel = outData.channel(chanIdx);
            outData.unitClassification = outData.unitClassification(chanIdx,:);
            
            outData.unitClassification = [table(outData.unitList, 'VariableNames', {'ID'}) outData.unitClassification];
    end
    
    
    for ievnt = 1:length(allEventIdx)
        switch inDataType
            case {'NCS','LFP','LFPb','MUAe'}
                outData.(allFields{allEventIdx(ievnt)}).Samples = inData.(allFields{allEventIdx(ievnt)}).Samples(chanIdx,:,:);
            case {'hash','unit','MUA_spont20'}
                outData.(allFields{allEventIdx(ievnt)}) = inData.(allFields{allEventIdx(ievnt)})(chanIdx,:);
        end
    end
end

%%% ------------------------------------------------------------------- %%%
function timeEpoch = selectTimeEpoch(trialdata, epochTable, analysisParam)
    % timeEpoch = selectTimeEpoch(trialdata, epochTable, analysisParam)
    % 
    % Extract time epoch to analyse
    % 
    assert(height(epochTable)==1, 'Can only process a single event')
    
    timeEpoch = NaN(length(trialdata),2);

    for itrial = 1:length(trialdata)
        
        switch epochTable.EventType{1}
            case 'NLX'
                events = trialdata(itrial).NLX_events;
            case 'CTX'
                events = trialdata(itrial).eventArray;
        end
        idx1 = find(events(:,2)==getEventCodes(epochTable.Event{1}, epochTable.EventType{1}),1,'first');
        idx2 = find(events(:,2)==getEventCodes(epochTable.NextEvent{1}, epochTable.EventType{1}),1,'first');
        
        if isempty(idx1)
            continue
        end
        
        % correct for alignment event, convert to seconds
        switch epochTable.EventType{1}
            case 'NLX'
                events(:,1) = (events(:,1) - events(idx1,1)) / 1e6;
            case 'CTX'
                events(:,1) = (events(:,1) - events(idx1,1)) / 1e3;
        end
        
        switch epochTable.Extraction{1}
            case 'full'
                if isempty(idx2)
                    continue
                end
                % get full epoch until the next event, plus any
                % timewindow
                timeEpoch(itrial,1) = events(idx1,1)  + epochTable.TimeWindow(1);
                timeEpoch(itrial,2) = events(idx2,1)  + epochTable.TimeWindow(2);
            case 'fixed'
                % get the exact time window around event
                timeEpoch(itrial,1) = events(idx1,1)    + epochTable.TimeWindow(1);
                timeEpoch(itrial,2) = events(idx1,1)    + epochTable.TimeWindow(2);
        end
        
    end
    
    % test whether the extracted time window is correct
    if any(timeEpoch(:,2) < timeEpoch(:,1))
        error('incorrect epoch extraction')
    end
    
    % test whether the extracted time window is appropriate for
%     % analysisType
%     switch analysisParam.analysisType
%         case {'rate','FF','noiseCorr'}
%             if (length(unique(timeEpoch(:,1)))>1) || (length(unique(timeEpoch(:,2)))>1)
%                 error('%s analysis needs fixed time epoch', analysisParam.analysisType)
%             end
%             timeEpoch = timeEpoch(1,:);
%     end
end

%%% ------------------------------------------------------------------- %%%
function HMM = getHMM(recInfo, analysisParam, event, condition, area)
    % load HMM data, computed in HMM_main.m
    
    % see whether there are specific analysisParam settings for loading HMM    
    [analysisParamHMM] = specifyAnalysisParam(analysisParam.HMMtype, recInfo);
    analysisParamHMM.dataType = 'hash';% force 
    
    if isfield(analysisParam, 'HMMEvent')
        fprintf('\tLoading HMM from different event: %s instead of %s\n', analysisParam.HMMEvent, event)
        event = analysisParam.HMMEvent;
    end
    if isfield(analysisParam, 'HMMbinSize')
        analysisParamHMM.binSize = analysisParam.HMMbinSize;
    end
    
    ievent = find(strcmpi(analysisParamHMM.event.Epoch, event));

    read_path = getAnalysisSaveDir(recInfo, analysisParamHMM, ievent, area, condition );
    
    readName = fullfile(read_path, sprintf('HMMfit_numState_%d.mat', analysisParam.numState));
    
    HMM = load(readName);
    HMM.event = event;
end

%%% ------------------------------------------------------------------- %%%
function [msu, label] = getHMM_transitionTimes(states, timeBins, numState, binSize, HMMtype)
    
    if ~strcmpi(HMMtype,'HMM')
        error('not implemented')
    end
    
    numArea = size(states,2);
    
    %%% what are the maximum number of bins that were used
    maxtwin     = cellfun(@length, timeBins(:, 1), 'UniformOutput', true);
%     bins        = ((0:max(maxtwin)) * binSize) + timeBins{1}(1); % centre bin
    bins        = ((0:max(maxtwin)) * binSize - binSize/2) + timeBins{1}(1); % bin edges
    
    %%% get up and down state switches
    %stateswitch             = cellfun(@(X)(find([0 (diff(X) == 1) | (diff(X) == -1)])), states, 'UniformOutput', false); % when did a stateswitch happen
    %stateswitch_times       = cellfun(@(X)(bins(X)), stateswitch, 'UniformOutput', false); % when did a stateswitch happen
    %stateswitch_label       = cellfun(@(X)( ismember( find([0 (diff(X) == 1) | (diff(X) == -1)]), find([0 (diff(X) == 1)]) ) + 1 ), states, 'UniformOutput', false); % were they down (1) or up (2) switches?
    stateswitch_times_up    = cellfun(@(X)(bins([false (diff(X) == 1)])), states, 'UniformOutput', false); % were they down (1) or up (2) switches?
    stateswitch_times_down  = cellfun(@(X)(bins([false (diff(X) == -1)])), states, 'UniformOutput', false); % were they down (1) or up (2) switches?
    
    %stateswitch_idx_up      = ~cellfun(@isempty, stateswitch_times_up, 'UniformOutput', true); % on which trials was there a state transition
    %stateswitch_idx_down    = ~cellfun(@isempty, stateswitch_times_down, 'UniformOutput', true); % on which trials was there a state transition  
    
    
    % store transition times [area 1 DU, area 1 UD, area 2 DU, area 2 UD, ...];
    msu = [];
    for iarea = 1:numArea
        msu = [msu stateswitch_times_up(:,iarea) stateswitch_times_down(:,iarea)];
    end
    
    % store transition times: DU (area 1, 2...), UD (area 1, 2...)
%     msu = [stateswitch_times_up stateswitch_times_down];
    label = {'DU', 'UD'};
end

    
%%% ------------------------------------------------------------------- %%%
function outData = convertMat2Cell(inData)

    allFields = fieldnames(inData);
    allEventIdx = find(contains(allFields, 'Align'));
    
    outData = inData;
    for ievnt = 1:length(allEventIdx)
        [numchan,numtrials,numtimes] = size(inData.(allFields{allEventIdx(ievnt)}).Samples);
        
        outData.(allFields{allEventIdx(ievnt)}).Samples = cell(numchan, numtrials);
        D = ones(numtrials,1);
        for ichan = 1:numchan
            outData.(allFields{allEventIdx(ievnt)}).Samples(ichan,:) = mat2cell(squeeze(inData.(allFields{allEventIdx(ievnt)}).Samples(ichan,:,:)), D);
        end
       
        outData.(allFields{allEventIdx(ievnt)}).Samples( (numchan+1),:) = repmat( {outData.(allFields{allEventIdx(ievnt)}).TimeStamps}, [1 numtrials]);
        
    end
end

%%% ------------------------------------------------------------------- %%%
function outData = convoluteSpikes(inData, time, sigma)

    allFields = fieldnames(inData);
    allEventIdx = find(contains(allFields, 'Align'));
    
    outData = inData;
    for ievnt = 1:length(allEventIdx)
        [numchan,numtrials] = size(inData.(allFields{allEventIdx(ievnt)}));
        
        outData.(allFields{allEventIdx(ievnt)}) = cell(numchan+1, numtrials);
        for ichan = 1:numchan
            for itrial = 1:numtrials
                times = time(itrial,1):(1/1000):time(itrial,2);
                outData.(allFields{allEventIdx(ievnt)}){ichan,itrial} = spike_convoluteGaussian(inData.(allFields{allEventIdx(ievnt)}){ichan,itrial}, times, sigma);
                outData.(allFields{allEventIdx(ievnt)}){numchan+1,itrial} = times;
            end
        end
       
        outData.SF = 1000;
    end
end


%%% ------------------------------------------------------------------- %%%
function outData = excludeUnits(baseData, inData, analysisParam, timeEpoch)
    
    allFields = fieldnames(inData);
    allEventIdx = find(contains(allFields, 'Align'));
    
    unit_numTrial = 0.1; %[unit_numTrial] proportion of trials used at start and end
    rate_threshold = 1;

    outData = inData;
    for ievnt = 1:length(allEventIdx)
        
        % get units to exclude based on stimulus induced activity
%         msu = baseData.(allFields{allEventIdx(ievnt)});
%         [numUnit,numTrial] = size(msu);
%         
        msu = baseData.StimAlign;
        [numUnit,numTrial] = size(msu);
        timeEpoch = repmat( [0.03 0.2], [numTrial, 1]);
        
        excludeUnit = false(numUnit,1);
        for iunit = 1:numUnit
            % get the number of spikes in the time window per trial
            count = cellfun(@(x,y,z) length(find(x>=y & x<z)), msu(iunit,:)', num2cell(timeEpoch(:,1)), num2cell(timeEpoch(:,2)), 'UniformOutput', true );
            
            % convert to firing rate per trial
            rate = count ./ (timeEpoch(:,2) - timeEpoch(:,1));

            % fit line to rate
            P = polyfit(1:numTrial,rate',1);
            
            % evaluate fit at start and end
            rate_win = polyval(P,[1 numTrial]);
            
            % define exclusion criteria
            numTrial_win = ceil(numTrial*unit_numTrial);
            
            cond1 = any(rate_win<0); % if either end of recording shows 0 spike count continuously, this should show up as a negative spike rate in the linear fit
            cond2 = isempty(find( count(1:numTrial_win) )); % if first trials don't have spikes
            cond3 = isempty(find( count((numTrial - numTrial_win):numTrial) )) ; % if last trials don't have spikes
            %cond4 = rate_mean < rate_threshold; % exclude units with average rate < rate_threshold
            
            if cond1 || cond2 || cond3 %|| cond4
                excludeUnit(iunit) = true;
            end
            
        end
        
        % exclude units
        % fprintf('Excluding %d/%d units\n', length(find(excludeUnit)), numUnit)
        outData.(allFields{allEventIdx(ievnt)}) = outData.(allFields{allEventIdx(ievnt)})(~excludeUnit,:);
        outData.unitList = outData.unitList(~excludeUnit);
        outData.channel = outData.channel(~excludeUnit);
        outData.unitClassification = outData.unitClassification(~excludeUnit,:);
        
        outData.area = outData.area(~excludeUnit,:);
        outData.numUnit = length(find(~excludeUnit));

%         outData.(allFields{allEventIdx(ievnt)}) = msu;

    end
end




