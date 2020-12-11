function pool_data(recordingList, analysisType, readDataFromServer, saveDataToServer, varargin)
% pool_data(recordingList, analysisType, readDataFromServer, saveDataToServer, varargin)
%
% Pool data from specific analysis
% Population data are stored in ./data/population/ANALYSISNAME/EVENTNAME/
%
% Parameters
% ----------
% recordingList : table
%     table with information about recordings, acquired via _getRecordingList 
% analysisType : string
%     analysis to run
% readDataFromServer : boolean
%     boolean indicating whether to read data from server location (see
%     computerSetup.m) 
% saveDataToServer : boolean
%     boolean indicating whether to save data to server location (see
%     computerSetup.m) 
% varargin : cell
%     specific settings, e.g. exampleRecording
% 
% Returns
% -------
% 

fprintf('Pooling data for %s\n', analysisType)

%%% check varargin
if nargin>=4
    varargin2check = {'exampleRecording','HMMtype','dataType'};
    for ivarargin = 1:length(varargin2check)
        idx = find(strcmpi(varargin, varargin2check{ivarargin}));
        if ~isempty(idx)
            eval([ varargin2check{ivarargin} ' = varargin{idx+1}']);
        end
    end
end

% get rid of 'analyse_spike_' in analysisType
analysisType = regexprep(analysisType, 'analyse_spike_', '');

% load analysisParam, overwrite specific settings
[analysisParam] = specifyAnalysisParam(analysisType);

if exist('HMMtype','var')
    analysisParam.HMMtype = HMMtype;
end
if exist('dataType','var')
    analysisParam.dataType = dataType;
end

% set filenames depending on the analysis
switch analysisType
    
    case 'rate'
        loadfilename = 'rate.mat';
        savefilename = sprintf('rate_%s_%s_matchRateAcrossCond%d.mat', analysisParam.dataType, analysisParam.windowType, analysisParam.matchRateAcrossCond);

    case 'FF'
        loadfilename = 'FF.mat';
        savefilename = sprintf('FF_%s_%s.mat', analysisParam.dataType, analysisParam.windowType);

    case 'noiseCorr'
        loadfilename = 'noiseCorr.mat';
        savefilename = loadfilename;

    case 'HMM'
        loadfilename{1} = sprintf('HMMfit_numState_%d.mat', analysisParam.maxNumState);
        loadfilename{2} = sprintf('crossValHMMfit_numState_%d.mat', analysisParam.maxNumState);
%         savefilename = loadfilename{1};
        
        savefilename = sprintf('HMMfit_numState_%d_%s_includeBorder%d_matchRate%d_exclMs%d.mat', ...
            analysisParam.maxNumState, analysisParam.dataType, analysisParam.inclBorder, analysisParam.matchRateAcrossCond, analysisParam.excludeMicrosaccades);
                
    case 'HMMselect'
        loadfilename = sprintf('selectHMM_maxNumState_%d.mat', analysisParam.maxNumState);
        savefilename = sprintf('selectHMM_maxNumState_%d_%s.mat', analysisParam.maxNumState, analysisParam.dataType);
        
    case 'HMM_2area'
        loadfilename = sprintf('HMMfit_numState_%d.mat', analysisParam.maxNumState);
        savefilename = loadfilename;
        
    case 'HMM_crossCorrelation'
        
        switch analysisParam.CC_type
            case 'timeSeries'
                loadfilename = 'crossCorrelation_analog.mat';
            case 'cpt'
                loadfilename = 'HMM_TTA_cpt.mat';
            case 'jitter'
                loadfilename = sprintf('HMM_TTA_jitter%1.3f.mat',analysisParam.jitterBinSize);
        end
        
        savefilename = sprintf('crossCorrelation_%s_%s_exms%d.mat', analysisParam.HMMtype, analysisParam.CC_type, analysisParam.excludeMicrosaccades);
                
    case 'HMM_cvStateBased'
        loadfilename = sprintf('crossValStateBased_numState_%d.mat', analysisParam.maxNumState);
%         savefilename = loadfilename{1};
        
        savefilename = sprintf('crossValStateBased_numState_%d_%s.mat', analysisParam.maxNumState, analysisParam.dataType);
        
    case 'HMM_microsaccades'
        loadfilename{1} = sprintf('HMM_MS_%s.mat', analysisParam.CC_type);
        loadfilename{2} = 'HMM_MS_rate.mat';
        savefilename = sprintf('HMM_MS_%s_%s.mat', analysisParam.HMMtype, analysisParam.CC_type);
        
    case {'HMM_pupil'} % analyses with either HMM or HMM_2area
        loadfilename = 'HMM_corr.mat';
        savefilename = sprintf('HMM_corr_%s.mat', analysisParam.HMMtype);
        
    case {'HMM_rate'} % analyses with either HMM or HMM_2area
        loadfilename = 'HMM_rate.mat';
        savefilename = sprintf('HMM_rate_%s_%s_matchRate%d.mat', analysisParam.HMMtype, analysisParam.dataType, analysisParam.matchRateAcrossCond);
        
    case 'HMM_RT' 
        %         switch analysisParam.HMMtypes
        %             case 'HMM'
        %                 numState = 2;
        %             case 'HMM_2area'
        %                 numState = 4;
        %         end
        %         loadfilename = sprintf('%s_%s_HMMfit_numState_%u.mat', analysisParam.RT2use, analysisParam.RT2transform, numState);
        savefilename = sprintf('HMM_RT.mat');
                                            
    case {'HMM_spectrogram'} % analyses with either HMM or HMM_2area
        loadfilename = 'HMM_spectrogram.mat';
        savefilename = sprintf('HMM_spectrogram_%s.mat', analysisParam.HMMtype);
        
    case 'HMM_transitionTime'
        loadfilename = 'HMM_transitionTimeStats.mat';
        savefilename = sprintf('HMM_transitionTimeStats_%s_binSize%1.3f.mat', analysisParam.HMMtype, analysisParam.binSize);        
        
    case {'HMM_transitionTriggeredAverage'}
        switch analysisParam.CC_type
            case 'cpt'
                loadfilename = sprintf('HMM_TTA_%s.mat', analysisParam.CC_type);
            case 'jitter'
                loadfilename = sprintf('HMM_TTA_%s%1.3f.mat', analysisParam.CC_type, analysisParam.jitterBinSize);
        end
        
        if analysisParam.convoluteSpikes
            dataType = [analysisParam.dataType '_conv'];
        else
            dataType = [analysisParam.dataType];
        end
        savefilename = sprintf('%s_%s.mat', loadfilename(1:end-4), dataType);   
        
    case 'microsaccade_analyse_parameters'
        loadfilename = 'msStat.mat';
        savefilename = sprintf('msStat.mat');
        
    otherwise
        error('unknown analysis')
end

% get recInfo and paths
[recInfo] = getRecordingInfo(recordingList, 'RecIdx', 1);
[recInfo, paths] = set_paths(recInfo,[],'processed',readDataFromServer,saveDataToServer);
path_population = [paths.target paths.data 'population' filesep 'ANALYSISNAME' filesep 'EVENTNAME' filesep];

% get num events
switch analysisType
    case 'HMM_RT'
        nEvents = 1; % combination of events happens in HMM_RT.m
    otherwise
        nEvents = height(analysisParam.event);
end

%% loop over recordings once to get the recordingList, add RF, find maximum
% numArea, etc.
nareaTable = 0;
fprintf('\tammending recordingList\n')
nArea = NaN(height(recordingList),1);% keep track of how many areas we recorded from
for irec = 1:height(recordingList)
    % get recInfo and paths
    [recInfo] = getRecordingInfo(recordingList, 'RecIdx', irec);
    [recInfo] = set_paths(recInfo,[],'processed',readDataFromServer,saveDataToServer);
    
    % get areaInfo
    [areaInfo, areaInfoWarning] = getAreaInfo(recInfo, recInfo.path_read);% areaInfo
    [recInfo, areaInfo]         = getRFinfo(recInfo, areaInfo, recInfo.path_read);
    
    if areaInfoWarning
        warning('Could not load areaInfo')
        keyboard
    end
    nArea(irec) = length(areaInfo);
    
    % update recordingList with RF
    areaLabel = cell(1, length(areaInfo));
    for iarea = 1:nArea(irec)
        recordingList.(['RF_' areaInfo(iarea).name])(irec,:) = areaInfo(iarea).centerRF(:)';
        areaLabel{iarea} = areaInfo(iarea).name;
    end
    recordingList.RF_distance(irec) = recInfo.distanceRF;
    recordingList.RF_overlap(irec) = recInfo.overlapRF;
    
    if nArea(irec) > nareaTable
        nareaTable = nArea(irec);
        areaTable = array2table(1:length(areaInfo),'VariableNames',areaLabel);
    end        
end

%% make population file
% loop over events. We make a separate file for each event.

for ievent = 1:nEvents
    clear dataMat
    
    eventLabel = analysisParam.event.Epoch{ievent};
    fprintf('\tPooling data for event: %s\n', eventLabel)
    
    % loop over recordings
    rec2use = 0; %keep track of included recordings
    includedRecordings = false(height(recordingList),1);
    for irec = 1:height(recordingList)
        
        dispProgress(irec,height(recordingList),height(recordingList))
        
        % get recInfo and paths
        [recInfo] = getRecordingInfo(recordingList, 'RecIdx', irec);
        [recInfo] = set_paths(recInfo,[],'processed',readDataFromServer,saveDataToServer);
        
        % get areaInfo
        [areaInfo, areaInfoWarning] = getAreaInfo(recInfo, recInfo.path_read);% areaInfo
        [recInfo, areaInfo]         = getRFinfo(recInfo, areaInfo, recInfo.path_read);
        
        if areaInfoWarning
            warning('Could not load areaInfo')
            keyboard
        end
        
        % update analysisParam with recording specific info
        [analysisParam] = specifyAnalysisParam(analysisType, recInfo);
        
        if exist('HMMtype','var')
            analysisParam.HMMtype = HMMtype;
        end
        if exist('dataType','var')
            analysisParam.dataType = dataType;
        end
        
        % test whether we recorded from multiple areas, only for specified analyses
        switch analysisType
            case {'HMM_2area','HMM_crossCorrelation','HMM_crossCorrelationPupil'}
                if (nArea(irec)~=2)
                    continue
                end
            case 'HMM_spectrogram'
                switch analysisParam.HMMtype
                    case 'HMM_2area'
                        if (nArea(irec)~=2)
                            continue
                        end
                end
        end
        
        % check if pupil data is present for this recording
        switch analysisParam.analysisType
            case {'HMM_pupil','HMM_crossCorrelationPupil'}
                if ~exist( fullfile(recInfo.path_read,'pupil.mat'), 'file' )
                    continue
                end
        end
        
        % get some info on recording depth
        numChanArea = length(areaInfo(1).ChanIdx.all);
        allLayerAssignment = reshape([areaInfo.layerAssignment], [nArea(irec)*numChanArea, 1]);
        allDepth = reshape([areaInfo.Coordinates], [3*numChanArea, nArea(irec)]);
        allDepth = reshape(allDepth(1:numChanArea,:), [nArea(irec)*numChanArea, 1]);
        
        if nArea(irec)==1 && strcmpi(areaInfo(1).name,'V4')
            allLayerAssignment = [cell(numChanArea,1) ; allLayerAssignment];
            allDepth = [NaN(numChanArea,1) ; allDepth];
        end
        
        % start counting
        rec2use = rec2use+1;
        includedRecordings(irec) = true;

        % get number of conditions
        if (rec2use==1)
            [condLabel, numCond, ~] = getConditions(recInfo, analysisParam, analysisParam.event.ConcatenateCond(ievent));
            dataMat.condNames = condLabel;
        end
        
        % loop conditions
        for icond = 1:numCond
            clear tmpData
            
            %fprintf('\t\tPooling data for condition: %s (%d/%d)\n', condLabel{icond}, icond, numCond)

            switch analysisType
                
                    %%% =============================================================================================================================================================== %%%
                case {'rate','FF'}
                    
                    areas = [];
                    for iarea = 1:nArea(irec)
                        if iarea==1
                            areas = areaInfo(iarea).name;
                        else
                            areas = sprintf('%s_%s',areas,areaInfo(iarea).name);
                        end
                    end
                    
                    read_path = getAnalysisSaveDir(recInfo, analysisParam, ievent, areas, condLabel{icond} );
                    
                    loadfilename = sprintf('%s.mat',analysisParam.analysisType);
                    readfilename = [read_path loadfilename];
                    
                    if ~exist(readfilename, 'file')
                        warning('file does not exist: %s', readfilename)
                        continue
                    end
                    tmpData = load(readfilename);
                    
                    if (irec==1) && (icond==1)
                        % Recording, Event, Cond, area
                        dataMat.readfilename = cell(height(recordingList), numCond);
                        
                        dataMat.dataTable = [];
                        
                        dataMat.time = tmpData.time;
                        dataMat.dataType = tmpData.dataType;
                        
                        try
                            dataMat.window = tmpData.window;
                        catch
                        end
                    end
                    dataMat.readfilename{irec, icond} = readfilename;
                                        
                    % get layer info
                    if istable(tmpData.msUnitList)
                        unit = (tmpData.msUnitList.Channel);
                    else
                        unit = (tmpData.msUnitList);
                    end
                    chanDepth = [allDepth(unit)];
                    chanLayer = [allLayerAssignment(unit)];
                    
                    % make table with unit/channel info
                    if istable(tmpData.msUnitList)
                        unitTable = [tmpData.msUnitList table(unit, chanDepth, chanLayer, 'VariableNames', {'channel', 'channelDepth', 'chanLayer'})];
                    else
                        unitTable = table(unit, chanDepth, chanLayer, 'VariableNames', {'channel', 'channelDepth', 'chanLayer'});                        
                    end
                    
                    numChan = length(find(unit));
                    
                    % find channels/units in areas and make areaIdx
                    areaIdx = zeros(numChan,1);
                    for iarea = 1:nArea(irec)
                        chanIdx = strcmpi(string(tmpData.areaList), areaInfo(iarea).name);

                        areaIdx = areaIdx + chanIdx * iarea;
                    end
                    
                    switch analysisType
                        case 'rate'
                            tmpTable = ([...
                                table(repmat(irec, [numChan 1]), 'VariableNames', {'Recording'}) , ...
                                table(repmat(icond, [numChan 1]), 'VariableNames', {'Condition'}) , ...
                                table(repmat(condLabel(icond), [numChan 1]), 'VariableNames', {'ConditionNames'}) , ...
                                table(areaIdx, 'VariableNames', {'area_idx'}) , ...
                                table(tmpData.areaList, 'VariableNames', {'area'}) , ...
                                unitTable, ...
                                table(tmpData.rate, 'VariableNames', {'rate'}) , ...
                                table(tmpData.sRate, 'VariableNames', {'sRate'})
                                
                                ]);
                        case 'FF'
                            tmpTable = ([...
                                table(repmat(irec, [numChan 1]), 'VariableNames', {'Recording'}) , ...
                                table(repmat(icond, [numChan 1]), 'VariableNames', {'Condition'}) , ...
                                table(repmat(condLabel(icond), [numChan 1]), 'VariableNames', {'ConditionNames'}) , ...
                                table(areaIdx, 'VariableNames', {'area_idx'}) , ...
                                table(tmpData.areaList, 'VariableNames', {'area'}) , ...
                                unitTable, ...
                                table(tmpData.FF, 'VariableNames', {'FF'}) , ...
                                table(tmpData.sFF, 'VariableNames', {'sFF'}), ...
                                table(tmpData.meanCount, 'VariableNames', {'meanCount'}) , ...
                                table(tmpData.varCount, 'VariableNames', {'varCount'})
                            
                                ]);
                            
                    end
                    
                    try
                    dataMat.dataTable = [...
                        dataMat.dataTable ; ...
                        tmpTable];
                    catch
                        keyboard
                    end
                    %
            
                    
                    %%% =============================================================================================================================================================== %%%
                case 'noiseCorr'
                    
                    areas = [];
                    for iarea = 1:nArea(irec)
                        if iarea==1
                            areas = areaInfo(iarea).name;
                        else
                            areas = sprintf('%s_%s',areas,areaInfo(iarea).name);
                        end
                    end
                    
                    read_path = getAnalysisSaveDir(recInfo, analysisParam, ievent, areas, condLabel{icond} );

                    loadfilename = sprintf('%s.mat',analysisParam.analysisType);
                    readfilename = [read_path loadfilename];
                    
                    if ~exist(readfilename, 'file')
                        warning('file does not exist: %s', readfilename)
                        continue
                    end
                    tmpData = load(readfilename);
                    
                    if (rec2use==1) && (icond==1)
                        % Recording, Event, Cond, area
                        dataMat.readfilename = cell(height(recordingList), numCond);
                        % Rec, Cond, Area
                        dataMat.noiseCorr = [];
                        
                        
                        
                        dataMat.time = tmpData.time;
                        dataMat.window = tmpData.window;
                    end
                    dataMat.readfilename{irec, icond} = readfilename;
                                        
                    for iarea = 1:nArea(irec)
                        for iarea_in = 1:nArea(irec)
                            
                            areaCombination = {areaInfo(iarea).name, areaInfo(iarea_in).name};
                            
                            % find channels/units in areas
                            unit = areaInfo(iarea).ChanIdx.all;
                            chan_in = areaInfo(iarea_in).ChanIdx.all;
                            
                            % find indices of channel/unit combination
                            cond1 = ismember(tmpData.channelCombinations(:,1), unit);
                            cond2 = ismember(tmpData.channelCombinations(:,2), chan_in);
                            
                            combIdx = cond1 & cond2;
                            chanComb = squeeze(tmpData.channelCombinations(combIdx,:));      
                            
                            % get layer info
                            chanDepth = [allDepth(chanComb(:,1)) allDepth(chanComb(:,2))];
                            chanLayer = [allLayerAssignment(chanComb(:,1)) allLayerAssignment(chanComb(:,2))];
                            
                            % prep some variables
                            meancount = squeeze(tmpData.meanCount(combIdx,:,:,:));
                            varcount = squeeze(tmpData.varCount(combIdx,:,:,:));
                            
                            if size(meancount,1) ~= length(find(combIdx))% necessary to transpose if n==1
                                meancount = meancount';
                                varcount = varcount';
                            end
                            
                            numComb = length(find(combIdx));
                            tmpTable = ([...
                                table(repmat(irec, [numComb 1]), 'VariableNames', {'Recording'}) , ...
                                table(repmat(icond, [numComb 1]), 'VariableNames', {'Condition'}) , ...
                                table(repmat(condLabel(icond), [numComb 1]), 'VariableNames', {'ConditionNames'}) , ...
                                table(repmat([iarea iarea_in], [numComb 1]), 'VariableNames', {'areaCombination'}) , ...
                                table(repmat(areaCombination, [numComb 1]), 'VariableNames', {'areaCombinationLabel'}) , ...
                                table(chanComb, 'VariableNames', {'channelCombination'}) , ...
                                table(chanDepth, 'VariableNames', {'channelDepth'}) , ...
                                table(chanLayer, 'VariableNames', {'chanLayer'}) , ...
                                table(squeeze(tmpData.noiseCorr(combIdx,:,:)), 'VariableNames', {'noiseCorr'}) , ...
                                table(squeeze(tmpData.sNoiseCorr(combIdx,:,:)), 'VariableNames', {'sNoiseCorr'}) , ...
                                table(squeeze(tmpData.noiseCorrShuffle(combIdx,:,:)), 'VariableNames', {'noiseCorrShuffle'}) , ...
                                table(squeeze(tmpData.sNoiseCorrShuffle(combIdx,:,:)), 'VariableNames', {'sNoiseCorrShuffle'}), ...
                                table(meancount, 'VariableNames', {'meanCount'}), ...
                                table(varcount, 'VariableNames', {'varCount'})
                                ]);

                            dataMat.noiseCorr = [...
                                dataMat.noiseCorr ; ...
                                tmpTable];
                            
                            
                        end
                    end
                    
                    
                    %%% =============================================================================================================================================================== %%%
                case {'HMM'}
                    

                    %%% HMM individually for different areas
                    for iarea = 1:nArea(irec)
                        
                        areaLabel = areaInfo(iarea).name;
                        areaIdx = getAreaIdx(areaTable, areaLabel);
                        
                        %%% -----------------------------------------------
                        %%% first load HMM fit data
                        read_path = getAnalysisSaveDir(recInfo, analysisParam, ievent, areaLabel, condLabel{icond} );
                        readfilename = [read_path loadfilename{1}];
                        if ~exist(readfilename, 'file')
                            warning('file does not exist: %s', readfilename)
                            continue
                        end
                        tmpData = load(readfilename);
                        
                        if (rec2use==1) && (icond==1) && (iarea==1)
                            % Rec, Cond, Area, Pt, Pt+1
                            dataMat.readfilename = cell(height(recordingList), numCond, max(nArea), 2); % filename for HMM fit and cv

                            dataMat.estTrans = NaN(height(recordingList), numCond, max(nArea), analysisParam.maxNumState, analysisParam.maxNumState);
                            dataMat.residTimePDF = NaN(height(recordingList), numCond, max(nArea), size(tmpData.residTimePDF{1},1), size(tmpData.residTimePDF{1},2), analysisParam.maxNumState);

                            % transformed data
                            dataMat.meanEpochDuration = NaN(height(recordingList), numCond, max(nArea), analysisParam.maxNumState);
                            dataMat.meanFractTimeState = NaN(height(recordingList), numCond, max(nArea), analysisParam.maxNumState);
                            
                            % table
                            dataMat.cv = [];
                        end
                        
                        dataMat.readfilename{irec, icond, areaIdx, 1} = readfilename;

                        if any(isnan(tmpData.estTrans(:)))
                            error('NaN values in transition matrix')
                        end
                        dataMat.estTrans(irec, icond, areaIdx, :, :) = tmpData.estTrans;
                        
                        %%% get epoch duration and fraction of time spent in
                        %%% states
                        
                        %%% get the time epochs for each trial for each state
                        [timeEpochState, totalDurationState, numTransition] = extractUpDownTimeEpoch( tmpData.states, tmpData.timeBins, tmpData.numState, tmpData.binSize, analysisParam.inclBorder);
                        % timeEpochState: raw start and end times of states per trial, difference will give epoch duration
                        % totalDurationState: total time spent in a state per trial
                        % numTransition: number of state transitions per trial
                        
                        [dataMat.meanEpochDuration(irec, icond, areaIdx, :), dataMat.meanFractTimeState(irec, icond, areaIdx, :), residTimePDF] = HMM_extractEpochDuration(timeEpochState, tmpData.numState);
                        
                        if 0
                            % residTimePDF that is computed in
                            % hmmAnalysisCrossVal does not have the same
                            % x-axis across recordings
                            for istate=1:analysisParam.maxNumState
                                dataMat.residTimePDF(ipen, icond, areaIdx, :, :, istate) = tmpData.residTimePDF{istate};
                            end
                        else
                            if (rec2use==1) && (icond==1) && (iarea==1)
                                % reset time size
                                dataMat.residTimePDF = NaN(height(recordingList), numCond, max(nArea), size(residTimePDF,1), size(tmpData.residTimePDF{1},2), analysisParam.maxNumState);
                            end
                            dataMat.residTimePDF(irec, icond, areaIdx, :, :, :) = residTimePDF;
                        end
                        
                        if (rec2use==1)
                            dataMat.binSize = tmpData.binSize;
                            dataMat.numState = tmpData.numState;
                        end


                        %%% -----------------------------------------------
                        %%% now load HMM cv data
                        read_path = getAnalysisSaveDir(recInfo, analysisParam, ievent, areaLabel, condLabel{icond} );
                        readfilename = [read_path loadfilename{2}];
                        if ~exist(readfilename, 'file')
                            warning('file does not exist: %s', readfilename)
                            continue
                        end
                        tmpData = load(readfilename);
                        
                        dataMat.readfilename{irec, icond, iarea, 2} = readfilename;
                        
                        numChan = size(tmpData.cvError,1);
                        
                        chans = areaInfo(iarea).ChanIdx.all(areaInfo(iarea).ChanIdx.include);
                        
                        assert(numChan == length(chans), 'Channel numbers do not match')
                        
                        tmpTable = ([...
                            table(repmat(string(recInfo.Subject{1}), [numChan 1]), 'VariableNames', {'Subject'}) , ...
                            table(repmat(recInfo.Date, [numChan 1]), 'VariableNames', {'Date'}) , ...
                            table(repmat(irec, [numChan 1]), 'VariableNames', {'Recording'}) , ...
                            table(repmat(icond, [numChan 1]), 'VariableNames', {'Condition'}) , ...
                            table(repmat(condLabel(icond), [numChan 1]), 'VariableNames', {'ConditionNames'}) , ...
                            table(repmat(areaInfo(iarea).name, [numChan 1]), 'VariableNames', {'area'}) , ...
                            table(repmat(iarea, [numChan 1]), 'VariableNames', {'area_idx'}) , ...
                            table(chans(:), 'VariableNames', {'Channel'}) , ...
                            table(squeeze(mean(tmpData.cvError,3)), 'VariableNames', {'cvError'}) , ...
                            table(squeeze(mean(tmpData.looError,3)), 'VariableNames', {'looError'}) , ...
                            table(squeeze(mean(tmpData.veError,3)), 'VariableNames', {'veError'}) , ...
                            table(squeeze(mean(tmpData.mrError,3)), 'VariableNames', {'mrError'}) , ...
                            table(squeeze(mean(tmpData.R2,3)), 'VariableNames', {'R2'}) , ...
                            table(squeeze(mean(tmpData.veR2,3)), 'VariableNames', {'veR2'}) , ...
                            table(squeeze(mean(tmpData.looR2,3)), 'VariableNames', {'looR2'}) , ...
                            table(squeeze(tmpData.meanRate), 'VariableNames', {'meanRate'}) , ...
                            ]);
                        
                        dataMat.cv = [...
                            dataMat.cv ; ...
                            tmpTable];
                                       
                        if (rec2use==1)
                            dataMat.timeWindow = tmpData.timeWindow;
                        end
                        
                    end
                    
                    %%% =============================================================================================================================================================== %%%
                case 'HMMselect'
                                        
                    %%% HMMselect individually for different areas
                    for iarea = 1:nArea(irec)
                        
                        areaLabel = areaInfo(iarea).name;
                        areaIdx = getAreaIdx(areaTable, areaLabel);

                        read_path = getAnalysisSaveDir(recInfo, analysisParam, ievent, areaLabel, condLabel{icond} );
                        readfilename = [read_path loadfilename];
                        if ~exist(readfilename, 'file')
                            warning('file does not exist: %s', readfilename)
                            continue
                        end
                        tmpData = load(readfilename);
                        
                        if (rec2use==1) && (icond==1) && (iarea==1)
                            % Rec, Cond, area, ...
                            dataMat.readfilename = cell(height(recordingList), numCond, max(nArea));
    
                            dataMat.cvError = NaN(height(recordingList), numCond, max(nArea), analysisParam.maxNumState, length(tmpData.timeWindow));
                            dataMat.looError = NaN(height(recordingList), numCond, max(nArea), analysisParam.maxNumState, length(tmpData.timeWindow));
                            dataMat.veError = NaN(height(recordingList), numCond, max(nArea), analysisParam.maxNumState, length(tmpData.timeWindow));
                            
                            dataMat.timeWindow = tmpData.timeWindow;
                        end
                        
                        dataMat.readfilename{irec, icond, areaIdx} = readfilename;

                        % sum cross validation error over channels and cross validations
                        dataMat.cvError(irec, icond, areaIdx, :, :) = squeeze(sum(sum(tmpData.cvError,1),4));
                        dataMat.looError(irec, icond, areaIdx, :, :) = squeeze(sum(sum(tmpData.looError,1),4));
                        dataMat.veError(irec, icond, areaIdx, :, :) = squeeze(sum(sum(tmpData.veError,1),4));

                    end
                                        
                    %%% =============================================================================================================================================================== %%%
                case {'HMM_2area'}
                    
                    areaLabel = sprintf('%s_%s',areaInfo(1).name,areaInfo(2).name);
                    
                    read_path = getAnalysisSaveDir(recInfo, analysisParam, ievent, areaLabel, condLabel{icond} );
                    readfilename = [read_path loadfilename];
                    if ~exist(readfilename, 'file')
                        warning('file does not exist: %s', readfilename)
                        continue
                    end                    
                    tmpData = load(readfilename);
                    
                    if (rec2use==1) && (icond==1)
                        % Recording, Event, Cond, Pt, Pt+1
                        dataMat.readfilename = cell(height(recordingList), numCond);
                        
                        dataMat.estTrans = NaN(height(recordingList), numCond, analysisParam.maxNumState, analysisParam.maxNumState);
                        dataMat.pen = NaN(height(recordingList), numCond);
                        dataMat.cond = NaN(height(recordingList), numCond);
                        
                        % transformed data
                        dataMat.meanEpochDuration = NaN(height(recordingList), numCond, analysisParam.maxNumState);
                        dataMat.meanFractTimeState = NaN(height(recordingList), numCond, analysisParam.maxNumState);
                    end
                    
                    dataMat.readfilename{irec, icond} = readfilename;

                    
                    if any(isnan(tmpData.estTrans(:)))
                        error('NaN values in transition matrix')
                    end
                    dataMat.estTrans(irec, icond, :, :) = tmpData.estTrans;
                    
                    %%% get epoch duration and fraction of time spent in
                    %%% states
                    
                    %%% get the time epochs for each trial for each state
                    [timeEpochState, totalDurationState, numTransition] = extractUpDownTimeEpoch( tmpData.states, tmpData.timeBins, tmpData.numState, tmpData.binSize, analysisParam.inclBorder);
                    % timeEpochState: raw start and end times of states per trial, difference will give epoch duration
                    % totalDurationState: total time spent in a state per trial
                    % numTransition: number of state transitions per trial
                    
                    [dataMat.meanEpochDuration(irec, icond, :), dataMat.meanFractTimeState(irec, icond, :)] = HMM_extractEpochDuration(timeEpochState, tmpData.numState);
                    
                    if (rec2use==1)
                        dataMat.binSize = tmpData.binSize;
                        dataMat.numState = tmpData.numState;
                    end
                                        
                    %%% =============================================================================================================================================================== %%%
                case 'HMM_crossCorrelation'
                    
                    areaLabel = sprintf('%s_%s',areaInfo(1).name,areaInfo(2).name);
                    
                    read_path = getAnalysisSaveDir(recInfo, analysisParam, ievent, areaLabel, condLabel{icond} );
                    readfilename = [read_path loadfilename];
                    
                    if ~exist(readfilename,'file')
                        fprintf('%s does not exist, skipping...\n', [read_path loadfilename])
                        continue
                    end
                    tmpData = load(readfilename);
                                                            
                    switch analysisParam.CC_type
                        case 'timeSeries'
                            
                            if (rec2use==1) && (icond==1)
                                % Recording, Cond, time
                                dataMat.readfilename = cell(height(recordingList), numCond);

                                dataMat.crossCorr = NaN(height(recordingList), numCond, length(tmpData.timeBinsLag));
                                dataMat.crossCorrShuffle = NaN(height(recordingList), numCond, length(tmpData.timeBinsLag));
                                dataMat.AUC = NaN(height(recordingList), numCond);
                                dataMat.AUC_side = NaN(height(recordingList), numCond, 2);
                            end
                            
                            dataMat.readfilename{irec, icond} = readfilename;

                            dataMat.crossCorr(irec, icond, :) = tmpData.crossCorr;
                            dataMat.crossCorrShuffle(irec, icond, :) = tmpData.crossCorr_shift;
                            dataMat.AUC(irec, icond) = tmpData.crossCorr_AU;
                            dataMat.AUC_side(irec, icond,:) = tmpData.crossCorr_AU_side;
                                                        
                            if (rec2use==1)
                                dataMat.timeBinsLag = tmpData.timeBinsLag;
                            end
                            
                        case {'cpt','jitter'}
                            
                            for istate = 1:analysisParam.numState
                                chanIdx = contains((tmpData.labelCombinations), [tmpData.label{istate} '_' tmpData.label{istate}]);
                                
                                if (rec2use==1) && (icond==1) && (istate==1)
                                    % Recording, Event, Cond, areaHMM, areaNLX, state, time
                                    dataMat.readfilename = cell(height(recordingList), numCond);
    
                                    dataMat.crossCorr = NaN(height(recordingList), numCond, analysisParam.numState, length(tmpData.timeBinsLag));
                                    dataMat.meanRate = NaN(height(recordingList), numCond, max(nArea));
                                    
                                    dataMat.pen = NaN(height(recordingList), numCond, analysisParam.numState);
                                    dataMat.cond = NaN(height(recordingList), numCond, analysisParam.numState);
                                    dataMat.state = NaN(height(recordingList), numCond, analysisParam.numState);
                                end
                                
                                switch analysisParam.CC_type
                                    case 'cpt'
                                        dataMat.crossCorr(irec, icond, istate, :) = squeeze(mean(tmpData.crossCorr(chanIdx,:),1));
                                        %                                 dataMat.meanRate(ipen, ievnt, icond, iareaHMM, iareaNLX) = squeeze(mean(tmpData.allTA.meanrate(chanIdx),1));
                                    case 'jitter'
                                        dataMat.crossCorr(irec, icond, istate, :) = squeeze( sum(tmpData.crossCorr(chanIdx,:),1)/sqrt(length(find(chanIdx))) );
                                end
                                
                            end
                            dataMat.readfilename{irec, icond} = readfilename;

                            if (rec2use==1)
                                dataMat.timeBinsLag = tmpData.timeBinsLag;
                                dataMat.label = tmpData.label;
                            end
                    end
                    
                    %%% =============================================================================================================================================================== %%%                    
                    
                case 'HMM_cvStateBased'                    
                    
                    %%% HMM individually for different areas
                    for iarea = 1:nArea(irec)
                        
                        areaLabel = areaInfo(iarea).name;
                        areaIdx = getAreaIdx(areaTable, areaLabel);
                        
                        %%% -----------------------------------------------
                        %%% first load HMM fit data
                        read_path = getAnalysisSaveDir(recInfo, analysisParam, ievent, areaLabel, condLabel{icond} );
                        readfilename = [read_path loadfilename];
                        if ~exist(readfilename, 'file')
                            warning('file does not exist: %s', readfilename)
                            continue
                        end
                        tmpData = load(readfilename);
                        
                        if (rec2use==1) && (icond==1) && (iarea==1)
                            % Rec, Cond, Area, Pt, Pt+1
                            dataMat.readfilename = cell(height(recordingList), numCond, max(nArea)); % filename for HMM  cv
                            
                            % table
                            dataMat.cv = [];

                            dataMat.binSize = analysisParam.binSize;
                            dataMat.numState = tmpData.numState;
                        end
                        
                        dataMat.readfilename{irec, icond, areaIdx} = readfilename;
                        
                        numChan = size(tmpData.cvError,1);
                        
                        tmpTable = ([...
                            table(repmat(irec, [numChan 1]), 'VariableNames', {'Recording'}) , ...
                            table(repmat(icond, [numChan 1]), 'VariableNames', {'Condition'}) , ...
                            table(repmat(condLabel(icond), [numChan 1]), 'VariableNames', {'ConditionNames'}) , ...
                            table(repmat(string(areaInfo(iarea).name), [numChan 1]), 'VariableNames', {'area'}) , ...
                            table(repmat(iarea, [numChan 1]), 'VariableNames', {'area_idx'}) , ...
                            table(string(tmpData.dataType), 'VariableNames', {'dataType'}) , ...
                            table(squeeze(mean(tmpData.cvError,3)), 'VariableNames', {'cvError'}) , ...
                            table(squeeze(mean(tmpData.veError,3)), 'VariableNames', {'veError'}) , ...
                            table(squeeze(mean(tmpData.mrError,3)), 'VariableNames', {'mrError'}) , ...
                            table(squeeze(mean(tmpData.R2,3)), 'VariableNames', {'R2'}) , ...
                            table(squeeze(mean(tmpData.veR2,3)), 'VariableNames', {'veR2'}) , ...
                            table(squeeze(tmpData.meanRate), 'VariableNames', {'meanRate'}) , ...
                            ]);
                        
                        dataMat.cv = [...
                            dataMat.cv ; ...
                            tmpTable];
                                       
                        if (rec2use==1)
                            dataMat.timeWindow = tmpData.timeWindow;
                        end
                        
                    end

                    
                    %%% =============================================================================================================================================================== %%%
                case 'HMM_microsaccades'
                    
                    switch analysisParam.HMMtype
                        case 'HMM'
                            nAreaHMM = nArea(irec);
                        case 'HMM_2area'
                            if nArea(irec)<2
                                continue
                            end
                            nAreaHMM = 1;
                    end
                    
                    %%% HMM individually for different areas
                    for iareaHMM = 1:nAreaHMM
                        
                        switch analysisParam.HMMtype
                            case 'HMM'
                                areaLabel = areaInfo(iareaHMM).name;
                                areaIdx = getAreaIdx(areaTable, areaLabel);

                            case 'HMM_2area'
                                areaLabel = sprintf('%s_%s',areaInfo(1).name,areaInfo(2).name);
                                areaIdx = 1;
                                
                        end
                          
                        % Process cross correlation
                        % --------------------------
                        read_path = getAnalysisSaveDir(recInfo, analysisParam, ievent, areaLabel, condLabel{icond} );
                        readfilename = [read_path loadfilename{1}];
                        
                        if ~exist(readfilename, 'file')
                            warning('file does not exist: %s', readfilename)
                            continue
                        end
                        tmpData = load(readfilename);
                        
                        if (rec2use==1) && (icond==1) && (iareaHMM==1)
                            % Recording, Cond, areaHMM
                            dataMat.readfilename = cell(height(recordingList), numCond, max(nArea), 2);

                            dataMat.crossCorr = NaN(height(recordingList), numCond, nAreaHMM, analysisParam.numState, length(tmpData.timeBinsLag));
                            dataMat.meanRate = NaN(height(recordingList), numCond, nAreaHMM, analysisParam.numState);
                            dataMat.meanRateTrigger = NaN(height(recordingList), numCond, nAreaHMM, 1);
                            
                            dataMat.timeBinsLag = tmpData.timeBinsLag;
                            dataMat.label = tmpData.labelCombinations;
                        end
                        
                        dataMat.readfilename{irec, icond, iareaHMM, 1} = readfilename;

                        dataMat.crossCorr(irec, icond, areaIdx, :, :) = tmpData.crossCorr;
                        dataMat.meanRate(irec, icond, areaIdx, :) = tmpData.meanRateDU;
                        
                        switch analysisParam.analysisType
                            case 'HMM_microsaccades'
                                dataMat.meanRateTrigger(irec, icond, areaIdx, :) = tmpData.meanRateMS;
                        end
                        
                        % Process microsaccade theta-HMM
                        % --------------------------                        
                        readfilename = [read_path loadfilename{2}];
                        
                        if ~exist(readfilename, 'file')
                            warning('file does not exist: %s', readfilename)
                            continue
                        end
                        tmpData = load(readfilename);

                        
                        if (rec2use==1) && (icond==1) && (iareaHMM==1)
                            % Recording, Cond, areaHMM, state
                            dataMat.msTable = [];
                        end
                        
                        dataMat.readfilename{irec, icond, iareaHMM, 2} = readfilename;

                        numElements = height( tmpData.msTable );
                        
                        tmpTable = ([...
                            table(repmat(string(recInfo.Subject{1}), [numElements 1]), 'VariableNames', {'Subject'}) , ...
                            table(repmat(recInfo.Date, [numElements 1]), 'VariableNames', {'Date'}) , ...
                            table(repmat(irec, [numElements 1]), 'VariableNames', {'Recording'}) , ...
                            table(repmat(icond, [numElements 1]), 'VariableNames', {'Condition'}) , ...
                            table(repmat(condLabel(icond), [numElements 1]), 'VariableNames', {'ConditionNames'}) , ...
                            table(repmat(areaInfo(iareaHMM).name, [numElements 1]), 'VariableNames', {'area'}) , ...
                            table(repmat(iareaHMM, [numElements 1]), 'VariableNames', {'area_idx'}) , ...
                            tmpData.msTable ...
                            ]);
                        
                        dataMat.msTable = [...
                            dataMat.msTable ; ...
                            tmpTable];
                        
                    end
                    
                    
                    %%% =============================================================================================================================================================== %%%
                case {'HMM_pupil'}
                    switch analysisParam.HMMtype
                        case 'HMM'
                            nAreaHMM = nArea(irec);
                        case 'HMM_2area'
                            if nArea(irec)<2
                                continue
                            end
                            nAreaHMM = 1;
                    end
                    
                    %%% HMM individually for different areas
                    for iareaHMM = 1:nAreaHMM
                        
                        switch analysisParam.HMMtype
                            case 'HMM'
                                areaLabel               = areaInfo(iareaHMM).name;
                                areaIdx                 = getAreaIdx(areaTable, areaLabel);
                                analysisParam.numState  = 2;
                            case 'HMM_2area'
                                areaLabel               = sprintf('%s_%s',areaInfo(1).name,areaInfo(2).name);
                                areaIdx                 = 1;
                                analysisParam.numState  = 4;
                        end
                        
                        read_path = getAnalysisSaveDir(recInfo, analysisParam, ievent, areaLabel, condLabel{icond} );
                        readfilename = [read_path loadfilename];
                        if ~exist(readfilename, 'file')
                            warning('file does not exist: %s', readfilename)
                            continue
                        end
                        tmpData = load(readfilename);
                        
                        if (rec2use==1) && (icond==1) && (iareaHMM==1)
                            % Recording, Cond, areaHMM, state
                            dataMat.readfilename = cell(height(recordingList), numCond, max(nArea));
                            
                            dataMat.R_epochDuration = NaN(height(recordingList), numCond, nAreaHMM, analysisParam.numState);
                            dataMat.R_epochDurationShuffle = NaN(height(recordingList), numCond, nAreaHMM, analysisParam.numState);
                            
                            dataMat.R_fractimeState = NaN(height(recordingList), numCond, nAreaHMM, analysisParam.numState);
                            dataMat.R_fractimeStateShuffle = NaN(height(recordingList), numCond, nAreaHMM, analysisParam.numState);
                            
                            dataMat.meanData = NaN(height(recordingList), numCond);
                            
                        end
                        dataMat.readfilename{irec, icond, iareaHMM} = readfilename;
                        
                        dataMat.R_epochDuration(irec, icond, areaIdx, :) = tmpData.R_epochDuration;
                        dataMat.R_epochDurationShuffle(irec, icond, areaIdx, :) = tmpData.R_epochDurationShuffle;
                        
                        dataMat.R_fractimeState(irec, icond, areaIdx, :) = tmpData.R_fractimeState;
                        dataMat.R_fractimeStateShuffle(irec, icond, areaIdx, :) = tmpData.R_fractimeStateShuffle;
                        
                        dataMat.meanData(irec, icond) = mean(tmpData.msu{1});
                        
%                         %%% example subject/recording
%                         if exist('exampleRecording','var') && strcmpi(exampleRecording.Subject, recInfo.Subject) && (exampleRecording.Date==recInfo.Date)
%                             
%                             
%                             if (icond==1) && (iareaHMM==1)
%                                 dataMat.example.info = exampleRecording;
%                                 dataMat.example.recInfo = recInfo;
%                                 dataMat.example.areaInfo = areaInfo;
%                             end
%                             
%                             dataMat.example.data(icond, areaIdx) = tmpData;
%                         end
%                         
                    end
                    
                    
                    %%% =============================================================================================================================================================== %%%
                case {'HMM_rate','HMM_FF'}
                    
                    switch analysisParam.HMMtype
                        case 'HMM'
                            nAreaHMM = nArea(irec);
                        case 'HMM_2area'
                            if nArea(irec)<2
                                continue
                            end
                            nAreaHMM = 1;
                    end
                    
                    stateLabel = cell(1,analysisParam.numState);
                    for istate = 1:analysisParam.numState
                        stateLabel{istate} = sprintf('State%u',istate);
                    end
                                       
                    
                    %%% HMM individually for different areas
                    for iareaHMM = 1:nAreaHMM
                        
                        switch analysisParam.HMMtype
                            case 'HMM'
                                areaLabel   = areaInfo(iareaHMM).name;
                                areaIdxHMM  = getAreaIdx(areaTable, areaLabel);

                            case 'HMM_2area'
                                areaLabel   = sprintf('%s_%s',areaInfo(1).name,areaInfo(2).name);
                                areaIdxHMM  = 1;
                        end
                        
                        read_path = getAnalysisSaveDir(recInfo, analysisParam, ievent, areaLabel, condLabel{icond} );
                        readfilename = [read_path loadfilename];
                        if ~exist(readfilename, 'file')
                            warning('file does not exist: %s', readfilename)
                            continue
                        end
                        tmpData = load(readfilename);
                        
                        for iareaNLX = 1:length(areaInfo)
                            areaIdxNLX  = getAreaIdx(areaTable, areaInfo(iareaNLX).name);

                            
                            if (rec2use==1) && (icond==1) && (iareaHMM==1) && (iareaNLX==1)
                                % Recording, Cond, areaHMM, areaNLX, state
                                dataMat.readfilename = cell(height(recordingList), numCond, max(nArea));
                                                                
                                dataMat.rate = [];
                                
                            end
                                                       
                            % find channels/units in areas
                            chanIdx = strcmpi(string(tmpData.areaList), areaInfo(iareaNLX).name);
                            unit = (tmpData.msUnitList(chanIdx,:));
                            
                            switch analysisParam.dataType
                                case 'unit'
                                otherwise
                                    unit = table(unit, 'VariableNames', {'unit'});
                            end

%                             if istable(tmpData.msUnitList)
%                                 unit = (tmpData.msUnitList.Channel);
%                             else
%                                 unit = (tmpData.msUnitList);
%                             end
%                             % get layer info
%                             chanDepth = [allDepth(unit)];
%                             chanLayer = [allLayerAssignment(unit)];
                                                        
                            numChan = length(find(chanIdx));
                            tmpTable = ([...
                                table(repmat(irec, [numChan 1]), 'VariableNames', {'Recording'}) , ...
                                table(repmat(icond, [numChan 1]), 'VariableNames', {'Condition'}) , ...
                                table(repmat(condLabel(icond), [numChan 1]), 'VariableNames', {'ConditionNames'}) , ...
                                table(repmat(iareaHMM, [numChan 1]), 'VariableNames', {'area_idx_HMM'}) , ...
                                table(repmat(string(areaLabel), [numChan 1]), 'VariableNames', {'area_HMM'}) , ...
                                table(repmat(areaIdxNLX, [numChan 1]), 'VariableNames', {'area_idx_rate'}) , ...
                                table(repmat(string(areaInfo(iareaNLX).name), [numChan 1]), 'VariableNames', {'area_rate'}) , ...
                                unit , ...
                                table(squeeze(tmpData.rate(chanIdx,:)), 'VariableNames', {'rate'}) , ...
                                table(squeeze(tmpData.sRate(chanIdx,:)), 'VariableNames', {'sRate'}) 

                                ]);

                            dataMat.rate = [...
                                dataMat.rate ; ...
                                tmpTable];
                            
                        end
                        dataMat.readfilename{irec, icond, iareaHMM} = readfilename;

                    end
                    
                    %%% =============================================================================================================================================================== %%%
                case 'HMM_RT'
                    
                    
                    for iHMM = 1:length(analysisParam.HMMtypes)
                        % test whether we recorded from multiple areas, only for specified
                        % analyses
                        if strcmpi(analysisParam.HMMtypes{iHMM},'HMM_2area') && (nArea(irec)~=2)
                            %fprintf('This recording does not have 2 areas, skipping...\n')
                            continue
                        end
                        
                        [analysisParamHMM] = specifyAnalysisParam(analysisParam.HMMtypes{iHMM}, recInfo);
                        
                        
                        switch analysisParam.HMMtypes{iHMM}
                            case {'HMM','HMMselect'}
                                %%% HMM individually for different areas
                                for iarea = 1:length(areaInfo)
                                                                                                        
                                    if (rec2use==1) && (icond==1) && (iarea==1)
                                        dataMat.readfilename = cell(height(recordingList), numCond, max(nArea));

                                        dataMat.HMM.RT = NaN(height(recordingList), numCond, max(nArea), analysisParamHMM.maxNumState);
                                        
                                        dataMat.HMM.rec = NaN(height(recordingList), numCond, max(nArea), analysisParamHMM.maxNumState);
                                        dataMat.HMM.cond = NaN(height(recordingList), numCond, max(nArea), analysisParamHMM.maxNumState);
                                        dataMat.HMM.area = NaN(height(recordingList), numCond, max(nArea), analysisParamHMM.maxNumState);
                                        dataMat.HMM.state = NaN(height(recordingList), numCond, max(nArea), analysisParamHMM.maxNumState);
                                    end
                                    
                                    areaLabel = areaInfo(iarea).name;
                                    areaIdx  = getAreaIdx(areaTable, areaLabel);

                                    read_path = getAnalysisSaveDir(recInfo, analysisParam, [1 2], areaLabel, condLabel{icond} );
                                    loadfilename = sprintf('%s_%s_HMMfit_numState_%u.mat', analysisParam.RT2use, analysisParam.RT2transform, analysisParamHMM.maxNumState);
                                    readfilename = [read_path loadfilename];
                                    if ~exist(readfilename, 'file')
                                        warning('file does not exist: %s', readfilename)
                                        continue
                                    end
                                    dataMat.readfilename{irec, icond, iarea} = readfilename;
    
                                    tmpRT = load(readfilename, 'RT_HMM');
                                                                                                            
                                    % allRT.stateArea(ipen, ievnt, icond, iarea, :) = tmpRT.RT_HMM;
                                    
                                    dataMat.HMM.RT(irec, icond, areaIdx, :) = tmpRT.RT_HMM;
                                    
                                    dataMat.HMM.rec(irec, icond, areaIdx, :) = repmat(irec, [1 analysisParamHMM.maxNumState]);
                                    dataMat.HMM.cond(irec, icond, areaIdx, :) = repmat(icond, [1 analysisParamHMM.maxNumState]);
                                    dataMat.HMM.area(irec, icond, areaIdx, :) = repmat(areaIdx, [1 analysisParamHMM.maxNumState]);
                                    dataMat.HMM.state(irec, icond, areaIdx, :) = 1:analysisParamHMM.maxNumState;
                                    
                                end
                                
                            case {'HMM_2area'}
                                
                                if (rec2use==1) && (icond==1)
                                    dataMat.readfilename = cell(height(recordingList), numCond);

                                    dataMat.HMM_2area.RT = NaN(height(recordingList), numCond, 4);
                                    dataMat.HMM_2area.rec = NaN(height(recordingList), numCond, 4);
                                    dataMat.HMM_2area.cond = NaN(height(recordingList), numCond, 4);
                                    dataMat.HMM_2area.state = NaN(height(recordingList), numCond, 4);
                                end
                                
                                areaLabel = sprintf('%s_%s',areaInfo(1).name,areaInfo(2).name);
                                    
                                read_path = getAnalysisSaveDir(recInfo, analysisParam, [1 2], areaLabel, condLabel{icond} );
                                loadfilename = sprintf('%s_%s_HMMfit_numState_%u.mat', analysisParam.RT2use, analysisParam.RT2transform, analysisParamHMM.maxNumState);
                                readfilename = [read_path loadfilename];
                                if ~exist(readfilename, 'file')
                                    warning('file does not exist: %s', readfilename)
                                    continue
                                end
                                dataMat.readfilename{irec, icond} = readfilename;

                                tmpRT = load(readfilename, 'RT_HMM');
                                
                                dataMat.HMM_2area.RT(irec, icond, :) = tmpRT.RT_HMM;
                                
                                dataMat.HMM_2area.rec(irec, icond, :) = repmat(irec, [1 analysisParamHMM.maxNumState]);
                                dataMat.HMM_2area.cond(irec, icond, :) = repmat(icond, [1 analysisParamHMM.maxNumState]);
                                dataMat.HMM_2area.state(irec, icond, :) = 1:analysisParamHMM.maxNumState;
                                
                        end
                    end
                    
                    %%% =============================================================================================================================================================== %%%
                case 'HMM_spectrogram'
                    switch analysisParam.HMMtype
                        case 'HMM'
                            nAreaHMM = nArea(irec);
                        case 'HMM_2area'
                            if nArea(irec)<2
                                continue
                            end
                            nAreaHMM = 1;
                    end
                    
                    %%% HMM individually for different areas
                    for iareaHMM = 1:nAreaHMM
                        
                        switch analysisParam.HMMtype
                            case 'HMM'
                                areaLabel = areaInfo(iareaHMM).name;
                                areaIdxHMM  = getAreaIdx(areaTable, areaLabel);
                                analysisParam.numState = 2;
                            case 'HMM_2area'
                                areaLabel = sprintf('%s_%s',areaInfo(1).name,areaInfo(2).name);
                                areaIdxHMM = 1;
                                analysisParam.numState = 4;
                        end
                        
                        read_path = getAnalysisSaveDir(recInfo, analysisParam, ievent, areaLabel, condLabel{icond} );
                        readfilename = [read_path loadfilename];
                        if ~exist(readfilename, 'file')
                            warning('file does not exist: %s', readfilename)
                            continue
                        end
                        tmpData = load(readfilename);
                                                
                        for iareaNLX = 1:length(areaInfo)
                            areaIdxNLX  = getAreaIdx(areaTable, areaInfo(iareaNLX).name);
                            
                            if (rec2use==1) && (icond==1) && (iareaHMM==1) && (iareaNLX==1)
                                % Recording, Cond, areaHMM, areaNLX, state
                                dataMat.readfilename = cell(height(recordingList), numCond, max(nArea));

                                dataMat.frequencies = tmpData.frequencies;
                                dataMat.SF = tmpData.SF;
                                
                                dataMat.spg = [];
                            end
                                                       
                            % find channels/units in areas
                            chanIdx = strcmpi(string(tmpData.areaList), areaInfo(iareaNLX).name);
                            unit = (tmpData.msUnitList(chanIdx));
                            
                            % get layer info
                            chanDepth = [allDepth(unit)];
                            chanLayer = [allLayerAssignment(unit)];
                            
                            numChan = length(find(chanIdx));
                            tmpTable = ([...
                                table(repmat(irec, [numChan 1]), 'VariableNames', {'Recording'}) , ...
                                table(repmat(icond, [numChan 1]), 'VariableNames', {'Condition'}) , ...
                                table(repmat(condLabel(icond), [numChan 1]), 'VariableNames', {'ConditionNames'}) , ...
                                table(repmat(iareaHMM, [numChan 1]), 'VariableNames', {'area_idx_HMM'}) , ...
                                table(repmat(areaLabel, [numChan 1]), 'VariableNames', {'area_HMM'}) , ...
                                table(repmat(areaIdxNLX, [numChan 1]), 'VariableNames', {'area_idx_spg'}) , ...
                                table(repmat(areaInfo(iareaNLX).name, [numChan 1]), 'VariableNames', {'area_spg'}) , ...
                                table(unit, 'VariableNames', {'channel'}) , ...
                                table(chanDepth, 'VariableNames', {'channelDepth'}) , ...
                                table(chanLayer, 'VariableNames', {'chanLayer'}) , ...
                                table(squeeze(tmpData.spg(chanIdx,:,:)), 'VariableNames', {'spg'}) 
                                
                                ]);
                            
                            dataMat.spg = [...
                                dataMat.spg ; ...
                                tmpTable];
                            
                        end
                        
                        dataMat.readfilename{irec, icond, iareaHMM} = readfilename;

                    end
                    
                    %%% =============================================================================================================================================================== %%%
                case 'HMM_transitionTime'
                    
%                     nIdx = 1:10;
%                     warning('only using select indices, only for testing!')
                    
                    %%% HMM individually for different areas
                    for iareaHMM = 1:length(areaInfo)
                        
                        areaLabel = areaInfo(iareaHMM).name;
                        areaIdxHMM  = getAreaIdx(areaTable, areaLabel);

                        read_path = getAnalysisSaveDir(recInfo, analysisParam, ievent, areaLabel, condLabel{icond} );
                        readfilename = [read_path loadfilename];
                        if ~exist(readfilename, 'file')
                            warning('file does not exist: %s', readfilename)
                            continue
                        end
                        tmpData = load(readfilename);
                    
                        if (rec2use==1) && (icond==1) && (iareaHMM==1)
                            % Recording, Cond, areaHMM, areaNLX
                            dataMat.readfilename = cell(height(recordingList), numCond, max(nArea));
                            
                            dataMat.dataTable = [];
                        end
                        
                        numState = 2;
                        tmpTable = ([...
                            table(repmat(irec, [numState 1]), 'VariableNames', {'Recording'}) , ...
                            table(repmat(icond, [numState 1]), 'VariableNames', {'Condition'}) , ...
                            table(repmat(condLabel(icond), [numState 1]), 'VariableNames', {'ConditionNames'}) , ...
                            table(repmat(iareaHMM, [numState 1]), 'VariableNames', {'area_idx'}) , ...
                            table(repmat(string(areaLabel), [numState 1]), 'VariableNames', {'area'}) , ...
                            table((1:numState)', 'VariableNames', {'state'}) , ...
                            table(squeeze(tmpData.transition_pdf), 'VariableNames', {'pdf'}), ...
                            table(squeeze(tmpData.transition_probability), 'VariableNames', {'probability'}), ...
                            table(squeeze(tmpData.transition_count), 'VariableNames', {'count'}), ...
                            table(squeeze(tmpData.transition_countDensity), 'VariableNames', {'countDensity'}), ...
                            table(squeeze(tmpData.transition_cumcount), 'VariableNames', {'cumcount'}), ...
                            table(squeeze(tmpData.transition_cdf), 'VariableNames', {'cdf'})
                            
                            ]);
                        
                        try
                        dataMat.dataTable = [...
                            dataMat.dataTable ; ...
                            tmpTable];
                        catch
                            keyboard
                        end
                        
                        
                        dataMat.readfilename{irec, icond, iareaHMM} = readfilename;
                        
                    end
                    
                    if (rec2use==1)
                        dataMat.edges = tmpData.edges;
                        dataMat.transition_label = tmpData.transition_label;
                    end
                        
                        
                    %%% =============================================================================================================================================================== %%%
                case 'HMM_transitionTriggeredAverage'
                    
                    %%% HMM individually for different areas
                    for iareaHMM = 1:length(areaInfo)
                        
                        areaLabel = areaInfo(iareaHMM).name;
                        areaIdxHMM  = getAreaIdx(areaTable, areaLabel);

                        read_path = getAnalysisSaveDir(recInfo, analysisParam, ievent, areaLabel, condLabel{icond} );
                        readfilename = [read_path loadfilename];
                        if ~exist(readfilename, 'file')
                            warning('file does not exist: %s', readfilename)
                            continue
                        end
                        tmpData = load(readfilename);

                        for iareaNLX = 1:length(areaInfo)
                            areaIdxNLX  = getAreaIdx(areaTable, areaInfo(iareaNLX).name);

                            for istate = 1:analysisParam.numState
                                
                                if (rec2use==1) && (icond==1) && (iareaHMM==1) && (iareaNLX==1) && (istate==1)
                                    % Recording, Cond, areaHMM, areaNLX, state, time
                                    dataMat.readfilename = cell(height(recordingList), numCond, max(nArea));
                                    
                                    dataMat.crossCorr = [];
                                end
                                
                                % find channels/units in areas & with state
                                % label
                                cond1 = contains((tmpData.labelCombinations), areaInfo(iareaNLX).name);
                                cond2 = contains((tmpData.labelCombinations), tmpData.label{istate});
                                chanIdx = cond1 & cond2;
                                
                                a = cellfun(@(x) regexp(x,'(?<=_).*(?=_)','match'), tmpData.labelCombinations(chanIdx),'UniformOutput',false);
                                unit = cellfun(@(x) (str2double(x{1})), a, 'UniformOutput', true);
                                
                                % get layer info
                                chanDepth = [allDepth(unit)];
                                chanLayer = [allLayerAssignment(unit)];
                                
                                numChan = length(find(chanIdx));
                                tmpTable = ([...
                                    table(repmat(irec, [numChan 1]), 'VariableNames', {'Recording'}) , ...
                                    table(repmat(icond, [numChan 1]), 'VariableNames', {'Condition'}) , ...
                                    table(repmat(condLabel(icond), [numChan 1]), 'VariableNames', {'ConditionNames'}) , ...
                                    table(repmat(iareaHMM, [numChan 1]), 'VariableNames', {'area_idx_HMM'}) , ...
                                    table(repmat(string(areaLabel), [numChan 1]), 'VariableNames', {'area_HMM'}) , ...
                                    table(repmat(areaIdxNLX, [numChan 1]), 'VariableNames', {'area_idx'}) , ...
                                    table(repmat(string(areaInfo(iareaNLX).name), [numChan 1]), 'VariableNames', {'area'}) , ...
                                    table(repmat(istate, [numChan 1]), 'VariableNames', {'state_idx'}) , ...
                                    table(repmat(string(tmpData.label{istate}), [numChan 1]), 'VariableNames', {'state_label'}) , ...
                                    table(unit, 'VariableNames', {'channel'}) , ...
                                    table(chanDepth, 'VariableNames', {'channelDepth'}) , ...
                                    table(chanLayer, 'VariableNames', {'chanLayer'}) , ...
                                    table(squeeze(tmpData.crossCorr(chanIdx,:)), 'VariableNames', {'crossCorr'})
                                    
                                    ]);
                                
                                dataMat.crossCorr = [...
                                    dataMat.crossCorr ; ...
                                    tmpTable];
                                
                            end
                        end
                        
                        dataMat.readfilename{irec, icond, iareaHMM} = readfilename;

                    end
                    
                    if (rec2use==1)
                        dataMat.timeBinsLag = tmpData.timeBinsLag;
                        dataMat.label = tmpData.label;  
                    end
                    
                case 'microsaccade_analyse_parameters'
                    read_path = getAnalysisSaveDir(recInfo, analysisParam, ievent, areaLabel, condLabel{icond} );
                    readfilename = [read_path loadfilename];
                    if ~exist(readfilename, 'file')
                        warning('file does not exist: %s', readfilename)
                        continue
                    end
                    tmpData = load(readfilename);
                    
                    if (rec2use==1) && (icond==1)
                        numStim = length(tmpData.msRateStim);
                        
                        % Recording, Cond
                        dataMat.readfilename = cell(height(recordingList), numCond);
                        
                        dataMat.msPropTrial = NaN(height(recordingList), numCond);
                        dataMat.msRate = NaN(height(recordingList), numCond);
                        dataMat.msRateQuad = NaN(height(recordingList), numCond, 4);
                        dataMat.msRateStim = NaN(height(recordingList), numCond, numStim);
                                                
                        dataMat.gaze = NaN(height(recordingList), numCond, 2);
                        dataMat.gazeRotate = NaN(height(recordingList), numCond, 2);
                        
                        %                         dataMat.gaze_gaussfit = struct([]);
                        %                         dataMat.gazeRotate_gaussfit = struct([]);
                        dataMat.dataTable = [];
                    end

                    dataMat.readfilename{irec, icond} = readfilename;
                    
                    dataMat.msPropTrial(irec, icond) = tmpData.msPropTrial;
                    dataMat.msRate(irec, icond) = tmpData.msRate;
                    dataMat.msRateQuad(irec, icond, :) = tmpData.msRateQuad;
                    dataMat.msRateStim(irec, icond, :) = tmpData.msRateStim;
                    
                    dataMat.gaze(irec, icond, :) = squeeze(nanmean(tmpData.meanGaze,1));
                    dataMat.gazeRotate(irec, icond, :) = squeeze(nanmean(tmpData.meanGazeRotate,1));
                    
                    dataMat.gaze_noMs(irec, icond, :) = squeeze(nanmean(tmpData.meanGaze(~tmpData.msTrials,:),1));
                    dataMat.gazeRotate_noMs(irec, icond, :) = squeeze(nanmean(tmpData.meanGazeRotate(~tmpData.msTrials,:),1));
                    
                    [dataMat.gaze_gaussfit(irec, icond)] = deal(tmpData.gaze_gaussfit);
                    [dataMat.gazeRotate_gaussfit(irec, icond)] = deal(tmpData.gazeRotate_gaussfit);
                    
                    [dataMat.gaze_gaussfit_noMs(irec, icond)] = deal(tmpData.gaze_gaussfit_noMs);
                    [dataMat.gazeRotate_gaussfit_noMs(irec, icond)] = deal(tmpData.gazeRotate_gaussfit_noMs);
                    
            end
        end
    end
    
    % rename target path
    tmp_path_population = regexprep(path_population,'ANALYSISNAME',analysisType);
    tmp_path_population = regexprep(tmp_path_population,'EVENTNAME',eventLabel);
    
    % store output for each event
    if ~isfolder(tmp_path_population)
        mkdir(tmp_path_population)
    end
    recordingList.include = includedRecordings;
    
    save( fullfile(tmp_path_population,savefilename), ...
        '-struct', 'dataMat');
    save( fullfile(tmp_path_population,savefilename), ...
        '-append', 'recordingList');
    save( fullfile(tmp_path_population,savefilename), ...
        '-append', 'analysisParam');
    
end

end
%% extra functions

function areaIdx = getAreaIdx(areaTable, areaName)
   
    areaIdx = areaTable.(areaName);

end



