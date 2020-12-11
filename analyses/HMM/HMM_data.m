function HMM_data(recInfo, areaInfo, analysisParam, resultSave)
% HMM_data(penInfo, areaInfo, analysisParam, resultSave)
% 
% base function for testing relationship between HMM phases/states and data, e.g.
% transition triggered average
% 
% Parameters
% ----------
% recInfo : table
%     table with single row (one recording)
% areaInfo : struct
%     struct with info about each area, obtained by 'getAreaInfo'
% analysisParam : struct
%     struct with analysis parameters, obtained by 'specifyAnalysisParam
% resultSave : struct
%     struct with fields about saving data/figures
%


% check input
assert(height(recInfo)==1, 'height recInfo ~= 1')

if contains(analysisParam.dataType,'_conv')
    error('Come on mate, spike convolution is controlled by analysisParam.convolveSpikes, You know better than that')
end

nArea = length(areaInfo);

%%% check whether analysis is appropriate for this recording
switch analysisParam.HMMtype
    case 'HMM_2area'
        if nArea<=1
            warning('This recording does not have multiple areas, skipping...\n')
            return
        end
end

switch analysisParam.analysisType
    case 'HMM_pupil'
        if ~exist( fullfile(recInfo.path_read,'pupil.mat'), 'file' )
            warning('This recording does not have pupil data, skipping...\n')
            return
        end
end

channels2include = [];
switch analysisParam.dataType
    case {'LFPb'}
        for iarea = 1:nArea
            channels2include = [channels2include areaInfo(iarea).ChanIdx.all(areaInfo(iarea).ChanIdx.include_bip)];
        end
    otherwise
        for iarea = 1:nArea
            channels2include = [channels2include areaInfo(iarea).ChanIdx.all(areaInfo(iarea).ChanIdx.include)];
        end
end

%%% loop over events to analyse
nEvents = height(analysisParam.event);
for ievnt = 1:nEvents
    
    [condNames, nCond, ~] = getConditions(recInfo, analysisParam, analysisParam.event.ConcatenateCond(ievnt));
    
    %%% print settings on screen
    fprintf('\t------------------------------------------------------------------\n')
    fprintf('\tRunning %s analysis with the following settings:\n', analysisParam.analysisType)
    fprintf('\tEstimate across conditions: %d\n', analysisParam.event.ConcatenateCond(ievnt))
    fprintf('\tTime epoch to align : %s, window %1.3f-%1.3f\n', analysisParam.event.Epoch{ievnt}, analysisParam.event.TimeWindow(ievnt,1), analysisParam.event.TimeWindow(ievnt,2))
%     fprintf('States: %d-%d \n', analysisParam.numState, analysisParam.numState)
    
    for icond = 1:nCond
        
        fprintf('\tRunning condition: %s \n', condNames{icond})
        
        [readData, timeEpoch, HMM, positionStimuli] = prepareData(recInfo, areaInfo, analysisParam, ...
            'Event', analysisParam.event.Epoch{ievnt}, 'Condition', icond, 'Channels', channels2include);
        
        if isempty(readData)
            break
        end
        
        idx = find(strcmpi(HMM(1:2:end),'HMM'));
        assert( ~isempty(idx) && length(idx)==1, 'prepareData did not load HMM correctly')
        HMM = HMM{idx+1};

        nArea = size(HMM.states,2); % limit looping over multiple areas when using HMM_2area
        
        %%% HMM individually for different areas
        for iareaHMM = 1:nArea
            
            fprintf('\t\tRunning area: %s \n', areaInfo(iareaHMM).name)
            
            switch analysisParam.HMMtype
                case 'HMM'
                    areaString = areaInfo(iareaHMM).name;
                case 'HMM_2area'
                    areaString = sprintf('%s_%s',areaInfo(1).name,areaInfo(2).name);
            end
            
            %%% add HMM timebins and states as extra channels
            tmpAnalysisParam = analysisParam;
            switch analysisParam.dataType
                case {'hash','unit','MUA_spont20'}
                    msu = [readData.( [analysisParam.event.Epoch{ievnt} 'Align'] ); HMM.timeBins'; HMM.states(:,iareaHMM)' ];
                    SF = 1000;
                    
                    if analysisParam.convoluteSpikes
                        tmpAnalysisParam.dataType = [analysisParam.dataType '_conv'];
                    end
                    
                    switch analysisParam.dataType
                        case 'unit'
                            readData.unitList = [readData.unitClassification table(readData.channel,'VariableNames',{'Channel'})];
                    end
                    
                case {'NCS','LFP','LFPb','MUAe'}
                    
                    switch analysisParam.analysisType
                        case 'HMM_transitionTriggeredAverage'
                            msu = [readData.( [analysisParam.event.Epoch{ievnt} 'Align'] ).Samples; HMM.timeBins'; HMM.states(:,iareaHMM)' ];
                        otherwise
                            analogSig = readData.( [analysisParam.event.Epoch{ievnt} 'Align'] ) ;
                    end
                    SF = readData.SF;
                    
                case {'microsaccades'}
                    
                    areaIdx = (1:length(HMM.msulabel)) + (length(HMM.msulabel)*(iareaHMM-1));
                    
                    % check whether there are any microsaccades at all
                    idx_field = ~arrayfun(@(X)( isfield(X, 'StartTime')), readData.( [analysisParam.event.Epoch{ievnt} 'Align'] ).ms, 'UniformOutput', true);
                    idx_empty = arrayfun(@(X)( isempty(X.StartTime)), readData.( [analysisParam.event.Epoch{ievnt} 'Align'] ).ms, 'UniformOutput', true);
                    
                    if all(idx_field) || all(idx_empty)
                        continue
                    end
                        
                    ms_startTime = arrayfun(@(X)(X.StartTime), readData.( [analysisParam.event.Epoch{ievnt} 'Align'] ).ms, 'UniformOutput', false);
                    ms_theta = arrayfun(@(X)(X.Phi), readData.( [analysisParam.event.Epoch{ievnt} 'Align'] ).ms, 'UniformOutput', false);
                    
                    msu = [HMM.msu(:,areaIdx)' ; ms_startTime ];
                    msuLabel = [HMM.msulabel {'ms'}];
                                        
                    if iareaHMM==1
                        % get the stimulus position (theta)
                        idx = find(strcmpi(positionStimuli(1:2:end),'positionStimuli'));
                        assert( ~isempty(idx) && length(idx)==1, 'prepareData did not load positionStimuli correctly')
                        positionStimuli = positionStimuli{idx+1};
                        
                        % attendedStimulus = ~cellfun(@isempty, regexp( condNames{icond}, positionStimuli.Properties.VariableNames) );
                        attendedStimulus = ~cellfun(@isempty, regexp( 'RF', positionStimuli.Properties.VariableNames) );
                        assert(sum(attendedStimulus)==1, 'found multiple conditions')
                        
                        stimfields = positionStimuli.Properties.VariableNames;
                        
                        thetaStimulus = atan2(positionStimuli.(stimfields{attendedStimulus})(2), positionStimuli.(stimfields{attendedStimulus})(1));
                    end
                case 'pupil'
                    msu{1} = readData(:)';
                    msu{2} = [HMM.timeBins'; HMM.states(:,iareaHMM)' ];
                    
                case 'HMM'
                    msu = [HMM.timeBins'; HMM.states(:,iareaHMM)' ];
            end     
            
            % set paths
            target_path = getAnalysisSaveDir(recInfo, tmpAnalysisParam, ievnt, areaString, condNames{icond} );
            fig_path = [target_path 'figures' filesep ];
            if ~exist(fig_path,'dir')
                mkdir(fig_path)
            end
            
            resultSave.dataDirName = target_path;
            
            switch analysisParam.analysisType
                case 'HMM_microsaccades'
%                     switch analysisParam.CC_type
%                         case 'cpt'
                            CC_HMM_MS_cpt(msu, msuLabel, timeEpoch, analysisParam, resultSave);
%                         case 'jitter'
                            % CC_HMM_MS_jitter(msu, msuLabel, timeEpoch, analysisParam, resultSave);
%                     end

                    msu = [msu ; ms_theta ];
                    msuLabel = [msuLabel {'theta'}];

                    microsaccade_HMM_rate(msu, msuLabel, thetaStimulus, timeEpoch, analysisParam, resultSave);

                case 'HMM_pupil'
                    HMM_correlateEpochDuration(msu, analysisParam, tmpAnalysisParam.dataType, resultSave);
                case 'HMM_rate'
                    HMM_rate(msu, analysisParam, tmpAnalysisParam.dataType, readData.unitList, readData.area, resultSave);
                case 'HMM_transitionTriggeredAverage'
                    %timeEpoch = timeEpoch(1,:);
                    CC_HMM_TA(msu, timeEpoch, SF, analysisParam, tmpAnalysisParam.dataType, readData.unitList, readData.area, resultSave);
                case 'HMM_spectrogram'
                    HMM_spectrogram(analogSig, HMM.states(:,iareaHMM), HMM.timeBins, SF, analysisParam, tmpAnalysisParam.dataType, readData.unitList, readData.area, resultSave);
            end

        end
    end
end





