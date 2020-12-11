function HMM_crossCorrelation(recInfo, areaInfo, analysisParam, resultSave)
% HMM_crossCorrelation(recInfo, areaInfo, analysisParam, resultSave)
% 
% cross correlation between HMM applied to different areas
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
%

% check input
assert(height(recInfo)==1, 'height penInfo ~= 1')
nArea = length(areaInfo);

if (nArea~=2)
    fprintf('This recording does not have 2 areas, skipping...\n')
    return
end

switch analysisParam.analysisType
    case 'HMM_crossCorrelationPupil'
        if isempty(recInfo.EyeFileName{1})
            warning('This recording does not have pupil data, skipping...\n')
            return
        end
end

%%% loop over events to analyse
nEvents = height(analysisParam.event);

for ievnt = 1:nEvents
    
    [condNames, nCond, ~] = getConditions(recInfo, analysisParam, analysisParam.event.ConcatenateCond(ievnt));
    
    %%% print settings on screen
    fprintf('\n------------------------------------------------------------------\n')
    fprintf('Running %s analysis with the following settings:\n', analysisParam.analysisType)
    fprintf('\tEstimate across conditions: %d\n', analysisParam.event.ConcatenateCond(ievnt))
    fprintf('\tTime epoch to align : %s, window %1.3f-%1.3f\n', analysisParam.event.Epoch{ievnt}, analysisParam.event.TimeWindow(ievnt,1), analysisParam.event.TimeWindow(ievnt,2))
%     fprintf('States: %d-%d \n', analysisParam.numState, analysisParam.numState)
    fprintf('\t...\n')
    
    for icond = 1:nCond
        
        fprintf('\tRunning condition: %s \n', condNames{icond})
        
        %%% HMM individually for different areas
        for iarea = 1:length(areaInfo)
                        
            switch analysisParam.analysisType
                case 'HMM_crossCorrelationPupil'
                    [HMM, timeEpoch, argOut] = prepareData(recInfo, areaInfo, analysisParam, ...
                        'Event', analysisParam.event.Epoch{ievnt}, 'Condition', icond, 'Area', areaInfo(iarea).name);
                    
                    idx = find(strcmpi(argOut(1:2:end),'pupil'));
                    assert( ~isempty(idx) && length(idx)==1, 'prepareData did not load pupil correctly')
                    pupil = argOut{idx+1};
                    
                otherwise
                    [HMM, timeEpoch] = prepareData(recInfo, areaInfo, analysisParam, ...
                        'Event', analysisParam.event.Epoch{ievnt}, 'Condition', icond, 'Area', areaInfo(iarea).name);
            end

            if isempty(HMM)
                continue
            end
            
            if iarea==1
                states = cell(length(HMM.states),2);
                timeBins = cell(length(HMM.timeBins),2);
            end
            
            states(:,iarea) = HMM.states;
            timeBins(:,iarea) = HMM.timeBins;
            
        end
        
        if isempty(HMM)
            continue
        end
        
        % extract info about state switches
        stateswitch     = cellfun(@(X)(find([0 diff(X)])), states, 'UniformOutput', false); % when did a stateswitch happen
        stateswitch_idx = ~cellfun(@isempty, stateswitch, 'UniformOutput', true); % on which trials was there a transition
        
        trIdx           = (sum(stateswitch_idx, 2)==2); % get trials where at least one switch ocurred in both areas
        
        if isempty(find(trIdx)) || length(find(trIdx))==1
            continue
        end
        
        resultSave.dataDirName = getAnalysisSaveDir(recInfo, analysisParam, ievnt, sprintf('%s_%s',areaInfo(1).name,areaInfo(2).name), condNames{icond} );
        
        switch analysisParam.CC_type
            case 'timeSeries'
                
                switch analysisParam.analysisType
                    case 'HMM_crossCorrelation'
                        % perform cross correlation on timeseries
                        CC_timeSeries(states(trIdx, 2:end)', states(trIdx, 1)', 1/analysisParam.binSize, analysisParam.maxTimeLag, analysisParam, resultSave);

                    case 'HMM_crossCorrelationPupil'
                        
                        [~, binIdx, ~] = binVal(pupil, analysisParam.nPupilBin);

                        if ( length(binIdx)~=size(states,1) )
                            error('Unequal trial numbers')
                        end
                        
                        for ibin = 1:analysisParam.nPupilBin
                            
                            if isempty(find(trIdx & (binIdx==ibin))) || length(find(trIdx & (binIdx==ibin)))==1
                                continue
                            end
                            % perform cross correlation on timeseries, for bins
                            % of trials                            
                            CC_timeSeries(states(trIdx & (binIdx==ibin), 2:end)', states(trIdx & (binIdx==ibin), 1)', 1/analysisParam.binSize, analysisParam.maxTimeLag, analysisParam, resultSave);

                        end
                        filename = fullfile(resultSave.dataDirName, sprintf('HMM_crossCorrelation_numstate%d_exms%d_nbin%u.mat', analysisParam.numState, analysisParam.excludeMicrosaccades, analysisParam.nPupilBin));
                        save(filename, 'ccMat');
                end
            otherwise
                % perform cross correlation on transition times.
                % change around order of states, last channel is trigger, second to
                % last is timeBins
                clear tmpstates;
                dataType = 'HMM';
                timeEpoch = timeEpoch(1,:);
                states = states';
                tmpstates(1,:) = states(2,:); % take V4 as signal
                tmpstates(2,:) = timeBins(:,iarea)';
                tmpstates(3,:) = states(1,:); % take V1 as trigger
                msUnitList = {'V4'};
                areaList = {'V4'};
                                
                CC_HMM_TA(tmpstates, timeEpoch, [], analysisParam, dataType, msUnitList, areaList, resultSave)
        end
        
    end
end





