function HMM_RT(recInfo, areaInfo, analysisParam, resultSave)
% HMM_RT(recInfo, areaInfo, analysisParam, resultSave)
% 
% compute RT relationship to HMM states.
% Here we combine across multiple events, as the response does not always
% occur after the same event
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
% Returns
% -------
% RT_HMM : array
%     array of size (numState x 1) with the average RT for each state 
% RT_cond : array
%     array of size (numTrial x 1) with RT values
% lastState : array
%     array of size (numTrial x 1) with the state at the time of target
%     dimming
% trIdx : struct
%     struct with indices for incuded trials and condition 
%

% check input
assert(height(recInfo)==1, 'height penInfo ~= 1')
nArea = length(areaInfo);

fprintf('\tRunning HMM_RT\n')

for iHMM = 1:length(analysisParam.HMMtypes)
    
    if strcmpi(analysisParam.HMMtypes{iHMM},'HMM_2area') && (nArea~=2)
        fprintf('This recording does not have 2 areas, skipping...')
        continue
    end
    
    [analysisParamHMM] = specifyAnalysisParam(analysisParam.HMMtypes{iHMM}, recInfo);
    
    %%% loop over events to analyse
    nEvents = height(analysisParam.event);
    
    
    [condNames, nCond, ~] = getConditions(recInfo, analysisParam, analysisParam.event.ConcatenateCond(1));
        
    % get RT, concatenate across RT from multiple events
    RT = [];
    targetDim = [];
    targetDim2use = [];
    condIdx = [];
    for ievnt = 1:nEvents
        [tmpRT,~,tmpTargetDim,tmpCond] = prepareData(recInfo, areaInfo, analysisParam, ...
            'Event', analysisParam.event.Epoch{ievnt});
        RT = [RT; tmpRT]; % concatenate RT
        targetDim = [targetDim; tmpTargetDim]; % keep track of which target dimming occurred
        targetDim2use = [targetDim2use; ones(length(tmpRT),1)*ievnt];% which target dimming to use for which event
        condIdx = [condIdx; tmpCond];
    end
        
    % convert RT
    if analysisParam.concatenateAcrossEvents
        trIdx.include = (~isnan(RT)) & (targetDim==targetDim2use);% limit included trials to those that had the correct target dimming at the right event

        switch analysisParam.RT2transform
            case 'zscore_log'
                RT(trIdx.include) = zscore(log(RT(trIdx.include)));
        end
    else
        trIdx.include = (targetDim==targetDim2use);
    end

    if any(isnan(RT(trIdx.include)))
        error('NaN values in RT')
    end
    
    for icond = 1:nCond
        
        RT_cond = RT(condIdx==icond);
        trIdx.cond = trIdx.include(condIdx==icond);
        
        fprintf('\tRunning condition: %s \n', condNames{icond})
        
        switch analysisParam.HMMtypes{iHMM}
            case {'HMM'}
                %%% HMM individually for different areas
                for iarea = 1:length(areaInfo)
                    
                    fprintf('\t\tRunning area: %s \n', areaInfo(iarea).name)
                    
                    lastState = [];
                    for ievnt = 1:nEvents
                                                
                        read_path = getAnalysisSaveDir(recInfo, analysisParam, ievnt, areaInfo(iarea).name, condNames{icond}, analysisParam.HMMtypes{iHMM} );
                        HMMFilename = sprintf('HMMfit_numState_%u.mat', analysisParamHMM.maxNumState);
                        HMM = load([read_path HMMFilename], 'states','numState');
                        lastState = [lastState; cellfun(@(X) (X(end)), HMM.states)];
                        
                    end
                    
                    RT_HMM = compute_RT_HMM(RT_cond(trIdx.cond), lastState(trIdx.cond), HMM.numState, analysisParam);
                    
                    % set paths
                    target_path = getAnalysisSaveDir(recInfo, analysisParam, [1 2], areaInfo(iarea).name, condNames{icond} );
                    if ~exist(target_path,'dir')
                        mkdir(target_path)
                    end
                    
                    save([target_path sprintf('%s_%s_%s',analysisParam.RT2use,analysisParam.RT2transform,HMMFilename)], 'RT_cond','lastState','trIdx','RT_HMM')
                end
                
            case {'HMM_2area'}
                
                lastState = [];
                for ievnt = 1:nEvents
                       
                    read_path = getAnalysisSaveDir(recInfo, analysisParam, ievnt, sprintf('%s_%s',areaInfo(1).name,areaInfo(2).name), condNames{icond}, analysisParam.HMMtypes{iHMM} );
                    HMMFilename = sprintf('HMMfit_numState_%u.mat', analysisParamHMM.maxNumState);
                    HMM = load([read_path HMMFilename], 'states','numState');
                    lastState = [lastState; cellfun(@(X) (X(end)), HMM.states)];

                end
                
                RT_HMM = compute_RT_HMM(RT_cond(trIdx.cond), lastState(trIdx.cond), HMM.numState, analysisParam);

                target_path = getAnalysisSaveDir(recInfo, analysisParam, [1 2], sprintf('%s_%s',areaInfo(1).name,areaInfo(2).name), condNames{icond} );
                if ~exist(target_path,'dir')
                    mkdir(target_path)
                end
                
                save([target_path sprintf('%s_%s_%s',analysisParam.RT2use,analysisParam.RT2transform,HMMFilename)], 'RT_cond','lastState','trIdx','RT_HMM')
        end
    end
    
end


function RT_HMM = compute_RT_HMM(RT_cond, lastState, numState, analysisParam)

    if length(RT_cond) ~= length(lastState)
        error('unequal number of trials')
    end
    
    RT_HMM = NaN(numState, 1);
    for istate = 1:numState
        
        if (length(find(lastState==istate)) < analysisParam.nTrialThreshold)
            continue
        end
        RT_HMM(istate) = mean(RT_cond(lastState==istate));
    end
    
