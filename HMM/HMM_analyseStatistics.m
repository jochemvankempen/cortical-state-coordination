function HMM_analyseStatistics(recInfo, areaInfo, analysisParam, resultSave)
% HMM_analyseStatistics(penInfo, areaInfo, analysisParam, resultSave)
% 
% base function for testing statistics of HMM phases/states
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
        
        %%% HMM individually for different areas
        for iarea = 1:length(areaInfo)
            
            fprintf('\t\tRunning area: %s \n', areaInfo(iarea).name)
            areaString = areaInfo(iarea).name;

            [HMM, timeEpoch] = prepareData(recInfo, areaInfo, analysisParam, ...
                'Event', analysisParam.event.Epoch{ievnt}, 'Condition', icond, 'Area', areaInfo(iarea).name);
            
            if isempty(HMM)
                continue
            end
            
            states = HMM.states;
            timeBins = HMM.timeBins;
                        
            % set paths
            target_path = getAnalysisSaveDir(recInfo, analysisParam, ievnt, areaString, condNames{icond} );
            fig_path = [target_path 'figures' filesep ];
            if ~exist(fig_path,'dir')
                mkdir(fig_path)
            end
            
            resultSave.dataDirName = target_path;
            
            switch analysisParam.analysisType
                case 'HMM_transitionTime'
                    HMM_analyseTransitionTimes(states, timeBins, analysisParam, areaInfo(iarea).name, resultSave);
            end

        end
    end
end





