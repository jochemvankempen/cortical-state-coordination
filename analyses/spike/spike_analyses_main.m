function spike_analyses_main(recInfo, areaInfo, analysisParam, resultSave)
% spike_analyses_main(penInfo, areaInfo, analysisParam, resultSave)
% 
% Base function for spike analyses
%
% - spike_rateAnalysis
% - spike_FFanalysis
% - spike_noiseCorrAnalysis
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
if height(recInfo)~=1
    error('height penInfo ~= 1')
end
nArea = length(areaInfo);

%%% save param
resultSave.data = true;
resultSave.overWriteExistingFiles = false;
global HPC
if HPC
    resultSave.figures = false;
else
    resultSave.figures = true;
end
resultSave.numTrialToPlot = 10;
resultSave.axisFontSize = 12;
resultSave.labelFontSize = 12;

channels2include = [];
areaLabel = [];

for iarea = 1:nArea
    channels2include = [channels2include areaInfo(iarea).ChanIdx.all(areaInfo(iarea).ChanIdx.include)];
    
    if iarea==1
        areaLabel = [areaInfo(iarea).name];
    else
        areaLabel = [areaLabel '_' areaInfo(iarea).name];
    end
end
% channels2include = boolean(channels2include);

%%% loop over events to analyse
nEvents = height(analysisParam.event);

for ievnt = 1:nEvents
    
    [condNames, nCond, ~] = getConditions(recInfo, analysisParam, analysisParam.event.ConcatenateCond(ievnt));
    
    %%% print settings on screen
    fprintf('\n------------------------------------------------------------------\n')
    fprintf('Running %s analysis with the following settings:\n', analysisParam.analysisType)
    fprintf('\tEstimate across conditions: %d\n', analysisParam.event.ConcatenateCond(ievnt))
    fprintf('\tTime epoch to align : %s, window %1.3f-%1.3f\n', analysisParam.event.Epoch{ievnt}, analysisParam.event.TimeWindow(ievnt,1), analysisParam.event.TimeWindow(ievnt,2))
    fprintf('...\n')
    
    for icond = 1:nCond
        
        fprintf('\tRunning condition: %s \n', condNames{icond})
        
        [sps, timeEpoch] = prepareData(recInfo, areaInfo, analysisParam, ...
            'Event', analysisParam.event.Epoch{ievnt}, 'Condition', icond, 'Channels',channels2include);
        
        % set paths
        target_path = getAnalysisSaveDir(recInfo, analysisParam, ievnt, areaLabel, condNames{icond} );
        fig_path = [target_path 'figures' filesep ];
        if ~exist(fig_path,'dir')
            mkdir(fig_path)
        end
        resultSave.dataDirName = target_path;
        resultSave.figuresDirName = fig_path;
        
        switch analysisParam.dataType
            case 'unit'
                sps.unitClassification = [sps.unitClassification table(sps.channel, 'VariableNames', {'Channel'})];
                sps.unitList = sps.unitClassification;
        end
        
        switch analysisParam.analysisType
            case 'rate'
                spike_rateAnalysis( sps.([analysisParam.event.Epoch{ievnt} 'Align']), timeEpoch, analysisParam, analysisParam.dataType, sps.unitList, sps.area, resultSave );
            case {'FF','FF_state'}
                spike_FFanalysis( sps.([analysisParam.event.Epoch{ievnt} 'Align']), timeEpoch, analysisParam, analysisParam.dataType, sps.unitList, sps.area, resultSave )
            case 'noiseCorr'
                spike_noiseCorrAnalysis( sps.([analysisParam.event.Epoch{ievnt} 'Align']), timeEpoch, analysisParam, analysisParam.dataType, sps.unitList, sps.area, resultSave )
        end
        
    end
end





