function HMM_main(recInfo, areaInfo, analysisParam, resultSave)
% HMM_main(recInfo, areaInfo, analysisParam, resultSave)
% 
% Base function for HMM model fitting.
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
nArea = length(areaInfo);

switch analysisParam.analysisType
    case {'HMM_2area'}
        if (nArea~=2)
            fprintf('This recording does not have 2 areas, skipping...\n')
            return
        end
end

layerAssignment = [areaInfo.layerAssignment];
layerAssignment = reshape(layerAssignment, numel(layerAssignment), 1);

if length(areaInfo)==1 && strcmpi(areaInfo(1).name, 'V4')
    layerAssignment = vertcat( cellstr(repmat('none',[length(areaInfo(1).ChanIdx.all) 1])), layerAssignment);
end

% start parfor env
global HPC
if HPC
    SLURM_parpool_init
end

%%% save param
resultSave.data = true;
resultSave.overWriteExistingFiles = true;
if HPC
    resultSave.figures = false;
else
    resultSave.figures = true;
end
    resultSave.figures = false;

resultSave.numTrialToPlot = 10;
resultSave.axisFontSize = 12;
resultSave.labelFontSize = 12;


%%% loop over events to analyse
nEvents = height(analysisParam.event);

for ievnt = 1:nEvents
    
    [condNames, nCond, ~] = getConditions(recInfo, analysisParam, analysisParam.event.ConcatenateCond(ievnt));
    
    %%% print settings on screen
    fprintf('\n------------------------------------------------------------------\n')
    fprintf('Running analysis: %s with the following settings:\n', analysisParam.analysisType)
    fprintf('\tEstimate across conditions: %d\n', analysisParam.event.ConcatenateCond(ievnt))
    fprintf('\tTime epoch to align : %s, %s, window %1.3f-%1.3f\n', analysisParam.event.Epoch{ievnt}, analysisParam.event.Extraction{ievnt}, analysisParam.event.TimeWindow(ievnt,1), analysisParam.event.TimeWindow(ievnt,2))
    fprintf('\tStates: %d-%d \n', analysisParam.minNumState, analysisParam.maxNumState)
    fprintf('...\n')
    
    for icond = 1:nCond
        
        fprintf('\tRunning condition: %s \n', condNames{icond})
        
        switch analysisParam.analysisType
            case {'HMM','HMMselect','HMM_cvStateBased'}
                
                %%% HMM individually for different areas
                for iarea = 1:length(areaInfo)
                    
                    fprintf('\t\tRunning area: %s \n', areaInfo(iarea).name)

                    % get data
                    [sps, timeEpoch] = prepareData(recInfo, areaInfo, analysisParam, ...
                        'Event', analysisParam.event.Epoch{ievnt}, 'Condition', icond, 'Channels', areaInfo(iarea).ChanIdx.all(areaInfo(iarea).ChanIdx.include));
                    
%                     continue
                    % set paths
                    target_path = getAnalysisSaveDir(recInfo, analysisParam, ievnt, areaInfo(iarea).name, condNames{icond} );
                    fig_path = [target_path 'figures' filesep ];
                    if ~exist(fig_path,'dir')
                        mkdir(fig_path)
                    end
                    resultSave.dataDirName = target_path;
                    resultSave.figuresDirName = fig_path;
                    
                    recording.msUnitList    = num2cell(sps.unitList);
                    
                    switch analysisParam.analysisType
                        case 'HMM_cvStateBased'
                            % in case of HMM_cvStateBased, add dataType
                            recording.dataType = cellstr(sps.dataType);
                            
                            recording.layerData = layerAssignment(sps.channel);
                        otherwise
                            recording.layerData = layerAssignment(sps.unitList);
                    end
                    
                    % run analysis
                    switch analysisParam.analysisType
                        case 'HMM'
                            hmmAnalysisCrossVal( sps.([analysisParam.event.Epoch{ievnt} 'Align']), ...
                                analysisParam.binSize, timeEpoch, analysisParam.minNumState, analysisParam.maxNumState, analysisParam.numFold, recording, resultSave )
                        case 'HMMselect'
                            hmmModelSelection( sps.([analysisParam.event.Epoch{ievnt} 'Align']), ...
                                analysisParam.binSize, timeEpoch, analysisParam.maxNumState, analysisParam.numFold, recording, resultSave )
                        case 'HMM_cvStateBased'                            
                            crossValBasedOnState( sps.([analysisParam.event.Epoch{ievnt} 'Align']), ...
                                analysisParam.binSize, timeEpoch, analysisParam.numFold, recording, resultSave );
                    end
                end
                
            case {'HMM_2area'}
                % HMM across areas
                
                % get data
                [sps, timeEpoch] = prepareData(recInfo, areaInfo, analysisParam, ...
                    'Event', analysisParam.event.Epoch{ievnt}, 'Condition', icond, 'Channels', [areaInfo(1).ChanIdx.include areaInfo(2).ChanIdx.include]);
                         
                % set paths
                target_path = getAnalysisSaveDir(recInfo, analysisParam, ievnt, sprintf('%s_%s',areaInfo(1).name,areaInfo(2).name), condNames{icond} );
                fig_path = [target_path 'figures' filesep ];
                if ~exist(fig_path,'dir')
                    mkdir(fig_path)
                end
                resultSave.dataDirName = target_path;
                resultSave.figuresDirName = fig_path;
                
                recording.msUnitList    = num2cell(sps.unitList);
                recording.layerData     = [areaInfo(1).layerAssignment(areaInfo(1).ChanIdx.include) ; areaInfo(2).layerAssignment(areaInfo(2).ChanIdx.include)];
                
                % define from which area a channel was recorded
                chanPerArea = sum([strcmpi(string(sps.area),areaInfo(1).name), strcmpi(string(sps.area),areaInfo(2).name) * 2], 2);
                
                % run analysis
                hmmAnalysisCrossVal_constrE( sps.([analysisParam.event.Epoch{ievnt} 'Align']), ...
                    chanPerArea, analysisParam.binSize, timeEpoch, analysisParam.minNumState, analysisParam.maxNumState, analysisParam.numFold, recording, resultSave )

        end
    end
end





