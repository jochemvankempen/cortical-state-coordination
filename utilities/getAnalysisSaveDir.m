function target_path = getAnalysisSaveDir(recInfo, analysisParam, ievnt, areas, condition, analysisType )
% dirName = getAnalysisSaveDir(analysisParam, ievnt, areas, condition, analysisType )
%
% function retunrs the directory path where the analysis results will be
% saved, the new path contains the following sequence of labels:
%
% [local dir] / task / data / analyzed / analysisType / subject / recording / eventLabel / condLabel / epochLabel timeEpochLabel / paramLabel / area activity / includeLabel /
%
% Parameters
% ----------
% recInfo : 
% analysisParam :  struct
%     struct which fields contain analysis parameters 
% ievnt : int
%     idx which event to use
% areas : string
%     string with areas (e.g. "V1" or "V1_V4")
% condition : string 
%     string with condition name
% analysisType : string
%     string (optional) that overwrites analysisParam.analysisType
%
% Returns
% -------
% dirName : string
%     string that contains the specified path
% 

try, analysisType; catch, analysisType=analysisParam.analysisType; end

if ~isfield(analysisParam,'excludeMicrosaccades')
    analysisParam.excludeMicrosaccades = false;
end

switch analysisType
    case {'FF', 'noiseCorr'}
        switch analysisParam.windowType
            case {'single','varSize'}
                tmpSetting = sprintf('%s',analysisParam.windowType);
            case 'sliding'
                tmpSetting = sprintf('%s_%1.3f_%1.3f',analysisParam.windowType, analysisParam.windowSize, analysisParam.stepSize);
            case 'fullTimeEpoch'
                keyboard
            otherwise
                error('unknown setting for %s', analysisType)
        end
        
        dirName = [...
            sprintf('%s',analysisParam.event.Epoch{ievnt}) filesep ...eventLabel
            sprintf('%s',condition) filesep ...condLabel
            sprintf('%sEpoch_%1.3f_%1.3f',analysisParam.event.Extraction{ievnt}, analysisParam.event.TimeWindow(ievnt,1), analysisParam.event.TimeWindow(ievnt,2)) filesep ...epochLabel timeEpochLabel
            sprintf('%s',tmpSetting) filesep ...paramLabel
            sprintf('%s_%s', areas, analysisParam.dataType) filesep ...area activity
            sprintf('exMs%d', analysisParam.excludeMicrosaccades) filesep ... includeLabel
            ];
        
        
    case 'rate'     
        switch analysisParam.windowType
            case 'window'
                tmpSetting = sprintf('%s_%1.3f_%1.3f',analysisParam.windowType, analysisParam.windowSize, analysisParam.stepSize);
            case 'fullTimeEpoch'
                tmpSetting = sprintf('%s',analysisParam.windowType);
            otherwise
                error('unknown setting for %s', analysisType)
        end
        
        dirName = [...
            sprintf('%s',analysisParam.event.Epoch{ievnt}) filesep ...eventLabel
            sprintf('%s',condition) filesep ...condLabel
            sprintf('%sEpoch_%1.3f_%1.3f',analysisParam.event.Extraction{ievnt}, analysisParam.event.TimeWindow(ievnt,1), analysisParam.event.TimeWindow(ievnt,2)) filesep ...epochLabel timeEpochLabel
            sprintf('%s_matchRateCond%d',tmpSetting, analysisParam.matchRateAcrossCond) filesep ...paramLabel
            sprintf('%s_%s', areas, analysisParam.dataType) filesep ...area activity
            sprintf('exMs%d', analysisParam.excludeMicrosaccades) filesep ... includeLabel
            ];

        
        
    case 'HMM_RT'
        
        if length(ievnt)~=2
            error('directory currently coded for 2 events')
        end
        
        dirName = [...
            sprintf('%s_%s',analysisParam.event.Epoch{ievnt(1)},analysisParam.event.Epoch{ievnt(2)}) filesep ...% eventLabel
            sprintf('%s',condition) filesep ...% condLabel
            sprintf('numFold%u_binSize%1.3f',analysisParam.numFold, analysisParam.binSize) filesep ...paramLabel
            sprintf('%s_%s', areas, analysisParam.dataType) filesep ...area activity
            sprintf('exMs%d', analysisParam.excludeMicrosaccades) filesep ... includeLabel
            ];
        
    case {'HMM','HMMselect','HMM_2area','HMM_cvStateBased', 'HMM_transitionTime'}
        % cross validation fold. 
        
        if analysisParam.matchRateAcrossCond
            dirName = [...
                sprintf('%s',analysisParam.event.Epoch{ievnt}) filesep ... % eventLabel
                sprintf('%s',condition) filesep ... % condLabel
                sprintf('%sEpoch_%1.3f_%1.3f',analysisParam.event.Extraction{ievnt}, analysisParam.event.TimeWindow(ievnt,1), analysisParam.event.TimeWindow(ievnt,2)) filesep ... epochLabel timeEpochLabel
                sprintf('numFold%u_binSize%1.3f_matchRateCond%d',analysisParam.numFold, analysisParam.binSize, analysisParam.matchRateAcrossCond) filesep ... paramLabel
                sprintf('%s_%s', areas, analysisParam.dataType) filesep ... area activity
                sprintf('exMs%d', analysisParam.excludeMicrosaccades) filesep ... includeLabel
                ];
        else
            dirName = [...
                sprintf('%s',analysisParam.event.Epoch{ievnt}) filesep ... % eventLabel
                sprintf('%s',condition) filesep ... % condLabel
                sprintf('%sEpoch_%1.3f_%1.3f',analysisParam.event.Extraction{ievnt}, analysisParam.event.TimeWindow(ievnt,1), analysisParam.event.TimeWindow(ievnt,2)) filesep ... epochLabel timeEpochLabel
                sprintf('numFold%u_binSize%1.3f',analysisParam.numFold, analysisParam.binSize) filesep ... paramLabel
                sprintf('%s_%s', areas, analysisParam.dataType) filesep ... area activity
                sprintf('exMs%d', analysisParam.excludeMicrosaccades) filesep ... includeLabel
                ];
        end
%         
%         dirName = [...
%             sprintf('%s',analysisParam.event.Epoch{ievnt}) filesep ... % eventLabel
%             sprintf('%s',condition) filesep ... % condLabel
%             sprintf('%sEpoch_%1.2f_%1.2f_',analysisParam.event.Extraction{ievnt}, analysisParam.event.TimeWindow(ievnt,1), analysisParam.event.TimeWindow(ievnt,2))  ... epochLabel timeEpochLabel
%             sprintf('numFold%u_binSize%1.2f',analysisParam.numFold, analysisParam.binSize) filesep ... paramLabel
%             sprintf('%s_%s', areas, analysisParam.dataType) filesep ... area activity
%             ];
        
    case {'HMM_crossCorrelation'}
        dirName = [...
            sprintf('%s',analysisParam.event.Epoch{ievnt}) filesep ... % eventLabel
            sprintf('%s',condition) filesep ... % condLabel
            sprintf('%sEpoch_%1.3f_%1.3f',analysisParam.event.Extraction{ievnt}, analysisParam.event.TimeWindow(ievnt,1), analysisParam.event.TimeWindow(ievnt,2)) filesep ... epochLabel timeEpochLabel
            sprintf('numFold%u_binSize%1.3f',analysisParam.numFold, analysisParam.binSize) filesep ... paramLabel
            sprintf('%s_%s', areas, analysisParam.dataType) filesep ... area activity
            sprintf('exMs%d_minDur%1.3f', analysisParam.excludeMicrosaccades, analysisParam.minEpochDurationWin) filesep ... includeLabel
            ];
        
    case {'HMM_microsaccades'}
        dirName = [...
            sprintf('%s',analysisParam.event.Epoch{ievnt}) filesep ... % eventLabel
            sprintf('%s',condition) filesep ... % condLabel
            sprintf('%sEpoch_%1.3f_%1.3f',analysisParam.event.Extraction{ievnt}, analysisParam.event.TimeWindow(ievnt,1), analysisParam.event.TimeWindow(ievnt,2)) filesep ... epochLabel timeEpochLabel
            sprintf('numFold%u_binSize%1.3f_timeLag%1.3f', analysisParam.numFold, analysisParam.binSize, analysisParam.maxTimeLag) filesep ... paramLabel
            sprintf('%s_%s', areas, analysisParam.dataType) filesep ... area activity
            sprintf('exMs%d', analysisParam.excludeMicrosaccades) filesep ... includeLabel
            ];

    case 'HMM_rate'
        dirName = [...
            sprintf('%s',analysisParam.event.Epoch{ievnt}) filesep ...eventLabel
            sprintf('%s',condition) filesep ...condLabel
            sprintf('%sEpoch_%1.3f_%1.3f',analysisParam.event.Extraction{ievnt}, analysisParam.event.TimeWindow(ievnt,1), analysisParam.event.TimeWindow(ievnt,2)) filesep ...epochLabel timeEpochLabel
            sprintf('numFold%u_binSize%1.3f_matchRateCond%d',analysisParam.numFold, analysisParam.binSize, analysisParam.matchRateAcrossCond) filesep ... paramLabel
            sprintf('%s_minDuration%1.3f', analysisParam.HMMtype, analysisParam.minEpochDurationWin) filesep ... paramLabel2
            sprintf('%s_%s', areas, analysisParam.dataType) filesep ...area activity
            sprintf('exMs%d', analysisParam.excludeMicrosaccades) filesep ... includeLabel
            ];
        
    case {'HMM_spectrogram'}
        dirName = [...
            sprintf('%s',analysisParam.event.Epoch{ievnt}) filesep ...eventLabel
            sprintf('%s',condition) filesep ...condLabel
            sprintf('%sEpoch_%1.3f_%1.3f',analysisParam.event.Extraction{ievnt}, analysisParam.event.TimeWindow(ievnt,1), analysisParam.event.TimeWindow(ievnt,2)) filesep ...epochLabel timeEpochLabel
            sprintf('numFold%u_binSize%1.3f',analysisParam.numFold, analysisParam.binSize) filesep ...paramLabel
            sprintf('%s_minDuration%1.3f_tapers%d_%d_fpass%d_%d', analysisParam.HMMtype, analysisParam.minEpochDurationWin, analysisParam.tapers(1), analysisParam.tapers(2), analysisParam.fpass(1), analysisParam.fpass(2)) filesep ...paramLabel2
            sprintf('%s_%s', areas, analysisParam.dataType) filesep ...area activity
            sprintf('exMs%d', analysisParam.excludeMicrosaccades) filesep ... includeLabel            
            ];
    
    case 'HMM_transitionTriggeredAverage'
        dirName = [...
            sprintf('%s',analysisParam.event.Epoch{ievnt}) filesep ...% eventLabel 
            sprintf('%s',condition) filesep ...condLabel
            sprintf('%sEpoch_%1.3f_%1.3f',analysisParam.event.Extraction{ievnt}, analysisParam.event.TimeWindow(ievnt,1), analysisParam.event.TimeWindow(ievnt,2)) filesep ...epochLabel timeEpochLabel
            sprintf('numFold%u_binSize%1.3f_timeLag%1.3f',analysisParam.numFold, analysisParam.binSize, analysisParam.maxTimeLag) filesep ... paramLabel
            sprintf('%s_minDuration%1.3f_demean%d', analysisParam.HMMtype, analysisParam.minEpochDurationWin, analysisParam.demean) filesep ... paramLabel2
            sprintf('%s_%s', areas, analysisParam.dataType) filesep ...area activity
            sprintf('exMs%d', analysisParam.excludeMicrosaccades) filesep ... includeLabel
            ];
        
    case 'microsaccade_analyse_parameters'
        dirName = [...
            sprintf('%s',analysisParam.event.Epoch{ievnt}) filesep ...% eventLabel
            sprintf('%s',condition) filesep ...condLabel
            sprintf('%sEpoch_%1.3f_%1.3f',analysisParam.event.Extraction{ievnt}, analysisParam.event.TimeWindow(ievnt,1), analysisParam.event.TimeWindow(ievnt,2)) filesep ...epochLabel timeEpochLabel
            ];
        
        
    otherwise
        
        dirName = [...
            sprintf('%s',analysisParam.event.Epoch{ievnt}) filesep ...eventLabel
            sprintf('%s',condition) filesep ...condLabel
            sprintf('%sEpoch_%1.3f_%1.3f',analysisParam.event.Extraction{ievnt}, analysisParam.event.TimeWindow(ievnt,1), analysisParam.event.TimeWindow(ievnt,2)) filesep ...epochLabel timeEpochLabel
            sprintf('binSize%1.3f',analysisParam.binSize) filesep ...paramLabel
            sprintf('%s_%s', areas, analysisParam.dataType) filesep ...area activity
            sprintf('exMs%d', analysisParam.excludeMicrosaccades) filesep ... includeLabel     
            ];
        
        
end

target_path = [recInfo.path_target dirName];
target_path = regexprep(target_path, 'ANALYSISNAME', analysisType);
if ~exist(target_path,'dir')
    mkdir(target_path)
end
