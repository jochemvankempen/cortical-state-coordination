function perform_analyses(recInfo,Job)
% perform_analyses(recInfo,Job)
% 
% Wrapper function from which analyses are called
% 
% Parameters
% ----------
% recInfo : table
%     table with single row (one recording)
% Job : struct
%     structure with fields refering to individual analyses. Each field
%     is a boolean indicating whether the analysis will be executed
% 

tStart = tic;

%% check input & set paths

assert(height(recInfo)==1, 'height recInfo ~= 1')

%% print settings on screen
fprintf('\tRead path: %s\n', recInfo.path_read)
fprintf('\tTarget path: %s\n', recInfo.path_target)

%% get info about area: channel exclusion, RF etc. 

[areaInfo, areaInfoWarning] = getAreaInfo(recInfo, recInfo.path_read);% areaInfo
[recInfo, areaInfo]         = getRFinfo(recInfo, areaInfo, recInfo.path_read);

if areaInfoWarning
    warning('Could not load areaInfo')
%     return
end

%% parameters for the saved data and figures
resultSave.data = 1;    % whether to save data
resultSave.figures = 0; % whether to make figures
resultSave.figureEachChannel = 0;   % whether to make figures for each channel
resultSave.axisFontSize = 22;
resultSave.labelFontSize = 24;

%% Run Jobs

%%% ------------------------------------------------------------------- %%%
%%% Spike
%%% ------------------------------------------------------------------- %%%

% rate
if Job.analyse_spike_rate
    [analysisParam] = specifyAnalysisParam('rate', recInfo);
    spike_analyses_main(recInfo, areaInfo, analysisParam, resultSave)
end

% Fano Factor
if Job.analyse_spike_FF
    [analysisParam] = specifyAnalysisParam('FF', recInfo);
    spike_analyses_main(recInfo, areaInfo, analysisParam, resultSave)
end

%%% ------------------------------------------------------------------- %%%
%%% HMM
%%% ------------------------------------------------------------------- %%%

%%% ------------------------------------------------------------------- %%%
%%% HMM fit
if Job.HMM_fit
    [analysisParam] = specifyAnalysisParam('HMM', recInfo);
    HMM_main(recInfo, areaInfo, analysisParam, resultSave)
end

if Job.HMM_cvStateBased
    [analysisParam] = specifyAnalysisParam('HMM_cvStateBased', recInfo);
    HMM_main(recInfo, areaInfo, analysisParam, resultSave)
end

if Job.HMM_modelSelection
    [analysisParam] = specifyAnalysisParam('HMMselect', recInfo);
    HMM_main(recInfo, areaInfo, analysisParam, resultSave)
end

if Job.HMM_2area
    [analysisParam] = specifyAnalysisParam('HMM_2area', recInfo);
    HMM_main(recInfo, areaInfo, analysisParam, resultSave)
end



%%% ------------------------------------------------------------------- %%%
%%% HMM - HMM
if Job.HMM_crossCorrelation
    [analysisParam] = specifyAnalysisParam('HMM_crossCorrelation', recInfo);
    HMM_crossCorrelation(recInfo, areaInfo, analysisParam, resultSave)
end

%%% ------------------------------------------------------------------- %%%
%%% HMM - data
if Job.HMM_RT % special case
    [analysisParam] = specifyAnalysisParam('HMM_RT', recInfo);
    HMM_RT(recInfo, areaInfo, analysisParam, resultSave)
end

if Job.HMM_microsaccades 
    [analysisParam] = specifyAnalysisParam('HMM_microsaccades', recInfo);
    HMM_data(recInfo, areaInfo, analysisParam, resultSave)
end

if Job.HMM_pupil
    [analysisParam] = specifyAnalysisParam('HMM_pupil', recInfo);
    HMM_data(recInfo, areaInfo, analysisParam, resultSave)
end

if Job.HMM_rate
    [analysisParam] = specifyAnalysisParam('HMM_rate', recInfo);
    HMM_data(recInfo, areaInfo, analysisParam, resultSave)
end

if Job.HMM_transitionTime
    [analysisParam] = specifyAnalysisParam('HMM_transitionTime', recInfo);
    HMM_analyseStatistics(recInfo, areaInfo, analysisParam, resultSave)
end

if Job.HMM_transitionTriggeredAverage
    [analysisParam] = specifyAnalysisParam('HMM_transitionTriggeredAverage', recInfo);
    HMM_data(recInfo, areaInfo, analysisParam, resultSave)
end

if Job.HMM_spectrogram
    [analysisParam] = specifyAnalysisParam('HMM_spectrogram', recInfo);
    HMM_data(recInfo, areaInfo, analysisParam, resultSave)
end

% 
tElapsed = toc(tStart);
fprintf('TIC TOC: %g\n', tElapsed);
