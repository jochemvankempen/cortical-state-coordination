function batch_analyses
%
% perform analyses for electrophysiological data
%
% **Setup**
% 
% - Set your paths in `./analyses/computer_setup.m`
% - Indicate which analyses (jobs) to run below
% - Set the analysis parameters in ./analyses/utilities/specifyAnalysisParam.m
% - Press play :)
%


subjects = {'J','T','W'};

% Set whether to load/save data from server (input args to _set_paths.m)
readDataFromServer = 0; % boolean, 0: read data from local drive; 1: read data from server
saveDataToServer = 0; % boolean, 0: save data to local drive; 1: save data to server

poolDataForPopulation = 0;

%% addpaths

% name the folder that contains this repository (optional). Note: repositories
% hpc and lists should be in the same folder. If not provided,
% it is searched for.
repoPath = '';

if isempty(repoPath)
    fileName = mfilename('fullpath');
    folderName = 'analyses';
    repoIdx = strfind(fileName,folderName);
    repoPath = fileName(1:repoIdx(1)-1);
end

assert(isfolder(fullfile(repoPath, 'analyses')), 'Cannot find repository: analyses')

addpath(genpath( fullfile(repoPath, 'analyses') ))
addpath(genpath( fullfile(repoPath, 'hpc') ))
addpath(genpath( fullfile(repoPath, 'lists') ))

%% check whether we are on the hpc cluster

global HPC
HPC = false;
if isfolder('/mnt/nfs/home/')
    HPC = true;
end

%% get recording list

recordingList = getRecordingList(subjects);
fprintf('RecordingList contains %d recordings\n', height(recordingList) )

%% Set Jobs to run ---------------------------------------------------- %%%

% specifyAnalysisParam: settings of analysis params

% spike
Job.analyse_spike_rate              = 0; % spike_analyses_main         
Job.analyse_spike_FF                = 0; % spike_analyses_main         

% HMM
Job.HMM_fit                         = 1; % HMM_main
Job.HMM_cvStateBased                = 0; % HMM_main
Job.HMM_modelSelection              = 0; % HMM_main
Job.HMM_2area                       = 0; % HMM_main

% HMM-data
Job.HMM_pupil                       = 0; % HMM_data
Job.HMM_microsaccades               = 0; % HMM_data
Job.HMM_rate                        = 0; % HMM_data
Job.HMM_RT                          = 0; % HMM_RT, special case
Job.HMM_transitionTriggeredAverage  = 0; % HMM_data
Job.HMM_spectrogram                 = 0; % HMM_data

% HMM-statistics
Job.HMM_transitionTime              = 0; % HMM_analyseStatistics                   

% HMM-HMM
Job.HMM_crossCorrelation            = 0; % HMM_crossCorrelation

%% run analyses ------------------------------------------------------- %%%

for irec = 1:height(recordingList)

    if HPC
        % if running on HPC, select the arrayID (recording)
        [arrayID] = SLURM_getEnv;
        [recInfo] = getRecordingInfo(recordingList, 'RecIdx', arrayID);
    else
        [recInfo] = getRecordingInfo(recordingList, 'RecIdx', irec);
    end
    
    fprintf([...
        '----------------------------------------------------------------------------------------------------------------------\n',...
        'Running subject: %s, recording: %s\n'...
        '----------------------------------------------------------------------------------------------------------------------\n'],...
        recInfo.Subject{1}, recInfo.Date)
        
    % get paths, store in recInfo
    [recInfo] = set_paths(recInfo,[],'processed',readDataFromServer,saveDataToServer);
    
    % run analyses
    try
        perform_analyses(recInfo,Job);
    catch
        warning('Cannot perform analyses for subject: %s, recording: %s \n', recInfo.Subject{1}, recInfo.Date)
    end
%     
    if HPC
        %if HPC then only do one loop. comment out when running parfor
%         break
    end
end

%% Pool data for population analyses
if poolDataForPopulation && ~HPC
    perform_pooling(recordingList,Job,readDataFromServer,saveDataToServer)
end

