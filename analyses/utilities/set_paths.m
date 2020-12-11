function [recInfo, paths] = set_paths(recInfo,task,dataFolder,readDataFromServer,saveDataToServer)
% [paths] = setPaths(recInfo,task,processingStage,readDataFromServer,saveDataToServer)
% 
% Create folder structure, in the form of:
% read data : [drive path]/task/data/dataFolder/subject/rec/fig
% store results : [drive path]/task/data/analyzed/ANALYSISNAME/subject/rec/fig
%
% Note: ANALYSISNAME will later be replaced with the performed analysis
%
% This function calls 'computer_setup' to set the machine specific
% environment
% 
% Example
% -------
% ::
% 
%     [paths] = setPaths(recInfo,'gratc','processed',1,1)
%     load_data_path = [paths.read paths.rec];
%
%
% Parameters
% ----------
% recInfo : table
%     table (single row) with information about recordings
% task : string
%     string with task name (e.g. gratc), optional
% dataFolder : string
%     string with folder name (e.g. 'processed'), optional
% readDataFromServer : boolean
%     boolean indicating whether to read data from server location (see
%     computer_setup.m), default is 0
% saveDataToServer : boolean
%     boolean indicating whether to save data to server location (see
%     computer_setup.m), default is 0
%
% Returns
% -------
% recInfo : ammended with paths.read and paths.target
% paths : struct
%     structure with fields
%
%     * base : local drive
%     * server : server location
%     * read : drive where data will be read from
%     * target : drive where data will be stored to
%     * data : data folder
%     * subject : subject folder
%     * rec : recording folder
%     * fig : figure folder
%
% 

%% check input
assert(height(recInfo)==1, 'height recInfo ~= 1')

if isempty(task)
    % warning('task is empty, data will be read from/saved to [drive path]')
end
if isempty(dataFolder)
    % warning('processingStage is empty, data will be read from/saved to [drive path]')
end
if nargin<4
    readDataFromServer = 0;
end
if nargin<5
    saveDataToServer = 0;
end

%% set paths
[paths] = computer_setup;
paths.dataFolder = dataFolder;

% set read path to local or server location
if readDataFromServer
    paths.read = paths.server;
elseif ~readDataFromServer
    paths.read = paths.base;
end
if saveDataToServer
    paths.target = paths.server;
elseif ~saveDataToServer
    paths.target = paths.base;
end

% create folder structure for specific dataFolder and 'analysed' folder
paths.data = [char(task) filesep 'data' filesep ];
paths.readData = [paths.data dataFolder filesep];
paths.analysis = [paths.data 'analysed' filesep 'ANALYSISNAME' filesep];% ANALYSISNAME will be replaced later

% create folder structure for subject, pen and fig
paths.subject = [recInfo.Subject{1} filesep];
paths.rec = [paths.subject char(recInfo.Date(1)) filesep];
paths.fig = [paths.rec 'figures' filesep];

% set these paths in recInfo for easy access
recInfo.path_read = fullfile(paths.read, paths.readData, paths.rec);
recInfo.path_target = fullfile(paths.target, paths.analysis, paths.rec);

assert(isfolder(recInfo.path_read), 'read path does not exist')


