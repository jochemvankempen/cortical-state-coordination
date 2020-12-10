function [paradigmInfo] = getParadigmInfo(recInfo, varargin)
%
% get info about paradigms that were run for a specific recording. 
%
% compile paradigmlist from csv file thiele.Paradigm_[Subject].csv
% path of lists is defined in _getListPath
%
% Example
% -------
% ::
%     [paradigmInfo] = getParadigmInfo(recInfo, varargin)
%     [paradigmInfo] = getParadigmInfo(recInfo, 'Paradigm', 'Gratc')
% 
% Parameters
% ----------
% recordingList : table
%     table with information about recordings, acquired via _getRecordingList 
% varargin : cell
%     criteria for selection of penetration info:
%     
%     - e.g. 'Paradigm', 'Gratc'
%
% Return
% ------
% paradigmInfo: table
%     table with information about paradigms for specific recording
%
% Notes
% -----
% - Created, Jochem van Kempen, 16-03-2020
%

assert(height(recInfo)==1, 'height recInfo ~= 1')

path_list = getListPath(recInfo.Subject{1});

csvdata_paradigm = readtable(fullfile(path_list, sprintf('thiele.Paradigm_%s.csv', recInfo.Subject{1})), 'TextType', 'string');

% reformat date
csvdata_paradigm.Date = datetime(csvdata_paradigm.Date,'Format','yyyy-MM-dd');

% select the indices for this recording
idx_rec = ( csvdata_paradigm.Date == recInfo.Date );

paradigmInfo = csvdata_paradigm(idx_rec,:);

%%% select specific paradigms based on criteria
switch varargin{1}
    case 'Paradigm'       
        idx = strcmpi(paradigmInfo.Paradigm, varargin{2});
        paradigmInfo = paradigmInfo(idx,:);
    otherwise
        error('not coded yet')
end

% convert cells
% loop over all columns in table
columnNames = paradigmInfo.Properties.VariableNames;
for icol = 1:length(columnNames)
    try
        paradigmInfo.(columnNames{icol}) = eval(paradigmInfo.(columnNames{icol}){1});
    catch
    end
end


