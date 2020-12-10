function [recInfo] = getRecordingInfo(recordingList, varargin)
%
% select information from a specific recording
%
% Example
% -------
% ::
% 
%     [recInfo] = getRecordingInfo(recordingList, varargin)
%     [recInfo] = getRecordingInfo(recordingList, 'Subject', 'T', 'Date', '2017-03-28')
%     [recInfo] = getRecordingInfo(recordingList, 'RecIdx', 1)
% 
% Parameters
% ----------
% recordingList : table
%     table with information about recordings, acquired via _getRecordingList 
% varargin : cell
%     criteria for selection of recording info:
%     
%     - e.g. 'RecIdx', irec
%     - e.g. 'Subject', 'T', 'Date', '2017-03-28'
%
% Return
% ------
% recInfo: table (single row)
%     table with information about a single recording
%


%%% select recording based on criteria
switch varargin{1}
    case 'RecIdx'       
        recordingList = recordingList(varargin{2},:);
    case 'Subject'
        Subject = varargin{2};
        
        recordingList = recordingList( strcmpi(recordingList.Subject,Subject) ,:); % limit recordingList to specific Subject
        
        switch varargin{3}
            case 'RecIdx'
                recordingList = recordingList(varargin{4},:);
            case 'Date'
                recordingList = recordingList( recordingList.Date==varargin{4} ,:);
        end
    case 'HPC'
        listIdx = varargin{2};
        recordingList = recordingList(listIdx.irec,:);
        
        recordingList.HPC = listIdx;
end

if height(recordingList)~=1
    error('You need to select a single recording')
end
recInfo = recordingList;
clear recordingList
