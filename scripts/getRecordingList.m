function [recordingList] = getRecordingList(Subjects, varargin)
% [recordingList] = getRecordingList(Subjects, varargin)
% 
% compile recordinglist from csv file thiele.Session_[Subject].csv
% path of lists is defined in _getListPath
%
% Example
% -------
% ::
% 
%     [recordingList] = getRecordingList(Subjects, varargin)
%     [recordingList] = getRecordingList({'T'}, 'task', 'gratc', 'recType', 'Pen')
%     [recordingList] = getRecordingList({'T'}, 'recordings', { {'2016-07-13','2016-07-14'} })
%
%
% Parameters
% ----------
% Subjects : cell
%     cell array with subject names
% varargin : cell
%     select recordings based on 
%
%     - 'task': select all recordings with a certain task
%     - 'recordings': select specific recordings based on date
%     - 'recType': 'Pen' (neural recordings), 'Beh' (Behavioural recordings)
%
% Returns
% -------
% recordingList: table
%     table with information about recordings. 
%

if any( strcmpi(varargin, 'recType') )
    idx = find(strcmpi(varargin, 'recType'));
    recType = varargin{idx+1};
else
    recType = [];
end

recordingList = [];

for iSubject = 1:length(Subjects)
    
    path_list = getListPath(Subjects{iSubject});
    
    csvdata_session = readtable(fullfile(path_list, sprintf('Thiele.Session_%s.csv', Subjects{iSubject})), 'TextType', 'string');
      
    csvdata_session{:,'Subject'} = Subjects(iSubject);
    csvdata_session{:,'iSubject'} = iSubject;
       
    % select sessions 
    if (~isempty(varargin))
        if strcmpi(varargin{1}, 'recordings')
            idx_session = ismember( csvdata_session.Date, varargin{2}{iSubject} );
        elseif strcmpi(varargin{1}, 'task')
            idx_session = strcmpi( csvdata_session.Task, varargin{2} );
        else
            error('Not coded: selection criterium %s', varargin{2})
        end
    else
        idx_session = true( height(csvdata_session), 1);
    end

    if isempty(csvdata_session(idx_session,:))
        fprintf('No recordings found for %s\n', Subjects{iSubject})
        continue
    end
    
    % concatenate multiple subjects
    if (iSubject==1)
        recordingList = table2struct( csvdata_session(idx_session,:) );
    else
        recordingList = [recordingList; table2struct( csvdata_session(idx_session,:) )];
    end

end
recordingList = struct2table(recordingList);

recordingList = [recordingList(:,end-1:end) recordingList(:,1:end-2)];

% reformat date
recordingList.Date = datetime(recordingList.Date,'Format','yyyy-MM-dd');

if ~isempty(recType)
    switch recType
        case 'Pen'
            cond1 = ~isnan(recordingList.PenID);
            cond2 = ~(recordingList.PenID==0);
            idx = cond1 & cond2;
            recordingList = recordingList(idx,:);
        case 'Beh'
    end
end
    
