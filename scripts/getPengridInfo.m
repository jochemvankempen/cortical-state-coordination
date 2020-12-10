function [pengridInfo] = getPengridInfo(recInfo)
%
% get info from thiele.Pengrid_[Subject].csv
% path of lists is defined in _getListPath
%
% Example
% -------
% ::
% 
%     [pengridInfo] = getPengridInfo(recInfo)
%
%
% Parameters
% ----------
% recInfo : table
%     table (single row) with information about recordings
%
% Returns
% -------
% pengridInfo : struct
%     structure with information about recording parameters (channel info
%     etc) for each area 
% 

assert(height(recInfo)==1, 'height recInfo ~= 1')

% Get information from Pengrid csv file
try
    path_list = getListPath(recInfo.Subject{1});
    csvdata_pengrid = readtable(fullfile(path_list, sprintf('Thiele.PenGrid_%s.csv', recInfo.Subject{1})), ...
        'TextType', 'string');
catch
    csvdata_pengrid = [];
end

% reformat date
csvdata_pengrid.DATE = datetime(csvdata_pengrid.DATE,'Format','yyyy-MM-dd');

% select the indices for this penetration
idx_pengrid = find( csvdata_pengrid.DATE == recInfo.Date );

nArea = length(idx_pengrid);
for iarea = 1:nArea
    
    % select only one entry (row)
    csvdata_area = csvdata_pengrid(idx_pengrid(iarea),:);
    
    % extract the name of the area
    tmp_idx = strfind(csvdata_area.CHAMBERNAME{1}, '_');
    pengridInfo(iarea).NAME = csvdata_area.CHAMBERNAME{1}(1:(tmp_idx-1));
    
    % loop over all columns in table
    columnNames = csvdata_area.Properties.VariableNames;
    for icol = 1:length(columnNames)
        if isnumeric( csvdata_area.(columnNames{icol})(1) ) || isdatetime( csvdata_area.(columnNames{icol})(1) )
            pengridInfo(iarea).(columnNames{icol}) = csvdata_area.(columnNames{icol})(1);
            
        elseif isstring( csvdata_area.(columnNames{icol})(1) ) || iscell( csvdata_area.(columnNames{icol})(1) )
            
            try 
                pengridInfo(iarea).(columnNames{icol}) = eval(csvdata_area.(columnNames{icol}){1});
            catch
                pengridInfo(iarea).(columnNames{icol}) = csvdata_area.(columnNames{icol}){1};
            end
        else 
            warning('Undefined column in %s', sprintf('thiele.Pengrid_%s.csv', recInfo.Subject{1}) )
        end
    
    end
    
end

        
        



