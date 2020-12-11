function perform_pooling(recordingList,Job,readDataFromServer,saveDataToServer)
% perform_pooling(recordingList,Job,readDataFromServer,saveDataToServer)
%
% Wrapper function to call `pool_data.m` 
%
% Parameters
% ----------
% recordingList : table
%     table with information about recordings, acquired via _getRecordingList 
% Job : struct
%     struct with fields to specify which analyses to run
% readDataFromServer : boolean
%     boolean indicating whether to read data from server location (see
%     computerSetup.m) 
% saveDataToServer : boolean
%     boolean indicating whether to save data to server location (see
%     computerSetup.m) 
% 



allJobs = fields(Job);

for ijob = 1:length(allJobs)
    
    if ~Job.(allJobs{ijob})
        continue
    end
    
    switch allJobs{ijob}
        case 'HMM_fit'
            pool_data(recordingList, 'HMM', readDataFromServer, saveDataToServer);
            
        case 'HMM_modelSelection'
            pool_data(recordingList, 'HMMselect', readDataFromServer, saveDataToServer);
            
        otherwise
            
            pool_data(recordingList, allJobs{ijob}, readDataFromServer, saveDataToServer);
    end
end




