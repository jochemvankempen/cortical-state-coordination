function outData = selectStabilityWindow(recInfo, inData, inDataType)
% limit data to trial window determined using
% determine_recording_stability.m

readStabilityName = fullfile(recInfo.path_read, 'stability.mat');

load(readStabilityName, 'selectedTrialWindow');

switch inDataType
    case {'hash','unit','MUA_spont20',...
            'NCS','LFP','LFPb','MUAe',...
            'microsaccades',...
            }
        
        allFields = fieldnames(inData);
        allEventIdx = find(contains(allFields, 'Align'));

        outData = inData;
        
        for ievnt = 1:length(allEventIdx)
            switch inDataType
                case {'NCS','LFP','LFPb','MUAe'}
                    outData.(allFields{allEventIdx(ievnt)}).Samples = inData.(allFields{allEventIdx(ievnt)}).Samples(:,selectedTrialWindow(1):selectedTrialWindow(2),:);
                case {'hash','unit','MUA_spont20'}
                    outData.(allFields{allEventIdx(ievnt)}) = inData.(allFields{allEventIdx(ievnt)})(:,selectedTrialWindow(1):selectedTrialWindow(2));
                case {'microsaccades'}
                    outData.(allFields{allEventIdx(ievnt)}).ms = inData.(allFields{allEventIdx(ievnt)}).ms(selectedTrialWindow(1):selectedTrialWindow(2));                    
                    outData.(allFields{allEventIdx(ievnt)}).Samples = inData.(allFields{allEventIdx(ievnt)}).Samples(:,selectedTrialWindow(1):selectedTrialWindow(2),:);
                otherwise
                    error('unknown datatype')
            end
        end
        
    case {'trialData','pupil'}
        outData = inData(selectedTrialWindow(1):selectedTrialWindow(2));
    otherwise
        error('unknown datatype')
end

