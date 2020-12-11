function outData = selectTrials(inData, inDataType, trIdx)
% select trials 

% check input
if ~islogical(trIdx)
    error('trIdx needs to be logical vector')
end

switch inDataType
    case {'hash','unit','MUA_spont20',...
            'NCS','LFP','LFPb','MUAe'...
            'microsaccades'}       
        
        allFields = fieldnames(inData);
        allEventIdx = find(contains(allFields, 'Align'));
        
        outData = inData;
        
        for ievnt = 1:length(allEventIdx)
            switch inDataType
                case {'NCS','LFP','LFPb','MUAe'}
                    outData.(allFields{allEventIdx(ievnt)}).Samples = inData.(allFields{allEventIdx(ievnt)}).Samples(:,trIdx,:);
                case {'hash','unit','MUA_spont20'}
                    outData.(allFields{allEventIdx(ievnt)}) = inData.(allFields{allEventIdx(ievnt)})(:,trIdx);
                case {'microsaccades'}
                    outData.(allFields{allEventIdx(ievnt)}).ms = inData.(allFields{allEventIdx(ievnt)}).ms(trIdx);
                    outData.(allFields{allEventIdx(ievnt)}).Samples = inData.(allFields{allEventIdx(ievnt)}).Samples(:,trIdx,:);
                otherwise
                    error('unknown datatype')
            end
        end
        
    case {'trialData','pupil'}
        outData = inData(trIdx);
    case 'HMM'
        outData = inData(trIdx,:);
    otherwise
        error('unknown datatype')
        
end
