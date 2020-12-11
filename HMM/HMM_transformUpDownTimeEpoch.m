function [ timeEpochState_out ] = HMM_transformUpDownTimeEpoch( timeEpochState_in, numTransition )
%  [ timeEpochState_out ] = HMM_transformUpDownTimeEpoch( timeEpochState_in, numTransition )
%
% converts time epochs from cell array (with multiple entries per trial) to
% an array of size (numEpoch x time). Dim 1 is trial number, Dim 2 is timewin
%
% Parameters
% ----------
% timeEpochState : cell 
%     cell array of size (numTrial, numState), contains the time epochs
%     spent in each state on each trial 
% numTransition : array 
%     array of size (numTrial,numState), contains the total number of
%     transitions in+out of each state on each trial 
%
% Returns
% -------
% timeEpochState : cell 
%     cell array of size (1,numState), contains the time epochs spent in
%     each state 
% 
    
    [numTrial, numState] = size(timeEpochState_in);
    
    timeEpochState_out = cell(1, numState);
    
    for istate = 1:numState
        nepoch = 0;
        for itrial = 1:numTrial
            for itrans = 1:numTransition(itrial,istate)
                nepoch = nepoch+1;
                timeEpochState_out{istate}(nepoch,:) = [itrial, timeEpochState_in{itrial,istate}(itrans,:)];  
            end;
        end;
    end;
end

