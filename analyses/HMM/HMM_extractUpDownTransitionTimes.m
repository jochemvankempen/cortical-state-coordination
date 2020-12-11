function [ stateTransitionTimesPerState, label, timeBins ] = HMM_extractUpDownTransitionTimes( states, timeBins, numState, binSize, timeEpoch, minDuration, tableOutput )
%  [ stateTransitionTimesPerState, label ] = HMM_extractUpDownTransitionTimes( states, timeBins, numState, binSize, minDuration )
%
% extracts transition times, currently only for 2-state HMM
%
% Parameters
% ----------
% states : cell 
%     cell array of size (numTrial x 1) containing sequence of states for a
%     set of trials 
% timeBins : cell
%     cell  array of size (numTrial x 1) containing sequence of
%     corresponding time bins for each trials, cell array 
% numState : double
%     number of states in the model
% binSize : float
%     size of the bin used to fit the model
% minDuration : float 
%     the minimum time around a transition [-minDuration minDuration] where
%     no other transitions occur  
%
% Returns
% -------
% stateTransitionTimesPerState : cell 
%     cell array of size (numTrial, numState), contains the transition
%     times towards each state on each trial 
% label : cell
%     cell array of size (1 x numState) with transition labels
% 
    
    % check input 
    if ~exist('timeEpoch','var')
        timeEpoch = [];
    end
    if ~exist('minDuration','var')
        minDuration = 0;
    end
    if ~exist('tableOutput','var')
        tableOutput = 0;
    end
    
    % check for number of states 
    uniqueStates = cellfun(@(X) ( unique(X) ), states, 'UniformOutput', false);
    uniqueStates = cellfun(@(X) ( length(X) ), uniqueStates, 'UniformOutput', true);
    if any(uniqueStates>2) || (numState>2)
        error('Only works with 2-state model')
    end
    
    numTrial = size(states,1);
    
    % adjust time bins to cover the whole time epoch
    %binSize = timeBins{1}(2)-timeBins{1}(1);
    for iTrial = 1:numTrial
        timeBins{iTrial}(end+1) = timeBins{iTrial}(end)+binSize ; % add one more time bin
        timeBins{iTrial} = timeBins{iTrial} - 0.5*binSize; % shift by half a bin
    end;
    
    if isempty(timeEpoch)
        timeEpoch = cellfun(@(x) ([x(1) x(end)]), timeBins, 'UniformOutput', false);
        timeEpoch = cell2mat(timeEpoch);
    elseif (size(timeEpoch,1)==1)
        timeEpoch = repmat(timeEpoch, [numTrial, 1]);
    end
    
    [~,maxtwin]     = max(cellfun(@length, timeBins(:, 1), 'UniformOutput', true),[],1);
    maxTimeBins     = timeBins{maxtwin};
    
    stateTransition            = cellfun(@(X)(find([0 (diff(X) == 1) | (diff(X) == -1)])), states, 'UniformOutput', false); % when did a stateswitch happen
    stateTranistionTimes       = cellfun(@(X)(maxTimeBins(X)), stateTransition, 'UniformOutput', false); % when did a stateswitch happen
    stateTransitionLabel       = cellfun(@(X)( ismember( find([0 (diff(X) == 1) | (diff(X) == -1)]), find([0 (diff(X) == 1)]) ) + 1 ), states, 'UniformOutput', false); % were they down (1) or up (2) switches?
    label = {'UD','DU'};
    
    if (minDuration~=0)
                
        % if (minDuration==0), no need to check epoch lengths
        for itrial = 1:numTrial
            
            trialWin = timeEpoch(itrial,:);
            
            rmIdx = false(1, length(stateTranistionTimes{itrial})); % keep track of which transition indices to remove
            
            for istateChange = 1:length(stateTranistionTimes{itrial})
                
                stateChangeWin = stateTranistionTimes{itrial}(istateChange) + [-minDuration minDuration];
                
                % get timings of previous/next transition
                stateChangePrevNext = NaN(1,2);
                
                % if stateChange is the first of this trial, select the start
                % of the timeEpoch, else select the previous stateChangeTime
                if istateChange==1
                    stateChangePrevNext(1) = trialWin(1);
                else
                    stateChangePrevNext(1) = stateTranistionTimes{itrial}(istateChange-1);
                end
                
                % if stateChange is the last of this trial, select the end
                % of the timeEpoch, else select the next stateChangeTime
                if istateChange==length(stateTranistionTimes{itrial})
                    stateChangePrevNext(2) = trialWin(2);
                else
                    stateChangePrevNext(2) = stateTranistionTimes{itrial}(istateChange+1);
                end
                
                % is window outside of trial start/end?
                if any(stateChangeWin < trialWin(1)) || any(stateChangeWin > trialWin(2))
                    rmIdx(istateChange) = true;
                    continue
                end
                
                % is window outside of prev/next transition time?
                if any(stateChangeWin < stateChangePrevNext(1)) || any(stateChangeWin > stateChangePrevNext(2))
                    rmIdx(istateChange) = true;
                    continue
                end
                
            end
            
            % remove identified indices
            stateTransition{itrial}(rmIdx) = [];
            stateTranistionTimes{itrial}(rmIdx) = [];
            stateTransitionLabel{itrial}(rmIdx) = [];
            
        end
    end
    
    if tableOutput
        
        transitionMat = [];
        for itrial = 1:numTrial
            if ~isempty(stateTranistionTimes{itrial})
                for itrans = 1:length(stateTranistionTimes{itrial})
                    
                    if itrans==1
                        prevTrans = timeEpoch(itrial,1);
                    else
                        prevTrans = stateTranistionTimes{itrial}(itrans-1);
                    end
                    if itrans==length(stateTranistionTimes{itrial})
                        nextTrans = timeEpoch(itrial,2);
                    else
                        nextTrans = stateTranistionTimes{itrial}(itrans+1);
                    end
                    
                    transitionMat = [transitionMat ; [itrial, stateTransitionLabel{itrial}(itrans), stateTranistionTimes{itrial}(itrans), prevTrans, nextTrans]];
                end
            end
        end
        
        stateTransitionTimesPerState = array2table(transitionMat, 'VariableNames', {'trial','UD','transition','previous','next'});
        
    else
        stateTransitionTimesPerState = cell(numTrial,numState);
        for istate = 1:numState
            
            stateLabelIdx = cellfun(@(X)( X==istate ), stateTransitionLabel, 'UniformOutput', false); % were they down (1) or up (2) switches?
            
            for itrial = 1:numTrial
                if any(stateLabelIdx{itrial})
                    stateTransitionTimesPerState{itrial,istate} = stateTranistionTimes{itrial}(stateLabelIdx{itrial});
                end
            end
        end
    end
    
    
end

