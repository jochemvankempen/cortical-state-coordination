function [ timeEpochState, durationState, numTransition ] = extractUpDownTimeEpoch( states, timeBins, numState, binSize, inclBorder)
%  [ timeEpochState, durationState, numTransition ] = extractUpDownTimeEpoch( states, timeBins, numState, inclBorder )
%
% extracts time epochs spent in each state in a sequence of states
%
% arguments:
% ---------
%   states      -   cell array of length numTrial containing sequence of states for a set of trials
%   timeBins    -   sequence of corresponding time bins for each trials, cell array
%   numState    -   number of states in the model
%   binSize     -   size of the bin used to fit the model
%   inclBorder  -   boolean, true if the time epochs at the beginning and end
%                   of each sequence (before the first transition and after the last
%                   transition) should be included
%
% outputs:
% -------
%   timeEpochState  -   cell array of size (numTrial, numState), contains the
%                       time epochs spent in each state on each trial
%   durationState   -   array of size (numTrial, numState), contains total
%                       duration spent at each state on each trial
%   numTransition   -   array of size (numTrial,numState), contains the total number of
%                       transitions in+out of each state on each trial
    
    numTrial = size(states,1);
    % determine unique states
    %{
    maxLength = max( cellfun(@(x) length(unique(x)), states ) );
    temp = cell2mat( cellfun(@(x) [unique(x), NaN(1,maxLength-length(unique(x))) ], states,'UniformOutput',0 ) );
    uniqueState = sort( unique( temp(~isnan(temp)) ) );
    clear temp;
    numState = length( uniqueState );
    %}
    
    timeEpochState = cell(numTrial, numState);
    durationState = zeros(numTrial, numState);
    numTransition = zeros(numTrial, numState);
    
    % adjust time bins to cover the whole time epoch
    %binSize = timeBins{1}(2)-timeBins{1}(1);
    for iTrial = 1:numTrial
        timeBins{iTrial}(end+1) = timeBins{iTrial}(end)+binSize ; % add one more time bin
        timeBins{iTrial} = timeBins{iTrial} - 0.5*binSize; % shift by half a bin
    end;
    
    for iState = 1:numState
        thisState = iState; %uniqueState(iState);     
        for iTrial = 1:numTrial
            isThisState = states{iTrial}==thisState;
            % +1 -> entering state; -1 -> leaving state
            stateChange = [0 isThisState(2:end)-isThisState(1:end-1) 0];

            if (any(stateChange))
                % if first change is -1, change the first element to 1 if
                % inclBorder = true, otherwise change -1 to 0 (i.e. discard the border timeEpoch)
                k = find(stateChange,1,'first');
                if ( stateChange(k)==-1 )
                    if (inclBorder)
                        stateChange(1) = 1;
                    else
                        stateChange(k)= 0;
                    end;
                end;
                % if last change is 1, change the last element to -1 if
                % inclBorder = true, otherwise change 1 to 0 (i.e. discard the border timeEpoch)
                k = find(stateChange,1,'last');
                if ( stateChange(k)==1 )
                    if (inclBorder)
                        if (k==length(stateChange))
                            stateChange(end) = 0;
                        else
                            stateChange(end) = -1;
                        end;
                    else
                        stateChange(k) = 0;
                    end;
                end;
                assert(sum(stateChange==1)==sum(stateChange==-1),'Number of entering and leaving state transitions are unequal');
                if (any(stateChange))
                    timeEpochState{iTrial,iState} = [timeBins{iTrial}(stateChange==1)', timeBins{iTrial}(stateChange==-1)'];
                end;
            elseif (all(isThisState)) % all the time in this state
                timeEpochState{iTrial,iState} = [timeBins{iTrial}(1) timeBins{iTrial}(end)];
            end;
            
            %numTransition(iTrial, iState) = sum(abs(stateChange));
            numTransition(iTrial, iState) = sum(stateChange==1);
            
            % determine duration in this state on this trial
            if (isempty(timeEpochState{iTrial,iState}))
                durationState(iTrial,iState) = 0;
            else
                durationState(iTrial,iState) = sum( timeEpochState{iTrial,iState}(:,2)- timeEpochState{iTrial,iState}(:,1));
            end;
        end;   
    end;


end

