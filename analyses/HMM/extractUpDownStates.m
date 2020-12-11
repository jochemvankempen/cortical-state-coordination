function [ tauUp, tauDown, change ] = extractUpDownStates( states, timeEpoch, binSize )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    tauUp = [];
    tauDown = [];

    numState = 2;
    
    change = states(2:end)-states(1:end-1);
    goUp = find(change==1);
    goDown = find(change==-1);

    if (numel(goDown)>1 && numel(goUp)>1)    
        numUpState = sum( goDown>goUp(1) );
        numDownState = sum( goUp>goDown(1) );
        tauUp = [ tauUp, (goDown(end-numUpState+1:end) - goUp(1:numUpState))*binSize];
        tauDown = [ tauDown, (goUp(end-numDownState+1:end) - goDown(1:numDownState))*binSize];
    elseif (numel(goDown)==0 && numel(goUp)==0)
        if (states==1)  % up state
            tauUp = [ tauUp, timeEpoch(2)-timeEpoch(1)];
        else % down state
            tauDown = [ tauDown, timeEpoch(2)-timeEpoch(1)];
        end;
    elseif (numel(goDown)==0 && numel(goUp)==1)
        tauUp = [ tauUp, timeEpoch(2)-timeEpoch(1) - goUp(1)*binSize];
        tauDown = [ tauDown, goUp(1)*binSize ];
    elseif (numel(goDown)==1 && numel(goUp)==0)
        tauUp = [ tauUp, goDown(1)*binSize ];
        tauDown = [ tauDown, timeEpoch(2)-timeEpoch(1) - goDown(1)*binSize];
    end;
    


end

