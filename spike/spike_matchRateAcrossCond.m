function [outData] = spike_matchRateAcrossCond( inData, trialIdx, timeEpoch, thisCond )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

    numCond = length(unique(trialIdx));

    allFields = fieldnames(inData);
    allEventIdx = find(contains(allFields, 'Align'));
    
    assert(length(allEventIdx)==1, 'Specify single align field')
        
    numUnit = size(inData.(allFields{allEventIdx}),1);
    numTotalTrial = size(inData.(allFields{allEventIdx}),2);
    
    % calculate mean rate of each unit in all conditions
    meanRate = NaN(numUnit, numCond);
    for icond = 1:numCond
        
        msu = inData.(allFields{allEventIdx})(:, trialIdx==icond);
        
        if (size(timeEpoch,1)==1)
            newTimeEpoch = timeEpoch;
        elseif (size(timeEpoch,1)==numTotalTrial )
            newTimeEpoch = timeEpoch(trialIdx==icond, :);
        end
        
        numTrial = size(msu,2);
        
        meanFR = zeros(numUnit,1);
        if (size(newTimeEpoch,1)==1)
            for iUnit = 1:numUnit
                numSpike = sum( cellfun(@(x) length(x(x>=newTimeEpoch(1) & x<=newTimeEpoch(2))), msu(iUnit,:)) );
                meanFR(iUnit) = numSpike/(numTrial*(newTimeEpoch(2)-newTimeEpoch(1)));
            end;
        elseif (size(newTimeEpoch,1)==numTrial )
            for iUnit = 1:numUnit
                numSpike = sum( cellfun(@(x,y,z) length(x(x>=y & x<=z)), msu(iUnit,:), num2cell(newTimeEpoch(:,1)'), num2cell(newTimeEpoch(:,2)') ));
                norm = sum(newTimeEpoch,1);
                meanFR(iUnit) = numSpike/(norm(2)-norm(1));
            end;
        end;
        meanRate(:,icond) = meanFR;
        
    end;
    
    % calculate the rate which should be matched for each unit
    rateToMatch = min(meanRate,[],2);
  
    % which condition should be extracted, thisCond
    mu = inData.(allFields{allEventIdx})(:, trialIdx==thisCond);
    if (size(timeEpoch,1)==1)
        newTimeEpoch = timeEpoch;
    elseif (size(timeEpoch,1)==numTotalTrial )
        newTimeEpoch = timeEpoch(trialIdx==thisCond, :);
    end
    numTrial = size(mu,2);
    
    msu = cell(numUnit, numTrial);
    % cut mu to the time epoch
    for iTrial = 1:numTrial
        for iUnit = 1:numUnit
            spikes = mu{iUnit,iTrial};
            if (size(newTimeEpoch,1)==1)
                msu{iUnit,iTrial} = spikes( spikes>=newTimeEpoch(1) & spikes<=newTimeEpoch(2) );
            elseif (size(newTimeEpoch,1)==numTrial)
                msu{iUnit,iTrial} = spikes( spikes>=newTimeEpoch(iTrial,1) & spikes<=newTimeEpoch(iTrial,2) );
            end;
        end;
    end;
   
   
    % now match rates by deleting spikes
    for iUnit = 1:numUnit
        % compute mean rate on these trials meanRateThisCond
        if (size(newTimeEpoch,1)==1)
            numSpike = sum( cellfun(@(x) length(x(x>=newTimeEpoch(1) & x<=newTimeEpoch(2))), msu(iUnit,:)) );
            norm = numTrial*(newTimeEpoch(2)-newTimeEpoch(1));
        else
            numSpike = sum( cellfun(@(x,y,z) length(x(x>=y & x<=z)), msu(iUnit,:), ...
                num2cell(newTimeEpoch(:,1)'), num2cell(newTimeEpoch(:,2)') ));
            norm = sum(newTimeEpoch(:,2)-newTimeEpoch(:,1));
        end;
        meanRateThisCond = numSpike/norm;
        
        numSpikeDelete = floor( (meanRateThisCond-rateToMatch(iUnit))*norm );
        temp = (meanRateThisCond-rateToMatch(iUnit))*norm - numSpikeDelete ;
        p = rand;
        if (p<temp)
            numSpikeDelete = numSpikeDelete + 1 ;
        end;
        %assert(numSpikeDelete>=0,'Negative number of spikes to delete encountered');
        
        if (numSpikeDelete>0)
            indDelete = sort(randperm( sum(cellfun(@(x) length(x), msu(iUnit,:))), numSpikeDelete ));
            
            cumSum = 0;
            for iTrial = 1:numTrial
                spikes = msu{iUnit,iTrial};
                indDeleteThisTrial = indDelete(indDelete>cumSum & indDelete<=cumSum+length(spikes))-cumSum;
                cumSum = cumSum + length(spikes);
                spikes(indDeleteThisTrial) = [];
                msu{iUnit,iTrial} = spikes;
            end;
        end;
    end;
    
    % replace align field with msu in the trials belonging to this
    % condition
    outData = inData;
    outData.(allFields{allEventIdx})(:,trialIdx==thisCond) = msu; 

end

