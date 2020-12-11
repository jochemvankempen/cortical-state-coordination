function [STA,BinCentres,SpikeNum,SP,CentreBinNr] = CC_SpikeSpikeTA(s, sTrig, binSize, maxTimeLag, timeEpoch) 
% [STA,BinCentres,SpikeNum,SP,CentreBinNr] = CC_SpikeSpikeTA(s, sTrig, binSize, maxTimeLag, timeEpoch) 
% 
% spike triggered average of a spiking signal
% Based on code by Alwin Gieselmann 
%
% Example
% -------
% ::
% 
%     [STA,BinCentres,SpikeNum,SP,CentreBinNr] = CC_SpikeSpikeTA(analogSignal, sTrig, 0.01, 0.1, [0.4 1.2]) 
% 
% Parameters
% ----------
% s : cell
%     cell array of size (numTrial, 1) with spike times relative to event
% sTrig : cell
%     cell array of size (numTrial, 1) with spike times relative to event
% binSize : float
%     width of the CCH bins
% maxTimeLag : float
%     float to indicate the time lag
% timeEpoch : array
%     array of size (1, 2) with times between which to consider spikes
% 
% Returns
% -------
% STA : array
%     array of size (numTrials x numBins). Spike triggered average cross
%     correlation histogram 
% BinCentres : array
%     array of size (), centres of time bins
% SpikeNum : array
%     array of size (numTrial x 1) double indicating the number of spikes
% SP : array 
%     array of size (numTrials x numBins). shift/shuffle predictor cross
%     correlation histogram  
%
% 


%% get trial groups
nTr = size(s,1);
TrInd = [1:nTr]';
%TrIndShift = circshift([1:nTr]',[-1]);
TrIndShift = randperm(nTr)';

%% settings
BinCentres = [];
BinCentres = [maxTimeLag*(-1):binSize:maxTimeLag];
CentreBinNr = find(BinCentres==0);
BinIndex = [1:size(BinCentres,2)] - CentreBinNr;

%% initialise
STA = NaN(nTr, length(BinCentres));
SP = NaN(nTr, length(BinCentres));
SpikeNum = NaN(nTr,1);

%% loop
for itrial = 1:nTr
        
    cSpkTrig = sTrig{itrial}(:);
    cSpkTrig = cSpkTrig(cSpkTrig(:,1)>=timeEpoch(1) & cSpkTrig(:,1)<=timeEpoch(2));
    
    cSpkSig = s{itrial}(:);
    cSpkSigShifted = s{TrIndShift(itrial)}(:);
    
    if isempty(cSpkSig) || isempty(cSpkSigShifted) || isempty(cSpkTrig)
%         keyboard
    end
    
    [STA(itrial,:),~,SpikeNum(itrial)]  = ComputeCCH( cSpkSig,         cSpkTrig, timeEpoch, timeEpoch, BinCentres );
    [SP(itrial,:),~,~]                  = ComputeCCH( cSpkSigShifted,  cSpkTrig, timeEpoch, timeEpoch, BinCentres );
    
end

function [F, n1, n2] = ComputeCCH(S1,S2,TW1,TW2,BinCentres)
% Times are spikes in S1 relative to S2,
% i.e. times are S1-S2, mind the sign of subtraction
% i.e. CCH zeros are spikes in S2
% S1, S2 .... timestamps
% TW1,TW2 ... time windows
% histBinEdges ... 

n1=0; n2=0;

try
    i1 = S1(:,1)>=TW1(1) & S1(:,1)<=TW1(2);
    i2 = S2(:,1)>=TW2(1) & S2(:,1)<=TW2(2);
    n1 = sum(i1);
    n2 = sum(i2);
catch
%     warning('no spikes to select')
end

binWidth = median(diff(BinCentres));
% add flanking bins to make sure that counts outside bins are neglected
BinCentres = [BinCentres(1)-binWidth BinCentres BinCentres(end)+binWidth];

if any([n1 n2]==0)
    F = zeros(1,length(BinCentres));
else
    D = repmat(S1(i1),1,sum(i2)) -repmat(S2(i2)',sum(i1),1);
    F = hist(D(:),BinCentres);
end
F = F(2:end-1);

