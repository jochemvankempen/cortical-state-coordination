function [STA,BinCentres,SpikeNum,SP,CentreBinNr] = CC_SpikeTA(analogSignal,sTrig,maxTimeLag,timeEpoch,WinExtendFlag,trialList) 
% [STA,BinCentres,SpikeNum,SP,CentreBinNr] = CC_SpikeTA(analogSignal,sTrig,maxTimeLag,timeEpoch,WinExtendFlag) 
% 
% spike triggered average
% Based on code by Alwin Gieselmann 
%
% Example
% -------
% ::
% 
%     [STA,BinCentres,SpikeNum,SP,CentreBinNr] = CC_SpikeTA(analogSignal,sTrig, 0.1, [0.4 1.2], 0) 
% 
% Parameters
% ----------
% analogSignal : cell
%     cell array of size (numTrial, 1) with spike times relative to event
% sTrig : cell
%     cell array of size (numTrial, 1) with spike times relative to event
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

if nargin<5
    WinExtendFlag = false;
end
if nargin<6
    trialList=[];
end

%% get trial groups
nTr = size(sTrig,1);
TrInd = [1:nTr]';
if ~isempty(trialList)
    
    % get unique trial IDs
    TrIndList = unique(trialList);
    
    % permute trial index
    TrIndPerm = randperm(length(TrIndList));
    
    % create LUT between trial IDs and indices
    TrIndListPerm = NaN(max(trialList),1);
    TrIndListPerm(TrIndList) = TrIndPerm;
    
    % get trial shift indices
    TrIndShift = TrIndListPerm(trialList);
else
    %TrIndShift = circshift([1:nTr]',[-1]);
    TrIndShift = randperm(nTr)';
end

%% check input

assert( size(timeEpoch,1)==nTr, 'timeEpoch does not have correct number of trials')

%% set timing
t = analogSignal.TimeStamps;
t = t(:)';% make sure t is horizontal
Fs = analogSignal.Fs;

BinCentres = [fliplr(1000/Fs:1000/Fs:(maxTimeLag*1000)).*(-1) 0 [1000/Fs:1000/Fs:(maxTimeLag*1000)]];
BinIndex = round(BinCentres./(1000/Fs));
CentreBinNr = size(BinCentres,2)/2+1;
        
STA = NaN(nTr, length(BinCentres));
SpikeNum = NaN(nTr,1);
SP = NaN(nTr, length(BinCentres));

%% loop
for itrial = 1:nTr
    
    cSpkTrig = sTrig{itrial}(:);
    cSpkTrig = cSpkTrig(cSpkTrig(:,1)>=timeEpoch(itrial,1) & cSpkTrig(:,1)<=timeEpoch(itrial,2));

    if WinExtendFlag
        binWin = [];
    else
        binWin = t>=timeEpoch(TrInd(itrial),1) & t<=timeEpoch(TrInd(itrial),2);
    end
    
    [STA(itrial,:),SpikeNum(itrial)] = ComputeSTA( ...
        cSpkTrig, ...
        squeeze(analogSignal.Samples(1,TrInd(itrial),:)), ...
        t,BinIndex,binWin);
    
    SP(itrial,:) = ComputeSTA( ...
        cSpkTrig, ...
        squeeze(analogSignal.Samples(1,TrIndShift(itrial),:)), ...
        t,BinIndex,binWin);

end

function [F,nSpk] = ComputeSTA(spkSig,anaSig,t,binIndex,binWindow)
nB = length(binIndex);
nSpk = length(spkSig);
if isempty(spkSig)
    F = zeros(1,nB).*NaN;
    return;
end
anaSig = anaSig(:)';% make sure anaSig is horizontal
nT = length(t);
dT = repmat(t,nSpk,1) - repmat(spkSig,1,nT);
[minBinDiff,ZeroBinIndex] = min(abs(dT),[],2);
binIndex = repmat(binIndex,nSpk,1) + repmat(ZeroBinIndex,1,nB);
iS2 = sub2ind([1,nT],ones(nSpk,nB),binIndex);

if ~isempty(binWindow)
    anaSig(~binWindow) = NaN;
end
F = sum(anaSig(iS2),1)./sum(~isnan(anaSig(iS2)),1);

