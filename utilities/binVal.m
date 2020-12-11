function [binEdges, binIdx, binValues] = binVal(inpData, nBin, binType, edges)
% [binEdges, binIdx, binValues] = binVal(inpData, nBin, binType, edges)
%
% sort values (e.g. RTs) in bins. Ignores NaN values (note that zeros are
% used in the binning!)
%
% Parameters
% ----------
% inpData : array
%     values to sort
% nBin : double
%     number of bins to sort it in
% binType : string
%     'manual', 'median' or 'equal' (default) split. If manual, edges is required
% edges : array
%     the edges for the binning procedure in case manual is selected
%
% Returns
% -------
% binEdges : array
%     Edges of the bins
% binIdx : array
%     idices of bin for every original value
% binValues : cell
%     original values, binned
%
% Jochem van Kempen, 07/02/2017


[tmpDatOriginal,sortIdxOriginal] = sort(inpData(:)); % sort all values

nanIdx          = isnan(tmpDatOriginal);
tmpDat          = tmpDatOriginal(~nanIdx);
% sortIdx         = sortIdxOriginal(~nanIdx);
if nargin<3
    binType = 'equal';
end

if strcmpi(binType,'manual')
    if nargin<4
        error('Need to define edges when binning values manually')
    end
end

if nBin == 1
    binEdges    = [tmpDat(1)-1 tmpDat(end)+1];
    binIdx      = ones(length(tmpDat),1);
    binValues{1}= tmpDat;
    return
end

% [~, revIdx] = sort(sortIdxOriginal);%unsort your data

% if sum(tmpDatOriginal(revIdx) == inpData(:)) ~= length(inpData)
%     error('something went wrong sorting')
% end      

switch binType
    case 'median'
        if nBin ~= 2
            error('cannot do median split with number of bins unequal to two')
        end
        binEdges = [tmpDat(1) median(tmpDat) tmpDat(end)];
        
    case 'equal'
        binSize    = floor(length(tmpDat)/nBin);
        for iBin = 1:nBin-1
            binEdges(iBin) = tmpDat(binSize*iBin);
        end
        binEdges = [tmpDat(1) binEdges tmpDat(end)];
    case 'manual'
        binEdges = edges;
end

if nargout == 1 % no need for 
    binIdx = [];
    binValues = [];
    return
end

binIdx      = zeros(length(inpData),1);
binValues   = cell(nBin,1);

for iBin = 1:nBin
    if iBin == 1
        binIdx(inpData >= binEdges(iBin) & inpData <= binEdges(iBin+1))=iBin;
    elseif iBin == nBin
        binIdx(inpData > binEdges(iBin) & inpData <= binEdges(iBin+1))=iBin;
    else
        binIdx(inpData > binEdges(iBin) & inpData <= binEdges(iBin+1))=iBin;
    end
    binValues{iBin} = [inpData(binIdx==iBin)];
end

if ~isempty(isnan(inpData))
    binIdx(isnan(inpData)) = NaN;
end
    
    
    
    
    