function [recInfo, areaInfo, RF] = getRFinfo(recInfo, areaInfo, read_path)
% load RF info

%% check input
assert(height(recInfo)==1, 'height recInfo ~= 1')

%% get RF info
rffilename = fullfile(read_path, 'RF.mat');

if ~exist(rffilename, 'file')
    warning('%s does not exist', rffilename)
    return
end
RF = load(rffilename);

recInfo.distanceRF = RF.distanceRF;
recInfo.overlapRF = RF.overlapRF;
%recInfo.centerRF = RF.centerRF;

for iarea = 1:length(areaInfo)
    areaInfo(iarea).centerRF = RF.centerRF.([areaInfo(iarea).name]);
    
%     if any(isnan(RF.centerRF(iarea,:))) || any((RF.centerRF(iarea,:))==0)
%         keyboard
%     end
end







