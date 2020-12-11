function [xspikes, yspikes, yoffset] = spike_computeRaster( mu, chan, yoffset, timeEpoch, lineLength )
%
% create matrix with raster
%
% Example
% -------
% ::
% 
%     [xspikes, yspikes] = spike_computeRaster(sps, 3, 100, [0.4 1.2], 10);
%     plot( xpsikes(:), yspikes(:) )
% 
% Parameters
% ----------
% mu : cell
%     cell array of size (numChannel x numTrial) with spike times relative
%     to event
% chan : double
%     double indicating which channel to select
% yoffset : double
%     scalar that indicates the offset of the raster (y-axis), default is 0
% timeEpoch : array (optional)
%     array of size (numTrial x 2) with time epochs. If size (1 x 2), the
%     same time is used across trials
% lineLength : float (optional)
%     float indicating the length of the spike
% 
% Returns
% -------
% xspikes : array
%     array of size (3, numSpikes) with x coordinates of spikes
% yspikes : array
%     array of size (3, numSpikes) with y coordinates of spikes
% yoffset : float
%     updated yoffset
% 
% Notes
% -----
% - Created, Jochem van Kempen, 15-04-2020
% 

% check input
if nargin<3
    yoffset=0;
end
if nargin<4
    timeEpoch=[];
end
if nargin<5
    lineLength = 1;
end

[numChan, numTrial] = size(mu);

cond1 = size(timeEpoch,1)==1;
cond2 = numTrial>1;
if cond1 && cond2
    timeEpoch = repmat(timeEpoch, [numTrial, 1]);
end

xspikes = [];
yspikes = [];

for itrial = 1:numTrial
    
    sps = mu{chan, itrial};
    
    if ~isempty(timeEpoch)
        sps( sps < timeEpoch(itrial,1) | sps > timeEpoch(itrial,2) ) = [];
    end
    
    [xtmp, ytmp] = computeRaster(sps, (itrial * lineLength + yoffset), lineLength );
    
    xspikes = [xspikes xtmp];
    yspikes = [yspikes ytmp];
end

yoffset = max(yspikes(:));


function [xspikes, yspikes] = computeRaster(spike_times, y, y_length)
% computes a raster from spike times
% 
% INPUT:
% spike_times: vector with times at which spike occured
% y: y-position, where to plot this trials' spikes (e.g. trial/channel index)
% y_length: The length of the spike to be drawn
% 
% OUTPUT:
% xspikes, yspikes: x and y values of the spikes
%
% plot(xspikes, yspikes) plots the raster for the current trial
%
% Jochem van Kempen, 2017
 
if nargin<3
    y_length = 1;
end

% convert to row vector
spike_times = spike_times(:)';

xspikes = repmat(spike_times,3,1);
yspikes = nan(size(xspikes));
if ~isempty(yspikes)
    yspikes(1,:) = y-y_length;
    yspikes(2,:) = y;
end