function [trialdata] = spike_convoluteGaussian(spike_times, time, sigma)
% [trialdata] = spike_convoluteGaussian(spike_times, time, sigma)
% Convolve spike times with Gaussian, with smooting sigma (in ms)
%
% Example
% -------
% ::
% 
%     [trialdata] = spike_convoluteGaussian(spike_times, 0:0.001:1, 0.01)
%     plot(time, trialdata), plots the histogram for the current trial
%
% Parameters
% ----------
% spike_times : array
%     vector with times at which spike occured
% time : array
%     vector of timepoints (ms)
% sigma : float
%     smoothing factor (ms)
%
% Returns
% -------
% trialdata : array
%     convolved spikes
%
%

% convert to row vector
spike_times = spike_times(:)';

trialdata=zeros(1,length(time));
y=ones(1,length(time));
for ispike = 1:length(spike_times)
    spike=(y*(1/(sigma*sqrt(2*pi))).*exp(-(((time-spike_times(ispike)).^2)/(2*sigma^2))));
    trialdata(1,:)=trialdata(1,:)+spike;
end
