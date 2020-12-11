 function [convolvedSpikes]=spike_convoluteEPSP(spiketrain, rise, decay, convolvedSpikeDur, epochlength)
% [convolvedSpikes]=spike_convoluteEPSP(spiketrain, rise, decay, convolvedSpikeDur, epochlength)
% Convolve spike times with EPSP
%
% Parameters
% ----------
% spiketrain : array
%     vector with times at which spike occured
% rise : float
%     float indicating rise time
% decay : float
%     float indicating decay time
% convolvedSpikeDur : int
%     duration of EPSP
% epochlength : int
%     length of epoch
%
% Returns
% -------
% convolvedSpikes : array
%     convolved spikes
%

%  spiketrain=randperm(2000);
%  spiketrain=spiketrain(1:100);
%  rise=6;
%  decay=60;
convterm=[1:1:convolvedSpikeDur];
EPSP=(1-exp(-convterm./rise)).*(exp(-convterm./decay));

% convolvedSpikes=zeros(1, epochlength + length(EPSP));
% for jk=1:length(spiketrain)
%      convolvedSpikes(spiketrain(jk):spiketrain(jk)+299) = convolvedSpikes(spiketrain(jk):spiketrain(jk)+299)+EPSP;
% end
 
convolvedSpikes = zeros(1, epochlength + convolvedSpikeDur);
convolvedSpikes(spiketrain) = 1;
convolvedSpikes = conv(convolvedSpikes, EPSP);
convolvedSpikes = convolvedSpikes(convolvedSpikeDur + 1:epochlength + convolvedSpikeDur);
 %  clf
%  
%  test=subplot(1,1,1);
%  plot(spiketrain,1,'.k');
%  hold on
% 
%  plot([1:length(convolvedSpikes)],convolvedSpikes,'r');
%  
    

