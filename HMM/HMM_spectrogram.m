function HMM_spectrogram(analogSignal, states, timeBins, SF, analysisParam, dataType, msUnitList, areaList, resultSave)
% HMM_spectrogram(msu, SF, analysisParam, dataType, msUnitList, areaList, resultSave)
%
% compute spectrogram for each HMM state
%
% Parameters
% ----------
% analogSignal : struct 
%     struct with fields Samples (numChannel x numTrial x numTimeStamps)
%     and TimeStamps relative to event
% states : cell
%     cell array (numTrial x 1) with states as estimated by HMM
% timeBins : cell
%     cell array (numTrial x 1) with timeBins corresponding to states
% SF : float
%     Sampling Frequency
% analysisParam : struct
%     structure, fields contain analysis parameters
% dataType : cell 
%     cell array that contains data type string e.g. 'LFPb', 'MUAe'
% msUnitList : cell 
%     cell array that contains labels of each channel: [chanNum] for hash,
%     and [chanNum unitNum] for su and mu 
% areaList : cell 
%     cell array that contains which area the channel was recorded from
% resultSave : struct
%     structure that contains information relevant for saving and plotting
%     the result 
%
% Returns
% -------
% spg : array
%     array of size (numChannel, frequencies, numState), spectral
%     decomposition per channel and state 
% frequencies : array
%     frequencies
% SF : float
%     Sampling Frequency
% dataType : cell 
%     cell array that contains data type string e.g. 'LFPb', 'MUAe'
% msUnitList : cell 
%     cell array that contains labels of each channel: [chanNum] for hash,
%     and [chanNum unitNum] for su and mu 
% areaList : cell 
%     cell array that contains which area the channel was recorded from
%
% 
% **HMM_spectrogram.mat : file**
%     file that includes all the variables listed under *Returns*
% 
% 

fprintf('\tComputing spectrogram per HMM state\n')

binSize = analysisParam.binSize;

[numChannel, numTrial, numTimePoints] = size(analogSignal.Samples);

% convert timeBin centres to time bins
timeBins = cellfun(@(x)(x-(binSize/2)),timeBins,'UniformOutput',false);
timeBins = cellfun(@(x)([x x(end)+binSize]),timeBins,'UniformOutput',false);

[ timeEpochState, durationState, numTransition ] = extractUpDownTimeEpoch( states, timeBins, analysisParam.numState, binSize, analysisParam.inclBorder);
[ timeEpochState ] = HMM_transformUpDownTimeEpoch( timeEpochState, numTransition );

if ~isfield(analysisParam, 'nfft') || isempty(analysisParam.nfft)
    maxtwin                 = cellfun(@length, timeBins, 'UniformOutput', true);
    maxSegmentLength        = ceil(SF * max(maxtwin) * binSize);
    nfft                    = max(2^(nextpow2(maxSegmentLength)));
else
    nfft = analysisParam.nfft;
end
[frequencies, findx]    = getfgrid(SF, nfft, analysisParam.fpass);

% initialise
spg = NaN(numChannel, length(find(findx)), analysisParam.numState);
spgdb = NaN(numChannel, length(find(findx)), analysisParam.numState);

for istate = 1:analysisParam.numState

    % only use epochs with minimum duration
    epochIdx = ( (timeEpochState{istate}(:,3)-timeEpochState{istate}(:,2)) >= analysisParam.minEpochDurationWin );
    timeEpochState{istate} = timeEpochState{istate}(epochIdx,:);
    
    nEpoch = length(find(epochIdx));
    
    faultyEpoch = []; % keep track of bad segments (e.g. because of LFP saturation)
    epochSpectrum = NaN(numChannel, length(find(findx)), nEpoch);
    epochSpectrumLog = NaN(numChannel, length(find(findx)), nEpoch);
    for iepoch = 1:nEpoch
        
        trialIdx = timeEpochState{istate}(iepoch,1);
        epoch = timeEpochState{istate}(iepoch,2:3);
        
        %%% spectral decomposition of LFP
        timeIdx   = (epoch(1) <= analogSignal.TimeStamps & epoch(2) > analogSignal.TimeStamps);
        
        segment = squeeze(analogSignal.Samples(:, trialIdx, timeIdx))'; % chop out segment, transpose to (samples x channels)
        segment = segment - mean(segment,1); % get rid of DC offset 
        
        % spectral analysis
        tapers = dpsschk(analysisParam.tapers, length(segment), SF); % get tapers
        J = mtfftc(segment, tapers, nfft, SF);
        J = J(findx, :, :);
        %                     segmentSpectrum(:, :, iepoch) = 10*log10(squeeze(sum(conj(J).*J, 2)))';
        epochSpectrum(:, :, iepoch) = (squeeze(sum(conj(J).*J, 2)))';
        epochSpectrumLog(:, :, iepoch) = 10*log10(squeeze(sum(conj(J).*J, 2)))';
        
%         if ~isempty(find(isinf(epochSpectrum(:, :, iepoch)))) || ~isempty(find(any(epochSpectrum(:, :, iepoch)==0)))
        if ~isempty(find(isinf(epochSpectrumLog(:, :, iepoch)))) || ~isempty(find(any(epochSpectrumLog(:, :, iepoch)==0)))
            warning('ZEROS in SPG!!')
            epochSpectrum(:, :, iepoch) = NaN;
            %plot(segment);
            
            faultyEpoch = [faultyEpoch iepoch];
        end;
    end;
    
    spg(:,:,istate) = squeeze(nanmean(epochSpectrum, 3));
    spgdb(:,:,istate) = squeeze(nanmean(epochSpectrumLog, 3));
    
end;


% save data
dirName = resultSave.dataDirName ;
if (resultSave.data)
    if (~isfolder(dirName))
        mkdir(dirName);
    end;
    fileName = [dirName 'HMM_spectrogram.mat'];
    save(fileName, 'spg', 'spgdb', 'frequencies', 'SF', 'dataType', 'msUnitList', 'areaList');
end;







