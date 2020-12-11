function CC_HMM_MS_cpt(msu, label, timeEpoch, analysisParam, resultSave)
% CC_HMM_MS_cpt(msu, label, timeEpoch, analysisParam, resultSave)
% 
% Cross correlation histogram between HMM states and microsaccades, coincidences per trigger
% Outline of script is based on jitterPSTHAnalysis.m from Tatiana Engel
%
% Parameters
% ----------
% msu : cell 
%     cell array of size (numChannel, numTrial) which contains transition
%     and microsaccade start times 
% label : cell
%     cell array with label (string) for each channel
% timeEpoch : array
%     array [t_begin t_end] time epoch for correlations computation
% analysisParam : struct
%     struct, fields contain analysis parameters: binSize, maxTimeLag
% resultSave : struct
%     struct that contains information relevant for saving and plotting the result
%
% Returns
% -------
% crossCorr : array
%     array of size (numChannel*(numChannel-1)/2, numBin), cross-correlation function for each channel
% meanRateDU : array 
%     array of size (numChannel, 1), mean transition rate 
% meanRateMS : array 
%     array of size (numChannel, 1), mean microsaccade rate 
% timeBinsLag : array
%     time lag axis
% labelCombinations : cell
%     cell with string combinations of labels for trigger and triggered channel
%
%
% **HMM_MS_cpt.mat : file**
%     file that includes all the variables listed under *Returns*
% 


    if size(timeEpoch,1)>1
        timeEpoch = timeEpoch(1,:);
    end
    
    binSize = analysisParam.binSize;
    maxTimeLag = analysisParam.maxTimeLag;
    numTrial = size(msu,2);

    timeBinsLag = -maxTimeLag:binSize:maxTimeLag;

    numChannel = size(msu,1)-1;
    numBin = length(timeBinsLag);

    crossCorr = zeros(numChannel, numBin);
    
    meanRateDU = zeros(numChannel, 1);
    
    labelCombinations = cell(numChannel, 1);

    % compute cross-correlations
    for iChannel = 1:numChannel
        pickChannel = [numChannel+1 iChannel]; % the last channel (microsaccades) is always trigger.
        
        [STA,timeBinsLag,SpikeNum,SP,CentreBinNr] = CC_SpikeSpikeTA(msu(pickChannel(2),:)', msu(pickChannel(1),:)', binSize, maxTimeLag, timeEpoch) ;
        
        meanRate_iChan = sum( cellfun(@(x) length(x(x>=timeEpoch(1) & x<=timeEpoch(2))), msu(pickChannel(1),:) ) ) ...
            /(size(msu,2)*(timeEpoch(2)-timeEpoch(1)));
        meanRate_jChan = sum( cellfun(@(x) length(x(x>=timeEpoch(1) & x<=timeEpoch(2))), msu(pickChannel(2),:) ) ) ...
            /(size(msu,2)*(timeEpoch(2)-timeEpoch(1)));
        
        % normalize to Coincidences/spike
        %STA=STA*(1/numTrial)*(1/sqrt(meanRate_iChan*meanRate_jChan));
        %SP=SP*(1/numTrial)*(1/sqrt(meanRate_iChan*meanRate_jChan));

        % normalize to Coincidences/trigger
        STA=STA*(1/sum(SpikeNum));
        SP=SP*(1/sum(SpikeNum));
        
        % correct shift predictor
        crossCorr(iChannel,:) = sum(STA,1) - sum(SP,1);

        meanRateDU(iChannel) = sum( cellfun(@(x) length(x(x>=timeEpoch(1) & x<=timeEpoch(2))), msu(iChannel,:) ) ) ...
            /(size(msu,2)*(timeEpoch(2)-timeEpoch(1)));

        labelCombinations{iChannel,1} = [label{iChannel}];

    end;
    
    
    meanRateMS = sum( cellfun(@(x) length(x(x>=timeEpoch(1) & x<=timeEpoch(2))), msu(numChannel+1,:) ) ) ...
        /(size(msu,2)*(timeEpoch(2)-timeEpoch(1)));

    dirName = resultSave.dataDirName;
    % save data
    if (resultSave.data)
        if (~isdir(dirName))
            mkdir(dirName);
        end;
        fileName = fullfile(dirName, 'HMM_MS_cpt.mat');

        % save cross-correlations        
        save(fileName, 'crossCorr', ...
            'meanRateDU', 'meanRateMS', 'timeBinsLag', 'labelCombinations');
        
    end;

end










