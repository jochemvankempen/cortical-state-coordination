function spike_FFanalysis( mu, newTimeEpoch, analysisParam, dataType, msUnitList, areaList, resultSave  )
% perform Fano Factor analysis for a collection of channels, save and plot the
% result
%
% arguments:
% ----------
%   mu              -   cell array of size (numChannel, numTrial) which contains spike
%                       times
%   newTimeEpoch    -   [t_begin t_end] time epoch for rate computation
%   analysisParam   -   structure array that contains analysis parameters,
%                       in particular windowType, single, sliding or
%                       varSize
%   layerData       -   cell array that contains layer information for each
%                       channel
%   dataType        -   cell array that contains data type string for each
%                       channel, 'su' or 'mu'
%   msUnitList      -   cell array that contains number of each channel:
%                       [chanNum] for mu, and [chanNum unitNum] for su
%   resultSave      -   structure that contains information relevant for
%                       saving and plotting the result
% saved data:
% -----------
%   file:   FF.mat
%   ----------------
%       FF          -   array of size (numChannel, numBin), contains
%                       FF of each channel
%       sFF         -   array of size (numChannel, numBin), contains
%                       jack knife error bars for the FF of each channel
%       time        -   array of length numBin, contains time axis for the FF
%       window      -   array window size axis for the FF
%       meanCount   -   mean spike count in each time bin for each window
%                       size (numTimeBin, numWinSize)
%       varCount    -   variance of spike count in each time bin for each
%                       window size (numTimeBin, numWinSize)

    dirName = resultSave.dataDirName;
    fileName = [dirName 'FF.mat'];

    if exist(fileName,'file')
        fprintf('File already exists: %s\n',fileName)
        %return
    end

    numChannel = size(mu,1);
    % get the dimensions for the FF array
    [ time, window, numSubWindow, shiftTime ] = spike_getFFtimeWin( newTimeEpoch, analysisParam );

    FF = zeros(numChannel,length(time),length(window));
    sFF = zeros(numChannel,length(time),length(window));
    meanCount = zeros(numChannel,length(time),length(window));
    varCount = zeros(numChannel,length(time),length(window));

    for iChannel = 1:numChannel
        % compute FF
        [FF(iChannel,:,:), sFF(iChannel,:,:), meanCount(iChannel,:), varCount(iChannel,:)] = spike_FF(  mu, iChannel, time, window, numSubWindow );
        
%         % perform plotting
%         resultSave.subdirName = getSubDirSingle(dataType, iChannel);
%         label = getFileLabel(dataType{iChannel}, msUnitList{iChannel});
%         resultSave.fileLabel = ['_' label];
%         resultSave.channelLabel{1} = getAxisLabel(dataType{iChannel}, msUnitList{iChannel});
%         resultSave.windowType = analysisParam.windowType;
%         plotFF( FF(iChannel,:), sFF(iChannel,:), time+shiftTime, window, resultSave );
% 
    end;




   % save data
    if (resultSave.data)
        if (~isdir(dirName))
            mkdir(dirName);
        end;
        time = time + shiftTime;
        save(fileName, 'FF', 'sFF', 'meanCount', 'varCount', 'time', 'window', 'dataType', 'msUnitList', 'areaList');
    end;


end

