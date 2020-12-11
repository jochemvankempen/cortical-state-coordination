function rate = spike_rateAnalysis(  mu, newTimeEpoch, analysisParam, dataType, msUnitList, areaList, resultSave )
% perform spike rate analysis for a collection of channels, save and plot the
% result
%
% arguments:
% ----------
%   mu              -   cell array of size (numChannel, numTrial) which contains spike
%                       times
%   newTimeEpoch    -   [t_begin t_end] time epoch for rate computation
%   analysisParam   -   structure, fields contain analysis parameters
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
%   file:   rate.mat
%   ----------------
%       rate        -   array of size (numChannel, numTimeBin), contains
%                       firing rate of each channel
%       sRate       -   array of size (numChannel, numTimeBin), contains
%                       jack knife error bars for the firing rate of each
%                       channel
%       time        -   array of length numTimeBin, contains time axis for
%                       the firing rate

    if size(newTimeEpoch,1)>1
%         error('%s analysis needs fixed time epoch', analysisParam.analysisType)
    end

    numChannel = size(mu,1);
    
    % get time bins
    [time, window ] = spike_getRateTime( newTimeEpoch, analysisParam );
    rate = zeros(numChannel, length(time));
    sRate = zeros(numChannel, length(time));
    for iChannel = 1:numChannel
        % calculate rate
        [rate(iChannel,:), sRate(iChannel,:)] = spike_rate(  mu, iChannel, time, window );
        
%         % perform plotting
%         resultSave.subdirName = getSubDirSingle(dataType, iChannel);
%         label = getFileLabel(dataType{iChannel}, msUnitList{iChannel});
%         resultSave.fileLabel = ['_' label];
%         resultSave.channelLabel{1} = getAxisLabel(dataType{iChannel}, msUnitList{iChannel});
%         resultSave.windowType = analysisParam.windowType;
%         plotRate( time, rate(iChannel,:), sRate(iChannel,:), resultSave );
        
    end;

    % save data
    dirName = resultSave.dataDirName ;
    if (resultSave.data)
        if (~isdir(dirName))
            mkdir(dirName);
        end;
        fileName = [dirName 'rate.mat'];
        save(fileName, 'rate', 'sRate', 'time', 'dataType', 'msUnitList', 'areaList');
    end;


end

