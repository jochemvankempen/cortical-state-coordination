function [areaInfo, warningGiven] = getAreaInfo(recInfo, read_path)
%
% [areaInfo] = getAreaInfo(recInfo, read_path)
%
% get info on which channels are included for each recording. In- or
% exclusion is based on the criteria set in *set_cfg.m*, potentially using
% 
% - recording depth (based on CSD), 
% - recording stability (selected trial window defined in NLX_determineRecordingStability.m) 
% - signal-to-noise ratio (computed in NCS_SNR.m)
%
% 
% Parameters
% ----------
% recInfo : table
%     table with single row (one recording)
% read_path : string
%     string with data path where data is read from
%
% Returns
% -------
% areaInfo : struct
%     struct with info about channel inclusion etc for each area
% warningGiven : boolean
%     boolean indicating whether one or more settings couldn't be set, e.g.
%     because some files could not be loaded
% 
%

% check input
assert(height(recInfo)==1, 'height recInfo ~= 1')

% load info from pengrid.csv
[pengridInfo] = getPengridInfo(recInfo);
cfg = set_cfg(recInfo);

% initialise
areaInfo = struct([]);
warningGiven = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% load recording stability file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loads file stability.mat (created in NLX_determineRecordingStability.m). It checks whether any channels
% should be excluded based on recording stability.
stability_chanExclude = [];
if any( strcmpi(cfg.channelSelectionCriteria, 'stability') )
    loadfilename = fullfile(read_path, 'stability.mat');
    if exist(loadfilename,'file')
        stability_chanExclude = load(loadfilename,'chan2exclude');
        stability_chanExclude = stability_chanExclude.chan2exclude;
    else
        % this file is essential for these criteria, give error when it cannot be found
        warning('cannot find channels to exclude based on recording stability, run NLX_determineRecordingStability.m')
        warningGiven = true;
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loop over areas to set the channel indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iarea = 1:length(pengridInfo)
    areaInfo(iarea).name = pengridInfo(iarea).NAME;
    
    nchan = length(pengridInfo(iarea).PROBE_CONTACTS);
    
    %%% initialise channel indices
    areaInfo(iarea).ChanIdx.all             = pengridInfo(iarea).PROBE_CONTACTS;
    areaInfo(iarea).ChanIdx.include         = true(1,nchan);
    areaInfo(iarea).ChanIdx.include_bip     = false(1,(nchan-1)); % bipolar channels
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% set channels to in/exclude based on recording stability
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if any( strcmpi(cfg.channelSelectionCriteria, 'stability') )
        
        cond1 = ~(ismember(areaInfo(iarea).ChanIdx.all, stability_chanExclude));
        
        areaInfo(iarea).ChanIdx.include = ...
            areaInfo(iarea).ChanIdx.include & cond1;

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% get coordinates and set the limits of cortical gray matter, in
    %%% relation to layer 4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loads file [area]_Coordinates.mat (created in
    % summaryLaminarPlot.m and determineLayer4.m). In these scripts the CSD
    % transform is performed and the coordinates of layer 4 are stored.
    % Channels are excluded based on distance from layer 4.
    
    coordinates_name = fullfile(read_path, sprintf('%s_Coordinates.mat',areaInfo(iarea).name));

    if exist(coordinates_name,'file')
        load(coordinates_name,'Coordinates');
        areaInfo(iarea).Coordinates = Coordinates;
        
        %%% assign layer info to recording depth
        recDepth = areaInfo(iarea).Coordinates(:,1);
        areaInfo(iarea).layerAssignment = cell(nchan,1);
        areaInfo(iarea).layer4 = (areaInfo(iarea).Coordinates(:,2) == 0);
        
        switch areaInfo(iarea).name
            case 'V1'
                areaInfo(iarea).layerAssignment(recDepth > 1)                           = {'TS'}; % TOO SUPERFICIAL
                areaInfo(iarea).layerAssignment(recDepth >   .25 & recDepth <= 1)       = {'S'};  % SUPRAGRANULAR
                areaInfo(iarea).layerAssignment(recDepth >= -.25 & recDepth <= 0.25)    = {'G'};  % GRANULAR
                areaInfo(iarea).layerAssignment(recDepth >= -.75 & recDepth < -.25)     = {'I'};  % INFRAGRANULAR
                areaInfo(iarea).layerAssignment(recDepth <-.75)                         = {'TD'};  % TOO DEEP
            case 'V4'
                areaInfo(iarea).layerAssignment(recDepth > 1)                           = {'TS'}; % TOO SUPERFICIAL
                areaInfo(iarea).layerAssignment(recDepth >   .1  & recDepth <= 1)       = {'S'};  % SUPRAGRANULAR
                areaInfo(iarea).layerAssignment(recDepth >= -.1  & recDepth <= 0.1)     = {'G'};  % GRANULAR
                areaInfo(iarea).layerAssignment(recDepth >= -.75 & recDepth < -.1)      = {'I'};  % INFRAGRANULAR
                areaInfo(iarea).layerAssignment(recDepth <-.75)                         = {'TD'};  % TOO DEEP
        end
        
        % set channels to exclude based on recording depth
        if any( strcmpi(cfg.channelSelectionCriteria, 'depth') )
            
            cond1 = ~strcmpi(areaInfo(iarea).layerAssignment,'TS');
            cond2 = ~strcmpi(areaInfo(iarea).layerAssignment,'TD');
            
            areaInfo(iarea).ChanIdx.include = ...
            areaInfo(iarea).ChanIdx.include & (cond1(:)' & cond2(:)');

        end
    else
        warning('Cannot find selected coordinates, run summaryLaminarPlot.m and determineLayer4.m')
        warningGiven = true;
        areaInfo(iarea).Coordinates = NaN(nchan,1);
    end  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% get SNR (Signal to noise ratio)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loads file SNR.mat 
    
    loadfilename = [read_path 'SNR.mat'];
    if exist(loadfilename,'file')
        load(loadfilename,'SNR')
        areaInfo(iarea).SNR = SNR(string(SNR.area) == areaInfo(iarea).name, :);
    else
        warning([areaInfo(iarea).name ': cannot find channels to exclude based on SNR, run NCS_SNR.m'])
        warningGiven = true;
        areaInfo(iarea).SNR.SNR = zeros(1, nchan);
    end
    
    if any( strcmpi(cfg.channelSelectionCriteria, 'snr') )
        
        cond1 = areaInfo(iarea).SNR.SNR >= cfg.snr;
        
        areaInfo(iarea).ChanIdx.include = ...
            areaInfo(iarea).ChanIdx.include & cond1(:)' ;
    end    
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% set indices for bipolar channels
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tmp = flip(conv(double(flip(~areaInfo(iarea).ChanIdx.include)), [1 1])); % add one channel index on each side of an excluded channel to exclude from bipolar analyses
    tmp(1) = [];
    tmp(nchan:end) = [];
    areaInfo(iarea).ChanIdx.include_bip = ~(tmp>0);
    
end
