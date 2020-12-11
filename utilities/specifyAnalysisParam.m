function [analysisParam] = specifyAnalysisParam(analysisType, recInfo)
% specify parameters of the analysis

%% check input

if nargin==2
    assert(height(recInfo)==1, 'height recInfo ~= 1')
    analysisParam = set_cfg(recInfo);
else
    analysisParam.fixAlignTimeWindow = [NaN];
end

%% -------------------------------------------------------------------- %%%
%%% settings general
%%% ------------------------------------------------------------------- %%%

analysisParam.ignoreDirection = true;% boolean, collapse across stimulus direction

analysisParam.binSize = 0.01; % in seconds

% set datatype
analysisParam.dataType = 'hash';
% analysisParam.dataType = 'unit';
% analysisParam.dataType = 'LFPb';
% analysisParam.dataType = 'microsaccades'; 

% analysisParam.dataType = 'hash_conv';% only for population analyses

analysisParam.convoluteSpikes = 0;%spikes convolved with Gaussian, extends dataType with '_conv' when saving
analysisParam.convoluteSpikesSigma = 0.010;%spikes convolved with Gaussian (seconds)

analysisParam.event2analyse = {'Fix';'Stim';'Cue';'Dim1';'Dim2'};
analysisParam.event2analyse = {'Cue'};
%
analysisParam.analysisType = analysisType;

analysisParam.matchRateAcrossCond = 0; % boolean, if true, equate firing rate across conditions (to match condition with lowest rate)
analysisParam.excludeMicrosaccades = 0; % boolean, excludes trials with microsaccades during relevant time window
analysisParam.excludeUnits = 1; % exclude units based on drift in activity 
%% -------------------------------------------------------------------- %%%
%%% settings analysis specific params
%%% ------------------------------------------------------------------- %%%

switch analysisType

    %-----------------------
    % Firing rate
    %-----------------------
    case 'rate'
%         analysisParam.windowType = 'sliding'; % type of the analysis window ('sliding', 'fullTimeEpoch', 'fixEvIgn')
        analysisParam.windowType = 'fullTimeEpoch'; % type of the analysis window ('sliding', 'fullTimeEpoch', 'fixEvIgn')
        if (strcmp(analysisParam.windowType,'sliding'))
            analysisParam.stepSize = 0.005; % sliding step for the sliding window (second)
            analysisParam.windowSize = 0.05; % size of the analysis window (second)
        end
        
        analysisParam.event = {...
            'Fix',  true, 'fixed', 'NLX', 'NLX_SUBJECT_START',  'NLX_STIM_ON',  [analysisParam.fixAlignTimeWindow(1) 0.5];
            'Stim', true, 'fixed', 'NLX', 'NLX_STIM_ON',        'NLX_CUE_ON',   [-0.3 0.03];
            'Cue',  false, 'full', 'NLX', 'NLX_CUE_ON',        'NLX_DIMMING1', [0.400 0.03];
            'Dim1', false, 'fixed', 'NLX', 'NLX_DIMMING1',      'NLX_DIMMING2', [-1.000 0.03];
            'Dim2', false, 'fixed', 'NLX', 'NLX_DIMMING2',      'NLX_DIMMING3', [-1.000 0.03]};
        analysisParam.event = cell2table(analysisParam.event, 'VariableNames', {'Epoch','ConcatenateCond','Extraction','EventType','Event','NextEvent','TimeWindow'});

        %-----------------------
        % Fano factor
        %-----------------------
    case 'FF'
        analysisParam.windowType = 'single'; % type of the analysis window
%         analysisParam.windowType = 'sliding'; % type of the analysis window
        analysisParam.windowType = 'varSize'; % type of the analysis window
%         analysisParam.windowType = 'fullTimeEpoch'; % type of the analysis window
        % 'single', 'sliding', 'varSize', 'fullTimeEpoch'
        switch analysisParam.windowType
            case 'sliding'
                analysisParam.stepSize = 0.01; % sliding step for the sliding window (second)
                
                analysisParam.windowSize = 0.1; % size of the analysis window (second)
            case 'fullTimeEpoch'
                analysisParam.windowSize = 0.01; % size of the analysis window (second)
        end
        
        analysisParam.event = {...
            'Fix',  true, 'fixed', 'NLX', 'NLX_SUBJECT_START',  'NLX_STIM_ON',  [analysisParam.fixAlignTimeWindow(1) 0.5];
            'Stim', true, 'fixed', 'NLX', 'NLX_STIM_ON',        'NLX_CUE_ON',   [-0.3 0.03];
            'Cue',  false, 'fixed', 'NLX', 'NLX_CUE_ON',        'NLX_DIMMING1', [0.400 1.0];
            'Dim1', false, 'fixed', 'NLX', 'NLX_DIMMING1',      'NLX_DIMMING2', [-1.000 0.03];
            'Dim2', false, 'fixed', 'NLX', 'NLX_DIMMING2',      'NLX_DIMMING3', [-1.000 0.03]};
        analysisParam.event = cell2table(analysisParam.event, 'VariableNames', {'Epoch','ConcatenateCond','Extraction','EventType','Event','NextEvent','TimeWindow'});
                
        %-----------------------
        % Hidden Markov Model, individual area
        %-----------------------
    case {'HMM'}

        analysisParam.event = {...
            'Fix',  true, 'full', 'NLX', 'NLX_SUBJECT_START',  'NLX_STIM_ON',  [analysisParam.fixAlignTimeWindow(1) 0.03];...
            'Stim', true, 'full', 'NLX', 'NLX_STIM_ON',        'NLX_CUE_ON',   [-0.3 0.03];
            'Cue',  false, 'full', 'NLX', 'NLX_CUE_ON',        'NLX_DIMMING1', [0.400 0.03];
            'Dim1', false, 'fixed', 'NLX', 'NLX_DIMMING1',      'NLX_DIMMING2', [-1.000 0.03];
            'Dim2', false, 'fixed', 'NLX', 'NLX_DIMMING2',      'NLX_DIMMING3', [-1.000 0.03]};
        analysisParam.event = cell2table(analysisParam.event, 'VariableNames', {'Epoch','ConcatenateCond','Extraction','EventType','Event','NextEvent','TimeWindow'});
        
        analysisParam.numFold = 2;  % number of cross-validtaion folds
        analysisParam.minNumState = 2; % Minimum number of states the model tries to fit.
        analysisParam.maxNumState = 2; % Maximum number of states the model tries to fit.
        
        analysisParam.binSize = 0.01; % bin size to bin data for analysis (second)
        
        %%% population analysis param
        analysisParam.inclBorder = true;

        
        %-----------------------
        % Hidden Markov Model 2-area
        %-----------------------
    case {'HMM_2area'}
        analysisParam.event = {...
            'Fix',  true, 'full', 'NLX', 'NLX_SUBJECT_START',  'NLX_STIM_ON',  [analysisParam.fixAlignTimeWindow(1) 0.03];...
            'Stim', true, 'full', 'NLX', 'NLX_STIM_ON',        'NLX_CUE_ON',   [-0.1 0.03];
            'Cue',  false, 'full', 'NLX', 'NLX_CUE_ON',        'NLX_DIMMING1', [0.400 0.03];
            'Dim1', false, 'fixed', 'NLX', 'NLX_DIMMING1',      'NLX_DIMMING2', [-1.000 0.03];
            'Dim2', false, 'fixed', 'NLX', 'NLX_DIMMING2',      'NLX_DIMMING3', [-1.000 0.03]};
        analysisParam.event = cell2table(analysisParam.event, 'VariableNames', {'Epoch','ConcatenateCond','Extraction','EventType','Event','NextEvent','TimeWindow'});
        
        analysisParam.numFold = 2;  % number of cross-validtaion folds
        analysisParam.minNumState = 4; % 4 states with: (V1-On, V4-On), (V1-On,V4-Off), (V1-Off,V4-Off), (V1-Off,V4-Off)
        analysisParam.maxNumState = 4;
        
        analysisParam.binSize = 0.01; % bin size to bin data for analysis (second)
        
        analysisParam.constrainFiringRates    = 1; % Only for HMM that fits V1 and V4 data. If 1, this sets the mean emmission rates constant across the same state in each area. E.g. for the 2 states (V1 On - V4 On) and (V1 On - V4 Off), V1 will have the same firing rates
        
        
        %%% population analysis param
        analysisParam.inclBorder = true;
        
        %
        %-----------------------
        % Hidden Markov Model Selection
        %-----------------------
    case 'HMMselect'
%         analysisParam.event2analyse = {'Fix';'Stim';'Cue';'Dim1';'Dim2'}; % overwrite event2analyse

        analysisParam.event = {...
            'Fix',  true, 'full', 'NLX', 'NLX_SUBJECT_START',  'NLX_STIM_ON',  [analysisParam.fixAlignTimeWindow(1) 0.03];...
            'Stim', true, 'full', 'NLX', 'NLX_STIM_ON',        'NLX_CUE_ON',   [-0.1 0.03];
            'Cue',  false, 'full', 'NLX', 'NLX_CUE_ON',        'NLX_DIMMING1', [0.400 0.03];
            'Dim1', false, 'fixed', 'NLX', 'NLX_DIMMING1',      'NLX_DIMMING2', [-1.000 0.03];
            'Dim2', false, 'fixed', 'NLX', 'NLX_DIMMING2',      'NLX_DIMMING3', [-1.000 0.03]};
        analysisParam.event = cell2table(analysisParam.event, 'VariableNames', {'Epoch','ConcatenateCond','Extraction','EventType','Event','NextEvent','TimeWindow'});
        
        analysisParam.numFold = 4;  % number of cross-validtaion folds
        analysisParam.minNumState = 1; % Minimum number of states the model tries to fit.
        analysisParam.maxNumState = 8; % Maximum number of states the model tries to fit.
        
        analysisParam.binSize = 0.01; % bin size to bin data for analysis (second)
        
        analysisParam.inclusionThreshold = -0.1;
        analysisParam.inclusionTimeWindow = 0.2;
        analysisParam.inclusionErrorSelection = 'looError';
%         analysisParam.inclusionErrorSelection = 'veError';
%         analysisParam.inclusionErrorSelection = 'cvError';
        
        %-----------------------
        % cross-correlation of HMM phases, Hidden Markov Model, across Area
        %-----------------------
    case {'HMM_crossCorrelation'}
        
        analysisParam.event = {...
            'Fix',  true, 'full', 'NLX', 'NLX_SUBJECT_START',  'NLX_STIM_ON',  [analysisParam.fixAlignTimeWindow(1) 0.03];...
            'Stim', true, 'full', 'NLX', 'NLX_STIM_ON',        'NLX_CUE_ON',   [-0.1 0.03];
            'Cue',  false, 'full', 'NLX', 'NLX_CUE_ON',        'NLX_DIMMING1', [0.400 0.03];
            'Dim1', false, 'full', 'NLX', 'NLX_DIMMING1',      'NLX_DIMMING2', [-1.000 0.03];
            'Dim2', false, 'fixed', 'NLX', 'NLX_DIMMING2',      'NLX_DIMMING3', [-1.000 0.03]};
        analysisParam.event = cell2table(analysisParam.event, 'VariableNames', {'Epoch','ConcatenateCond','Extraction','EventType','Event','NextEvent','TimeWindow'});
                                
        % HMM settings
        analysisParam.numFold = 2;  % number of cross-validtaion folds
        analysisParam.HMMtype = 'HMM'; analysisParam.numState = 2; % number of states
        analysisParam.binSize = 0.01; % bin size to bin data for analysis (second)
        analysisParam.HMMbinSize = 0.01; % bin size to bin data for analysis (second)

        % Cross-correlation time series
        analysisParam.performShuffle    = true; % does a trial shuffling of data and recomputes cross correlation
        analysisParam.maxTimeLag        = 0.3; % in seconds
        %analysisParam.scaleOpt         = 'unbiased';
        analysisParam.scaleOpt          = 'coeff';
        analysisParam.normaliseInput    = 'demean';
        
        % cross-correlation transition times
        analysisParam.CC_type = 'timeSeries'; %jitter
%         analysisParam.CC_type = 'jitter'; %jitter
%         analysisParam.CC_type = 'cpt'; analysisParam.normalisation = 'cpt'; %coincidences per trigger
        analysisParam.minEpochDurationWin = 0;% minimum time window (seconds)
        
        analysisParam.jitterBinSize = 0.3;  % window size for spike jittering (second)
        
        %-----------------------
        % cv-err of HMM, Hidden Markov Model, Individual Area
        %-----------------------
    case 'HMM_cvStateBased'
        
        analysisParam.event = {...
            'Fix',  true, 'full', 'NLX', 'NLX_SUBJECT_START',  'NLX_STIM_ON',  [analysisParam.fixAlignTimeWindow(1) 0.03];...
            'Stim', true, 'fixed', 'NLX', 'NLX_STIM_ON',        'NLX_CUE_ON',   [-0.3 0.03];
            'Cue',  false, 'full', 'NLX', 'NLX_CUE_ON',        'NLX_DIMMING1', [0.400 0.03];
            'Dim1', false, 'fixed', 'NLX', 'NLX_DIMMING1',      'NLX_DIMMING2', [-1.000 0.03];
            'Dim2', false, 'fixed', 'NLX', 'NLX_DIMMING2',      'NLX_DIMMING3', [-1.000 0.03]};
        analysisParam.event = cell2table(analysisParam.event, 'VariableNames', {'Epoch','ConcatenateCond','Extraction','EventType','Event','NextEvent','TimeWindow'});
        
        analysisParam.numFold = 2;  % number of cross-validtaion folds
        analysisParam.minNumState = 2; % Minimum number of states the model tries to fit.
        analysisParam.maxNumState = 2; % Maximum number of states the model tries to fit.
        
        analysisParam.binSize = 0.01; % bin size to bin data for analysis (second)
        
        
        %-----------------------
        % Microsaccade cross-correlation with HMM param, Hidden Markov Model, Individual Area
        %-----------------------
    case 'HMM_microsaccades'
        
%         analysisParam.dataType = 'microsaccades30'; % overwrite dataType
        
        analysisParam.event = {...
            'Cue',  false, 'fixed', 'NLX', 'NLX_CUE_ON',        'NLX_DIMMING1', [0.400 2.0];
            'Dim1', false, 'fixed', 'NLX', 'NLX_DIMMING1',      'NLX_DIMMING2', [-1.000 0]
            'Dim2', false, 'fixed', 'NLX', 'NLX_DIMMING2',      'NLX_DIMMING3', [-1.000 0]};
        analysisParam.event = cell2table(analysisParam.event, 'VariableNames', {'Epoch','ConcatenateCond','Extraction','EventType','Event','NextEvent','TimeWindow'});
        
        analysisParam.CC_type = 'jitter'; %jitter
        analysisParam.CC_type = 'cpt'; %coincidences per trigger
        
        % HMM settings
        analysisParam.numFold = 2;  % number of cross-validtaion folds
        analysisParam.numState = 2; % number of states
        analysisParam.binSize = 0.01; % bin size to bin data for analysis (second)

        analysisParam.jitterBinSize = 0.3;  % window size for spike jittering (second)
        analysisParam.maxTimeLag = 0.3; % maximal lag of the correlation functions (second)

        analysisParam.HMMtype = 'HMM'; analysisParam.numState = 2; % number of states
%         analysisParam.HMMtype = 'HMM_2area'; analysisParam.numState = 4; % number of states

        analysisParam.rateWindow = [0 0.2]; % window in which a transition after a microsaccade is included
        
        analysisParam.smoothData = [4]; % smooth data before population average of number of timebins

        %-----------------------
        % Pupil diameter correlation with HMM param, Hidden Markov Model, Individual and across Area
        %-----------------------
    case 'HMM_pupil'
        
        analysisParam.dataType = 'pupil'; % overwrite dataType
        analysisParam.pupilType = 'baseline'; %
        
        analysisParam.event = {...
            'Stim', false, 'fixed', 'CTX', 'STIM_ON',          'CUE_ON',   [-0.3 0];};
        analysisParam.event = cell2table(analysisParam.event, 'VariableNames', {'Epoch','ConcatenateCond','Extraction','EventType','Event','NextEvent','TimeWindow'});
        
        analysisParam.HMMEvent = 'Cue';% use a different event for the HMM data 
        
        % HMM settings
        analysisParam.numFold = 2;  % number of cross-validtaion folds
        analysisParam.binSize = 0.01; % bin size to bin data for analysis (second)
        
        analysisParam.HMMtype = 'HMM'; analysisParam.numState = 2; % number of states
        analysisParam.HMMtype = 'HMM_2area'; analysisParam.numState = 4; % number of states
        
        analysisParam.inclBorder    = true;
               
        %-----------------------
        % Reaction times across HMM phases, Hidden Markov Model, Individual and across Area
        %-----------------------
    case 'HMM_RT'
        
        analysisParam.event = {...
            'Cue',  false, 'full', 'NLX', 'NLX_CUE_ON',        'NLX_DIMMING1', [0.400 0.03];
            'Dim2', false, 'fixed', 'NLX', 'NLX_DIMMING2',      'NLX_DIMMING3', [-1.000 0.03]};
        analysisParam.event = cell2table(analysisParam.event, 'VariableNames', {'Epoch','ConcatenateCond','Extraction','EventType','Event','NextEvent','TimeWindow'});

        
        analysisParam.concatenateAcrossEvents = true;
        analysisParam.nTrialThreshold = 2;

        analysisParam.HMMtypes = {'HMM','HMM_2area'};
        analysisParam.numFold = 2;  % number of cross-validtaion folds

        analysisParam.RT2use = 'RT_EPP';
        analysisParam.RT2transform = 'zscore_log';
%         analysisParam.RT2transform = 'none';
        
        analysisParam.RTlim = 0.7;
        
        %-----------------------
        % Rate per state, Hidden Markov Model Individual Area
        %-----------------------
    case {'HMM_rate'}
        analysisParam.event = {...
            'Fix',  true, 'full', 'NLX', 'NLX_SUBJECT_START',  'NLX_STIM_ON',  [analysisParam.fixAlignTimeWindow(1) 0.03];...
            'Stim', true, 'fixed', 'NLX', 'NLX_STIM_ON',        'NLX_CUE_ON',   [-0.3 0.03];
            'Cue',  false, 'full', 'NLX', 'NLX_CUE_ON',        'NLX_DIMMING1', [0.400 0.03];
            'Dim1', false, 'fixed', 'NLX', 'NLX_DIMMING1',      'NLX_DIMMING2', [-1.000 0.03];
            'Dim2', false, 'fixed', 'NLX', 'NLX_DIMMING2',      'NLX_DIMMING3', [-1.000 0.03]};
        analysisParam.event = cell2table(analysisParam.event, 'VariableNames', {'Epoch','ConcatenateCond','Extraction','EventType','Event','NextEvent','TimeWindow'});
        
        
        analysisParam.HMMtype = 'HMM';
%         analysisParam.HMMtype = 'HMM_2area';
        
        switch analysisParam.HMMtype
            case 'HMM'
                analysisParam.numFold = 2;  % number of cross-validtaion folds
                analysisParam.numState = 2; % Number of states
            case 'HMM_2area'
                analysisParam.numFold = 2;  % number of cross-validtaion folds
                analysisParam.numState = 4; % Number of states
        end
        
        analysisParam.binSize = 0.01; % bin size to bin data for analysis (second)
        
        analysisParam.inclBorder    = true;
        analysisParam.minEpochDurationWin   = 0.25; % minimum epoch duration
        
        %-----------------------
        % Analysis of HMM param, Hidden Markov Model Individual Area
        %-----------------------
    case 'HMM_transitionTime'
        
        analysisParam.event = {...
            'Fix',  true, 'full', 'NLX', 'NLX_SUBJECT_START',  'NLX_STIM_ON',  [analysisParam.fixAlignTimeWindow(1) 0.03];...
            'Stim', true, 'full', 'NLX', 'NLX_STIM_ON',        'NLX_CUE_ON',   [-0.1 0.03];
            'Cue',  false, 'full', 'NLX', 'NLX_CUE_ON',        'NLX_DIMMING1', [0.400 0.03];
            'Dim1', false, 'fixed', 'NLX', 'NLX_DIMMING1',      'NLX_DIMMING2', [-1.000 0.03];
            'Dim2', false, 'fixed', 'NLX', 'NLX_DIMMING2',      'NLX_DIMMING3', [-1.000 0.03]};
        analysisParam.event = cell2table(analysisParam.event, 'VariableNames', {'Epoch','ConcatenateCond','Extraction','EventType','Event','NextEvent','TimeWindow'});
        
        
        analysisParam.HMMtype = 'HMM';
        %         analysisParam.HMMtype = 'HMM_2area';
        
        switch analysisParam.HMMtype
            case 'HMM'
                analysisParam.numFold = 2;  % number of cross-validtaion folds
                analysisParam.numState = 2; % Number of states
            case 'HMM_2area'
                analysisParam.numFold = 2;  % number of cross-validtaion folds
                analysisParam.numState = 4; % Number of states
        end
        
        analysisParam.HMMbinSize = 0.01;
        analysisParam.binSize = 0.1; % bin size to bin data for analysis (second)
        analysisParam.inclBorder    = true;
        
        %-----------------------
        % Analysis of neural data relative to HMM transitions, Individual Area
        %-----------------------
    case 'HMM_transitionTriggeredAverage'
        
        analysisParam.event = {...
            'Fix',  true, 'full', 'NLX', 'NLX_SUBJECT_START',  'NLX_STIM_ON',  [analysisParam.fixAlignTimeWindow(1) 0.03];...
            'Stim', true, 'full', 'NLX', 'NLX_STIM_ON',        'NLX_CUE_ON',   [-0.1 0.03];
            'Cue',  false, 'full', 'NLX', 'NLX_CUE_ON',        'NLX_DIMMING1', [0.400 0.03];
            'Dim1', false, 'fixed', 'NLX', 'NLX_DIMMING1',      'NLX_DIMMING2', [-1.000 0.03];
            'Dim2', false, 'fixed', 'NLX', 'NLX_DIMMING2',      'NLX_DIMMING3', [-1.000 0.03]};
        analysisParam.event = cell2table(analysisParam.event, 'VariableNames', {'Epoch','ConcatenateCond','Extraction','EventType','Event','NextEvent','TimeWindow'});

        % HMM settings
        analysisParam.HMMtype = 'HMM';
        switch analysisParam.HMMtype
            case 'HMM'
                analysisParam.numFold = 2;  % number of cross-validtaion folds
                analysisParam.numState = 2; % Number of states
            case 'HMM_2area'
                error('not implemented yet')
        end
        analysisParam.binSize = 0.01; % bin size to bin data for analysis (second)
        analysisParam.HMMbinSize = 0.01; % bin size to bin data for analysis (second)
        
        analysisParam.CC_type = 'cpt'; analysisParam.normalisation = 'cpt'; %coincidences per trigger

        analysisParam.minEpochDurationWin = 0.10;% minimum time window (seconds)
        
        analysisParam.demean = 1; % useful for LFP to compensate for DC offset
        
        analysisParam.maxTimeLag = 0.1; % maximal lag of the correlation functions (second)

        %-----------------------
        % Spectral decomposition of neural data during HMM phases, Individual and across Area
        %-----------------------
    case 'HMM_spectrogram'
        analysisParam.event = {...
            'Fix',  true, 'full', 'NLX', 'NLX_SUBJECT_START',  'NLX_STIM_ON',  [analysisParam.fixAlignTimeWindow(1) 0.03];...
            'Stim', true, 'full', 'NLX', 'NLX_STIM_ON',        'NLX_CUE_ON',   [-0.1 0.03];
            'Cue',  false, 'full', 'NLX', 'NLX_CUE_ON',        'NLX_DIMMING1', [0.400 0.03];
            'Dim1', false, 'fixed', 'NLX', 'NLX_DIMMING1',      'NLX_DIMMING2', [-1.000 0.03];
            'Dim2', false, 'fixed', 'NLX', 'NLX_DIMMING2',      'NLX_DIMMING3', [-1.000 0.03]};
        analysisParam.event = cell2table(analysisParam.event, 'VariableNames', {'Epoch','ConcatenateCond','Extraction','EventType','Event','NextEvent','TimeWindow'});

        % HMM settings
        analysisParam.binSize = 0.01; % bin size to bin data for analysis (second)

        % spectrogram settings
        analysisParam.minEpochDurationWin   = 0.25; % minimum epoch duration
        analysisParam.tapers        = [4 7];
        analysisParam.fpass         = [4 200];
        analysisParam.inclBorder    = 1;
        analysisParam.nfft          = 1024; %force nfft, if not defined it is calculated.
        analysisParam.HMMtype = 'HMM';
%         analysisParam.HMMtype = 'HMM_2area';
        
        switch analysisParam.HMMtype
            case 'HMM'
                analysisParam.numFold = 2;  % number of cross-validtaion folds
                analysisParam.numState = 2; % Number of states
            case 'HMM_2area'
                analysisParam.numFold = 2;  % number of cross-validtaion folds
                analysisParam.numState = 4; % Number of states
        end
        
        %-----------------------
        % Microsaccades
        %-----------------------        
    case 'microsaccade_analyse_parameters'
        
        analysisParam.event = {...
            'Stim', false, 'full', 'NLX', 'NLX_STIM_ON',        'NLX_CUE_ON',   [-0.1 0.03];
            'Cue',  false, 'full', 'NLX', 'NLX_CUE_ON',        'NLX_DIMMING1', [0.400 0.03];
            'Dim1', false, 'fixed', 'NLX', 'NLX_DIMMING1',      'NLX_DIMMING2', [-1.000 0.03];
            'Dim2', false, 'fixed', 'NLX', 'NLX_DIMMING2',      'NLX_DIMMING3', [-1.000 0.03]};
        analysisParam.event = cell2table(analysisParam.event, 'VariableNames', {'Epoch','ConcatenateCond','Extraction','EventType','Event','NextEvent','TimeWindow'});

        analysisParam.polarResolution = 9.5; % polar histogram resolution of (pi / analysisParam.polarResolution)

        
    otherwise
        error(['Unknown predictor type: ' analysisType]);
end

%-------------------------------
%   'full'      full time epoch between two events, [t_beg, t_end] will
%               be added to the beginning and end of this epoch
%   'fixed'     fixed time epoch relative to the event onset

idx = ismember(analysisParam.event.Epoch, analysisParam.event2analyse);
analysisParam.event = analysisParam.event(idx,:);

end
