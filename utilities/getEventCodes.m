function [evntOutput] = getEventCodes(evntInput,type)
%
% convert event code to name, or the other way around
%
% Example
% -------
% ::
%
%     [evntOutput] = getEventCodes(evntInput,type);
%     [8] = getEventCodes('NLX_STIM_ON','NLX');
%     ['NLX_STIM_ON'] = getEventCodes(8,'NLX');
%
% Parameters
% ----------
% evntInput : string/double
%     string with event name, or double with event code
% type : string
%     string with type of eventcode
%
%     * 'NLX' : Neuralynx 8-bit eventcodes
%     * 'CTX' : Cortex eventcodes
%
% Returns
% -------
% evntOutput : string/double
%     string with event name, or double with event code
%


if isnumeric(evntInput) %if evntInput is numeric, find evnt name (string)
    getCode     = 'name';
    inpIdx      = 1;
    outIdx      = 2;
else
    getCode     = 'code';
    inpIdx      = 2;
    outIdx      = 1;
end



switch type
    case 'NLX'
        Codes = {...
            2       , 'NLX_RECORD_START'; % start recording eye data (grcjdru1), at same time as ctx event START_EYE_DATA
            4       , 'NLX_SUBJECT_START'; % fixation onset
            8       , 'NLX_STIM_ON';
            9       , 'NLX_EVENT_1';
            10      , 'NLX_EVENT_2';
            11      , 'NLX_EVENT_3';
            12      , 'NLX_EVENT_4';
            13      , 'NLX_EVENT_5';
            14      , 'NLX_EVENT_6';
            15      , 'NLX_EVENT_7';
            17      , 'NLX_TESTDIMMED';
            18      , 'NLX_DISTDIMMED';    
            20      , 'NLX_CUE_ON';
            21      , 'NLX_CUE_OFF';
            22      , 'NLX_DIST1DIMMED';
            23      , 'NLX_DIST2DIMMED';
            24      , 'NLX_SACCADE_START';
            25      , 'NLX_DIMMING1';
            26      , 'NLX_DIMMING2';
            27      , 'NLX_DIMMING3';
            32      , 'NLX_SUBJECT_END'; % end recording eye data (grcjdru1), at same time as ctx event END_EYE_DATA
            104     , 'NLX_TARG_REL';
            150     , 'NLX_RF'; % unique event code for RF paradigm
            151     , 'NLX_SACC'; % unique event code for SACCD paradigm
            152     , 'NLX_PUPIL'; % unique event code for PUPIL paradigm
            153     , 'NLX_SHAPES'; % unique event code for SHAPES paradigm
            249     , 'NLX_TASKPARAM_START'
            248     , 'NLX_TASKPARAM_END'
            250     , 'NLX_STIMPARAM_END';
            251     , 'NLX_STIMPARAM_START';
            252     , 'NLX_TRIALPARAM_END'; % start of trial parameters
            253     , 'NLX_TRIALPARAM_START'; % end of trial parameters
            254     , 'NLX_TRIAL_END'; % end of trial, last event
            255     , 'NLX_TRIAL_START'; % start of trial, first event
            
            };
        
    case 'CTX'
        Codes = {...
            8      , 'FIXATION_OCCURS';
            100    , 'START_EYE_DATA'; % start recording eye pos
            101    , 'END_EYE_DATA'; % end recording eye pos
            300    , 'TRIALPARAM_START';
            301    , 'TRIALPARAM_END';
            302    , 'STIMPARAM_START';
            303    , 'STIMPARAM_END';
            304    , 'STIM_SWITCH';
            305    , 'REWARDPARAM_START';
            306    , 'REWARDPARAM_END';
            4000   , 'CUE_ON';
            4001   , 'STIM_ON';
            4002   , 'TEST_DIMMED';
            4003   , 'DIST_DIMMED';
            4004   , 'DIST1_DIMMED';
            4005   , 'DIST2_DIMMED';
            4006   , 'DIMMING1';
            4007   , 'DIMMING2';
            4008   , 'DIMMING3';
            4009   , 'MICRO_STIM';
            4010   , 'FIXSPOT_OFF';
            4011   , 'STIM_OFF';
            4012   , 'CUE_OFF';
            5000   , 'BAR_RELEASE_BEFORE_TEST';
            5001   , 'EARLY_BAR_RELEASE';
            5002   , 'BAR_RELEASE_ON_TEST';
            5003   , 'BAR_RELEASE_ON_DIST';
            5004   , 'LATE_BAR_RELEASE_ON_TEST';
            5005   , 'LATE_BAR_RELEASE_ON_DIST';
            5006   , 'SACCADE_ONSET';
            5007   , 'NO_BAR_RELEASE';
            5008   , 'CHOICE_TARGET';
            5009   , 'CHOICE_DIST';
            5010   , 'CHOICE_NONE';
            6000   , 'BROKE_FIXATION';
            10000  , 'PARAMBASE';
            };
end

switch getCode
    case 'name'
        evntOutput = Codes{[Codes{:,1}] == evntInput, 2};
    case 'code'
        evntOutput = Codes{strcmpi(Codes(:,2), evntInput), 1};
end

