function [ time, window, numSubWindow, shiftTime ] = spike_getFFtimeWin( newTimeEpoch, analysisParam )
% extract from the analysis parameters the window size and time bins for
% the Fano Factor computation
%
% arguments:
% ----------
%   newTimeEpoch    -   [t_begin t_end], time epoch to perform analysis
%   analysisParam   -   structure, that contains analysis parameters
%
% outputs:
% --------
%   time            -   time bins for the Fano Factor computation
%   window          -   window size bins for the Fano Factor computation
%   numSubWindow    -   array that contains number of subwindows for the FF
%                       computation, for each time bin or for each window size
%   shiftTime       -   number that should be added to the time varible to
%                       align time to the time bin ceneters


    switch analysisParam.windowType
        case 'single'
            time = newTimeEpoch(1,1);
            window = newTimeEpoch(1,2)-newTimeEpoch(1,1);
            numSubWindow = 1;
            shiftTime = 0.5*(newTimeEpoch(1,2)-newTimeEpoch(1,1));
        case 'varSize'
            time = newTimeEpoch(1,1);
            minWin = 0.005;
            win = newTimeEpoch(1,2)-newTimeEpoch(1,1);
            window = [];
            numSubWindow = [];
            while (win >= minWin)
                numSubWindow = [numSubWindow; floor( 100*(newTimeEpoch(1,2)-newTimeEpoch(1,1))/(100*win) ) ];
                window = [window; win];
                win = win/2;
            end;
            shiftTime = 0.5*(newTimeEpoch(1,2)-newTimeEpoch(1,1));
        case 'sliding'
            time = newTimeEpoch(1,1):analysisParam.stepSize:newTimeEpoch(1,2)-analysisParam.windowSize;
            shiftTime = 0.5*analysisParam.windowSize;
            numSubWindow = 1;
            window = analysisParam.windowSize;
        case 'fullTimeEpoch'
            time = newTimeEpoch(1,1);
            window = analysisParam.windowSize;
            shiftTime = 0.5*mean(newTimeEpoch(:,2)-newTimeEpoch(:,1),1);
            numSubWindow = floor(100*(newTimeEpoch(:,2)-newTimeEpoch(:,1))/(100*analysisParam.windowSize))';
        otherwise
         error(['Unknow window type ' analysisParam.windowType]);
    end;


end

