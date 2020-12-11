function [time, window ] = spike_getRateTime( newTimeEpoch, analysisParam )
% generate time bins for the rate computation
%
% arguments:
% ----------
%   newTimeEpoch    -   [t_begin t_end], time epoch to perform analysis
%   windowSize      -   size of sliding window
%   stepSize        -   step size of sliding window
%
% outputs:
% --------
%   time            -   time bins for the rate computation
%   window          -   window for counting spikes in each time bin and on
%                       each trial (numTrial or 1, numTimeBin, 2), on the last
%                       dimension [t_beg t_end]


    switch analysisParam.windowType
        case 'sliding'
            assert( size(newTimeEpoch,1)==1, 'newTimeEpoch cannot have multiple entries')
                
            time = newTimeEpoch(1,1)+0.5*analysisParam.windowSize:analysisParam.stepSize:newTimeEpoch(1,2)-0.5*analysisParam.windowSize;
            window(1,:,:) = [time'-0.5*analysisParam.windowSize, time'+0.5*analysisParam.windowSize];
        case 'fullTimeEpoch'
            time = 0.5*mean(newTimeEpoch(:,2)+newTimeEpoch(:,1),1);
            window(:,1,:) = newTimeEpoch;
        otherwise
            error(['Unknow window type ' analysisParam.windowType]);
    end;

end

