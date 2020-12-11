function [condNames, numCond, trIdx] = getConditions(recInfo, cfg, ConcatenateCond, trialdata)
% [condNames, nCond, trIdx] = getConditions(recInfo, cfg, ConcatenateCond, trialdata)
% 
% Get condition names and trial indices
% 
% Parameters
% ----------
% recInfo : table
%     table with single row (one recording)
% cfg : struct
%     struct with condition selection parameters, specific to paradigm. Set
%     in 'set_cfg'
% ConcatenateCond : boolean
%     boolean indicating whether to concatinate across conditions. e.g.
%     when aligned to fixation, there is no difference between stimulus
%     conditions
% trialdata : struct
%     struct (optional) with trial info, obtained through 'prepareData'. If
%     set, then trIdx can be returned.
%
% Returns
% -------
% condNames : cell
%     cell array of size (1 x numCond) with condition name strings 
% numCond : double
%     double indicating number of conditions
% trIdx : array
%     array of size (numTrial x 1) with condition index per trial. Only
%     returned if nargin==3
%

switch recInfo.Task
    
    case 'gratc'
        if nargin<3
            ConcatenateCond = false;
            [condNames,~,numCond] = conditionNames_gratc(recInfo, cfg, ConcatenateCond);
            trIdx=[];
            return
        elseif nargin<4
            [condNames,~,numCond] = conditionNames_gratc(recInfo, cfg, ConcatenateCond);
            trIdx=[];
            return
        else
            [condNames,condCodes,numCond] = conditionNames_gratc(recInfo, cfg, ConcatenateCond);
        end
        
        if ~ConcatenateCond
            trIdx = zeros(length(trialdata),1);
            for icond = 1:size(condCodes,1)
                tmp = sum( ([trialdata.cond_num])' == condCodes(icond,:),2)==1;
                trIdx(tmp) = icond;
            end
        else
            trIdx = true(length(trialdata),1);
        end
        
end

%%% extra functions
function [condNames,condCodes,nCond] = conditionNames_gratc(recInfo, cfg, ConcatenateCond)
if cfg.ignoreDirection
    condNames = {'AttRF','AttOut1','AttOut2'};
    condCodes = [1 4; 2 5; 3 6];
else
    switch recInfo.Subject{1}
        case 'Taylor'
            if (recInfo.Date >= datetime('2017-07-27')) && ~ConcatenateCond
                % for later recordings we used stationary stimuli, so
                % collapse across stimulus direction
                condNames = {'AttRF','AttOut1','AttOut2'};
                condCodes = [1 4; 2 5; 3 6];
            else
                condNames = {'AttRF_Dir1','AttOut1_Dir1','AttOut2_Dir1','AttRF_Dir2','AttOut1_Dir2','AttOut2_Dir2'};
                condCodes = (1:6)';
            end
        case {'Wyman', 'Jones'}
            condNames = {'AttRF_Dir1','AttOut1_Dir1','AttOut2_Dir1','AttRF_Dir2','AttOut1_Dir2','AttOut2_Dir2'};
            condCodes = (1:6)';
        otherwise
            error('Subject unknown')
    end
end

if ConcatenateCond
    condNames = {'allCond'};
end
nCond = length(condNames);
