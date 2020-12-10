function [path_list] = getListPath(Subject)
%
% get the path of csv lists, based on location of repository. Exception is
% possible for individual subjects
%
% [local drive]\lists\csv
%
% Example
% -------
% ::
% 
%     [path_list] = getListPath
%     [path_list] = getListPath('T')
%
%
% Parameters
% ----------
% Subject: string
%     string with subject name
%
% Returns
% -------
% path_list: string
%     string with path name
%
% Notes
% -----
% - Created, Jochem van Kempen, 2020/03/16
% - Adjusted, now determines path based on location of Lists repository (https://gitlab.com/thiele-lab/lists)
% 

repoName = 'lists';

p = mfilename('fullpath');
repoIdx = strfind(p,repoName);

path_list = [p(1:repoIdx-1) repoName filesep 'csv' filesep];

if (nargin==1)
    switch Subject
        case {'J','S','T','W'}
        case {'xxx'}
            % set exceptions to general path
            path_list = [];
        otherwise
            error('List path not defined for %s', Subject)
    end
end