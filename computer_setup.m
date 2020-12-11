function [paths] = computer_setup
% [paths] = computer_setup
%
% sets path environment for specific computer.
%
% Returns
% -------
% paths : struct
%     structure with fields
%       
%     - base: local drive
%     - server: server location (optional, used only when (readDataFromServer || saveDataToServer) == 1
%


global HPC

% set base path depending on which computer the scripts are run on
switch upper(getComputerName)
    case 'YOUR-COMPUTER-NAME'
        paths.base ='[your-local-drive]\Thiele_attention_gratc_V1_V4_laminar\';
        paths.server = '[your-remote-drive]\Thiele_attention_gratc_V1_V4_laminar\';
    otherwise
        error('Computer not found')
end

if HPC
    % folder locations on HPC
    paths.base = '[your-HPC-drive]/Thiele_attention_gratc_V1_V4_laminar/';
    paths.server = '[your-HPC-drive]/Thiele_attention_gratc_V1_V4_laminar/';
end

