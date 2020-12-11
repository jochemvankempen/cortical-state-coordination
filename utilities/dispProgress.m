function dispProgress(iter,skipline,totalIter)
%
% dispProgress(iter,skipline,totalIter)
%
% display progress in analysis as dots in command window
% 
% Parameters
% ----------
% iter : double
%     the current iteration of loop
% skipline : double
%     after how many iterations is a new line initiated
% totalIter: double
%     the total number of iterations
%
% 
% Jochem van Kempen

if ~mod(iter,skipline) && (iter == totalIter)
    fprintf('.%d\\%d\n',iter,totalIter)
elseif ~mod(iter,skipline) 
    fprintf('.%d\\%d\n',iter,totalIter)
elseif iter == totalIter
    fprintf('.\n')
else
    fprintf('.')
end

