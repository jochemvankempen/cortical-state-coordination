function [fit, warningMessage] = fitj_exponentialDecay(x, y, beta0)
% [fit, warningMessage] = fitj_exponentialDecay(x, y, beta0)
%
% Parameters
% ----------
% x : array
%     array with x values (e.g. time)
% y : array 
%     array with y values
% beta0 : array
%     array of size (1 x 2) with initial parameter estimates
% 
% Returns
% -------
% fit : struct 
%     struct with model fit, fields : 
% 
%     - coefficients
%     - residTimePDF
%     - Tau : decay 
%


%%% fit exponential decay function 
modelfun = @(b,x) b(1) * exp(-b(2)*x(:, 1));
if ~exist('beta0','var') || isempty(beta0)
    beta0 = [.4, 10];
end

options = statset('MaxFunEvals', 1000);

lastwarn('')
tbl = table(x(:), y(:));
mdl = fitnlm(tbl, modelfun, beta0, 'options', options);
[LASTMSG, LASTID] = lastwarn;

fit.coefficients = mdl.Coefficients{:, 'Estimate'};
% Create smoothed/regressed data using the model:
fit.residTimePDF = fit.coefficients(1) * exp(-fit.coefficients(2)*x(:));
fit.Tau = 1/fit.coefficients(2);