function Z = statj_gauss(X, meanX, stdX)
%
% Z = statj_gauss(X, meanX, stdX)
%
% transform values from normal distribution to its corresponding values on
% a standard normal distribution
%

if nargin<3
    stdX = std(X);
end
if nargin<2
    meanX = mean(X);
end

Z = (X - meanX)./stdX;