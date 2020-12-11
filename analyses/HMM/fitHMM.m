function [estTrans, estEmis, estPi0, aic, bic] = fitHMM( emissionSeq, numState, binSize, numIter )
% fits the HMM model to the emission sequence data and returns for each trial the most
% likely sequence of hidden states and the state probabilities
%
% arguments:
% ---------
%   emissionSeq     -   cell array of length numTrial, each element
%                       contains an array of size (numChannel, numBin),
%                       where numChannel is the number of neurons and
%                       numBin is the number of time bins on that trial
%   numState        -   number of states for the HMM
%
% outputs:
% -------
%   estTrans        -   estimated transition matrix, size (numState, numState)
%   estEmis         -   estimated emission matrix, size (numState, numChannel)
%   estPi0          -   estimated initial probability distribution, size
%                       numChannel
%   aic             -   Akaike information criterion
%   bic             -   Bayesian information criterion

  
    numChannel = size(emissionSeq{1},1);
    
    meanRate = zeros(1, numChannel);
    for iChann = 1:numChannel
        meanRate(iChann) = sum( cellfun(@(x) sum( x(iChann,:)), emissionSeq ) )/sum( cellfun(@(x) size(x,2), emissionSeq ) );
    end;

    % initialization for transition matrix
    %if (numState>1)
    %    trans = 0.01*ones(numState,numState)./(numState-1) + ((0.99*numState-1)/(numState-1))*eye(numState);
    %else
    %    trans = 1;
    %end;

    allEstTrans = zeros(numIter, numState, numState);
    allEstEmis = zeros(numIter, numState, numChannel);
    allLogLiks = zeros(numIter, 1);
    allEstPi0 = zeros(numIter, numState);
    
    for iIter = 1:numIter

        % initialize pi0 from Dirichlet(1,1, ..., 1) distribution
        pi0 = gamrnd(1, 1, 1, numState);
        pi0 = pi0./sum(pi0);
        
        % initialize each row of the transition matrix using Dirichlet distribution
        if (numState>1)
            trans = gamrnd(ones(numState) + ((0.1/binSize-1)*(numState-1) -1)*eye(numState),1);
            trans = trans./repmat(sum(trans,2),1,numState);
        else
            trans = 1;
        end;
         
        % initialize emission matrix randomly
        emis = unifrnd(0,repmat(2*meanRate,numState,1));
        %emis = unifrnd(0,100, numState, numChannel);
        %emis = emis.*repmat(meanRate./mean(emis,1), numState,1 );
        emis(emis<=0) = 0.01;    % the default minimal guessed number of spikes
    
        % fit hmm
        [estTrans, estEmis, estPi0, logliks] = hmmTrainPoisson(emissionSeq,trans,emis, pi0, ...
            'tolerance', 1.0e-5, 'etol', 1.0e-3, 'trtol', 1.0e-3,'maxiterations', 1000); % ,'Verbose',true
        allEstTrans(iIter,:,:) = estTrans;
        allEstEmis(iIter,:,:) = estEmis;
        allEstPi0(iIter,:) = estPi0;
        allLogLiks(iIter) = logliks(end);
        
    end;
 
    [loglik,ind] = max(allLogLiks);
    estTrans = reshape( allEstTrans(ind,:,:), numState, numState );
    estEmis = reshape( allEstEmis(ind,:,:), numState, numChannel );
    estPi0 = reshape( allEstPi0(ind,:), 1, numState );
    [~,ind] = sort(sum(estEmis,2));
    estTrans = estTrans(ind,ind);
    estEmis = estEmis(ind,:);
    estPi0 = estPi0(ind);
    
    p = (numState-1)*numState + numState*numChannel; % number of parameters
    N = sum( cellfun(@length, emissionSeq )  ); % num data points
    aic = -2*loglik + 2*p;
    bic = -2*loglik + p*log(N);
    
     
end
