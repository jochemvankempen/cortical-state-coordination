function [estTrans, estEmis, estPi0, aic, bic] = fitHMM_constrE( emissionSeq, chanPerArea, numState, binSize, numIter )
% fits the HMM model to the emission sequence data and returns for each trial the most
% likely sequence of hidden states and the state probabilities
%
% arguments:
% ---------
%   emissionSeq     -   cell array of length numTrial, each element
%                       contains an array of size (numChannel, numBin),
%                       where numChannel is the number of neurons and
%                       numBin is the number of time bins on that trial
%   chanPerArea     -   array of size (numChannel, 1) with a double indicating
%                       to which area the channel belongs
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
    
    constrE = zeros(4,length(chanPerArea));
    nChan(1) = length(find(chanPerArea==1));
    nChan(2) = length(find(chanPerArea==2));
    constrE([1 3],1:nChan(1)) = 1;
    constrE([2 4],1:nChan(1)) = 2;
    constrE(1:2,nChan(1)+1:length(chanPerArea)) = 1;
    constrE(3:4,nChan(1)+1:length(chanPerArea)) = 2;
    
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
        %emis = unifrnd(0,repmat(2*meanRate,numState,1));
        emis = unifrnd(0,repmat(2*meanRate,numState,1));
        %emis = unifrnd(0,100, numState, numChannel);
        %emis = emis.*repmat(meanRate./mean(emis,1), numState,1 );
        emis(emis<=0) = 0.01;    % the default minimal guessed number of spikes
    
        % fit hmm
        [estTrans,estEmis,estPi0,logliks] = hmmTrainPoisson_constrE(emissionSeq,trans,emis, pi0, constrE, ...
            'tolerance', 1.0e-5, 'etol', 1.0e-3, 'trtol', 1.0e-3,'maxiterations', 1000);
        
        allEstTrans(iIter,:,:) = estTrans;
        allEstEmis(iIter,:,:) = estEmis;
        allEstPi0(iIter,:) = estPi0;
        allLogLiks(iIter) = logliks(end);
        
    end;
 
    [loglik,ind] = max(allLogLiks);
    estTrans = reshape( allEstTrans(ind,:,:), numState, numState );
    estEmis = reshape( allEstEmis(ind,:,:), numState, numChannel );
    estPi0 = reshape( allEstPi0(ind,:), 1, numState );

    [~,ind] = sort(sum(estEmis,2), 'ascend');
       
    % ind(1), V1Off-V4Off, this will be state 1
    % ind(4), V1On-V4On, this will be state 4
    ind2sort = ind(2:3);
    
    [~,indV1] = max(sum(estEmis(ind2sort,chanPerArea==1),2));
    [~,indV4] = max(sum(estEmis(ind2sort,chanPerArea==2),2));
    
    if indV1 == indV4
        warning('Cannot sort V1-V4 activity\n')
    else
        ind(2) = ind2sort(indV1); % V1 On, will be state 2
        ind(3) = ind2sort(indV4); % V4 On, will be state 3
    end
        
    estTrans = estTrans(ind,ind);
    estEmis = estEmis(ind,:);
    estPi0 = estPi0(ind);
 
    p = (numState-1)*numState + numState*numChannel; % number of parameters
    N = sum( cellfun(@length, emissionSeq )  ); % num data points
    aic = -2*loglik + 2*p;
    bic = -2*loglik + p*log(N);
    
     
end
