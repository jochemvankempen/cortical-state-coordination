function logEPoiss = logPoissonEmmission( seq, e )
% for each state calculate the probability of firing pattern observed in each bin (given spike counts)
% assume Poisson firing distribution with the rate of each neuron specified by the emission
% matrix e
% the function computes log of the probability, and also it neglects the
% factorial 

    L = size(seq, 2);
    N = size(seq, 1);
    numStates = size(e, 1);
    
    checkN = size(e, 2);
    if checkN ~= N
        error('Input size mismatch');
    end
    
    logE = log(e);

   
    
    logEPoiss = zeros(numStates, L);
    for iBin = 1:L
        for iState = 1:numStates
            % sum over neurons
            logEPoiss(iState, iBin) = logEPoiss(iState, iBin) - sum(e(iState,:)) + logE(iState, :)*seq(:,iBin);
        end;
    end;
    
end

