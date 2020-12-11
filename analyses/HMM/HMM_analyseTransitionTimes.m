function HMM_analyseTransitionTimes(states, timeBins, analysisParam, area, resultSave)
% HMM_analyseTransitionTimes(msu, analysisParam, dataType, resultSave)
%
% perform analyses on times of transitions, e.g. hazard rate
% 
% Parameters
% ----------
% states : cell 
%     cell array of size (numTrial, 1) which contains HMM latent timeseries
% timeBins : cell 
%     cell array of size (numTrial, 1) which contains corresponding time
%     bins
% analysisParam : struct
%     structure, fields contain analysis parameters
% area: string
%     string with name of recording area
% resultSave : struct
%     structure that contains information relevant for saving and plotting
%     the result 
%
% Returns
% -------
% HMM_transitionTimes.mat : file
%     file that contains all the variables below
% edges : array
%     array of size (1 x numEdges) with the edges for the histogram
% transitionPDF : array
%     array of size (numState x numEdges) with transition pdf
%
%
% **HMM_transitionTimes.mat : file**
%     file that includes all the variables listed under *Returns*
% 
% 

fprintf('\t\tComputing HMM transition time statistics\n')

binSize = analysisParam.binSize;


%%% Get transition probability across time
[ transition_times, transition_label, timeBins ] = HMM_extractUpDownTransitionTimes( states, timeBins, analysisParam.numState, analysisParam.HMMbinSize);

[~,maxtwin]     = max(cellfun(@length, timeBins(:, 1), 'UniformOutput', true),[],1);
maxTimeBins     = timeBins{maxtwin};

edges = maxTimeBins(1):binSize:maxTimeBins(end);
% if maxTimeBins(end)-edges(end) > 0
%     edges = [edges edges(end)+binSize];
% end
% edges(end) = [];%timeBins{1}(end);

for istate = 1:analysisParam.numState
   
    allTimes = [transition_times{:, istate}];
    
    [n] = histcounts(allTimes, edges, 'Normalization', 'pdf');
    
    if istate==1
        transition_pdf          = NaN(analysisParam.numState, length(edges)-1);
        transition_probability  = NaN(analysisParam.numState, length(edges)-1);
        transition_count        = NaN(analysisParam.numState, length(edges)-1);
        transition_countDensity = NaN(analysisParam.numState, length(edges)-1);
        transition_cumcount     = NaN(analysisParam.numState, length(edges)-1);
        transition_cdf          = NaN(analysisParam.numState, length(edges)-1);
    end
    
    [transition_pdf(istate,:)]           = histcounts(allTimes, edges, 'Normalization', 'pdf');
    [transition_probability(istate,:)]   = histcounts(allTimes, edges, 'Normalization', 'probability');
    [transition_count(istate,:)]         = histcounts(allTimes, edges, 'Normalization', 'count');
    [transition_countDensity(istate,:)]  = histcounts(allTimes, edges, 'Normalization', 'countdensity');
    [transition_cumcount(istate,:)]      = histcounts(allTimes, edges, 'Normalization', 'cumcount');
    [transition_cdf(istate,:)]           = histcounts(allTimes, edges, 'Normalization', 'cdf');
    
%     n = histc(allTimes, edges)
%     transitionPDF(istate,:) = n/sum(n);
end


% save data
dirName = resultSave.dataDirName ;
if (resultSave.data)
    if (~isfolder(dirName))
        mkdir(dirName);
    end
    fileName = [dirName 'HMM_transitionTimeStats.mat'];
    save(fileName, 'edges', 'transition_pdf', 'transition_probability', 'transition_count', 'transition_countDensity', 'transition_cumcount', 'transition_cdf', 'transition_label');
end







