function [lnL AIC BIC] = maximumLogLikelihood(parms,binEdges,obsrts,prdrts)
% MAXIMUMLOGLIKELIHOOD Computes the maximum log-likelihood 
%
% DESCRIPTION 
% This function is used to organise the observed and predicted response 
% times into the group distribution retrieved using "vincentAveraging.m" 
% and specified in the input argument "binEdges", which contains the 
% aggregated response time quantiles.
%..........................................................................


% Get RT data and predictions in a suitable format
%==========================================================================

nPars  = max(obsrts(:,1)); % Number of participants (for BIC calc below)
obsrts = transdata(obsrts(:,2:end),parms);
prdrts = transdata(prdrts,parms);

% Create quantiles for defining observed response time bins
%==========================================================================

logLs = zeros(1,parms.ll * parms.ll); 

for i = 1:parms.ll * parms.ll
    
    % Remove zeros and NaN in data and predictions - ensure column vectors
    %----------------------------------------------------------------------
    yObs         = sort(nonzeros(nonnans(obsrts(:,i))));
    yPrd         = sort(nonzeros(nonnans(prdrts(:,i)))); 
 
    % Get response frequencies for observed response times
    %----------------------------------------------------------------------
    histCountObs = histc(yObs,[-Inf,binEdges(:,i)',Inf]);
    histCountObs = histCountObs(1:end-1);
    frObs        = histCountObs;
    frObs        = frObs(:);
    
    % Get probability masses for predicted response times
    %----------------------------------------------------------------------
    histCountPrd = histc(yPrd,[-Inf,binEdges(:,i)',Inf]);
    histCountPrd = histCountPrd(1:end-1);
    pmPrd        = histCountPrd/parms.nTrials;
    pmPrd        = pmPrd(:); pmPrd(pmPrd<eps) = eps;
       
    % Compute log-likelihood
    %---------------------------------------------------------------------- 
    logLs(i) = sum(frObs .* log(pmPrd)); 
    
end

% Compute "joint" log-likelihood statistic, AIC, and BIC
%==========================================================================

lnL = sum(logLs);
N   = ceil(sum(sum(obsrts>0)/nPars')); % Average observations per subject
AIC = -2*lnL + 2 * length(parms.LB); 
BIC = -2*lnL + length(parms.LB) * log(N);
