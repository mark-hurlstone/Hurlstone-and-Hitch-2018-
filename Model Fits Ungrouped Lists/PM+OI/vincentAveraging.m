function binEdges = vincentAveraging(parms,obsrts)
% VINCENT AVERAGING Computes the vincent averaged RT quantiles
%
% DESCRIPTION 
% This function is used to organise the observed response times into bins. 
% The edges of these bins are based on response time quantiles corresponding 
% to the cumulative probabilities defined in the variable cumprob. The method
% employed here uses "vincent avergaging" which produces a group distribution,
% by arranging the data from each participant in ascending order, and
% computing quantiles for each participant. The quantiles are then averaged
% across participants to obtain a group distribution. Vincent avergaging
% preserves the shape of individual distributions of response times
% (Ratcliff, 1979).
%
% REFERENCES
% Ratcliff, R. (1978). A theory of memory retrieval. Psychological 
% Review,85,59-108.
%..........................................................................

nPars = max(obsrts(:,1));

% Create quantiles for defining observed response time bins
%==========================================================================

cumprob = .5;
binEdges = zeros(length(cumprob),parms.ll * parms.ll); 

for i = 1:nPars
    Obs = transdata(obsrts(obsrts(:,1) == i,2:end),parms); 
    for j = 1:parms.ll * parms.ll
        
        % Get rid of zeros and NaN in data and ensure column  vector
        %------------------------------------------------------------------
        yObs    = sort(nonzeros(nonnans(Obs(:,j)))); 
        
        % Get response time bin edges of observed response times
        %------------------------------------------------------------------
        if(isnan(quantile(yObs,cumprob)') == 0); 
            binEdges(:,j) = binEdges(:,j) + quantile(yObs,cumprob)';
        end
        
    end
end
% Vincent averaged quantiles
binEdges = binEdges./nPars;
