function markers = combine(parms)

% Retrieve primacy gradient
primacy = primgrad(parms);

% Retrieve position markers
markers = createmarkers(parms); 

% Weighting of primacy gradient and item markers
markers = (1-parms.Mix).*markers + parms.Mix.*repmat(primacy,parms.ll,1);