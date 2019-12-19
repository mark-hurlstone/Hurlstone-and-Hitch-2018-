function markers = combine(parms,fminparms)
% COMBINE Combines the primacy gradient and position markers
%..........................................................................

% Retrieve primacy gradient
%--------------------------------------------------------------------------
primacy = primgrad(parms);

% Retrieve position markers
%--------------------------------------------------------------------------
markers = createmarkers(parms,fminparms); 

% Weight primacy gradient and position markers
%--------------------------------------------------------------------------
markers = (1-parms.Mix).*markers + parms.Mix.*repmat(primacy,parms.ll,1);