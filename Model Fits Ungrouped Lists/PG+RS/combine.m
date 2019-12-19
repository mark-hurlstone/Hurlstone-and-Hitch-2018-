function markers = combine(parms,fminparms)
% COMBINE Combines the primacy gradient and position markers
%..........................................................................

% Retrieve primacy gradient
%--------------------------------------------------------------------------
primacy = primgrad(parms,fminparms);

% Retrieve position markers
%--------------------------------------------------------------------------
markers = createmarkers(parms); 

% Weight primacy gradient and position markers
%--------------------------------------------------------------------------
markers = (1-parms.Mix).*markers + parms.Mix.*repmat(primacy,parms.ll,1);