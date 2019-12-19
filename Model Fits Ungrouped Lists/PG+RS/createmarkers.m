function markers = createmarkers(parms) 
% CREATEMARKERS Creates a set of position markers
%..........................................................................

markers = zeros(parms.ll,parms.ll);  
for i=1:parms.ll
    for j=1:parms.ll
        markers(i,j) = parms.ItemDistinct.^abs(i-j);
    end
end

% Weighting and normalisation
%--------------------------------------------------------------------------
markers = parms.ItemWeight*(markers./repmat(sum(markers'),parms.ll,1)');