function markers = createmarkers(parms,fminparms) 
% CREATEMARKERS Creates a set of position markers
%..........................................................................

% Create position-of-group markers
%--------------------------------------------------------------------------
pogMarkers = zeros(parms.ll,parms.ll);
r = 0;
for g=1:length(parms.Grouping)
    for p=1:parms.Grouping(g) 
        r=r+1; j=0;
        for l=1:length(parms.Grouping)  
            for i=1:parms.Grouping(l) 
                j=j+1;
                pogMarkers(r,j) = fminparms.ItemWeight * ...
                    fminparms.ItemDistinct.^abs(g-l);
            end
        end
    end
end

% Create position-within-sequence markers
%--------------------------------------------------------------------------
% Position-within-sequence markers
pwsMarkers = zeros(parms.ll,parms.ll);
for i=1:parms.ll
    for j=1:parms.ll
        pwsMarkers(i,j) = fminparms.ItemWeight * ...
            fminparms.ItemDistinct.^abs(i-j);
    end
end

% Weighting and normalisation
%--------------------------------------------------------------------------

% Weight position-within-sequence and position-within-group markers
markers = (1-parms.Sg).*pogMarkers+parms.Sg.*pwsMarkers;

% Normalize markers
markers = markers./repmat(sum(markers'),parms.ll,1)';
