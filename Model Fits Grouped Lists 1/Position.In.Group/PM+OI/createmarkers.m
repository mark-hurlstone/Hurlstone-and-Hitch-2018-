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

% Create position-within-group markers
%--------------------------------------------------------------------------
pwgMarkers = zeros(parms.ll,parms.ll);  
r = 0;
for i=1:length(parms.Grouping)
    for p=1:parms.Grouping(i) 
        r=r+1; j=0;
        for k=1:length(parms.Grouping) 
            for l=1:parms.Grouping(k) 
                j=j+1;
                pwgMarkers(r,j) = fminparms.ItemWeight * ...
                    fminparms.ItemDistinct.^abs(p-l);
            end
        end
    end
end

% Weighting and normalisation
%--------------------------------------------------------------------------

% Weight position-within-sequence and position-within-group markers
markers = (1-parms.Sg).*pogMarkers+parms.Sg.*pwgMarkers;

% Normalize markers
markers = markers./repmat(sum(markers'),parms.ll,1)';
