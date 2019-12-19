function markers = createmarkers(parms) 
markers = zeros(parms.ll,parms.ll);  
for i=1:parms.ll
    for j=1:parms.ll
        markers(i,j) = parms.ItemDistinct.^abs(i-j);
    end
end
markers = parms.ItemWeight * (markers./repmat(sum(markers'),parms.ll,1)');