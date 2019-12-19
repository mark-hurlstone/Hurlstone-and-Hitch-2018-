function W = weightMatrix(parms)
W = zeros(parms.ll,parms.ll);
for i=1:parms.ll
	for j=1:parms.ll
		if (i==j)
			W(i,j) = parms.ExciteWeight;
		else
			W(i,j) = parms.InhibitWeight;
		end
	end
end
