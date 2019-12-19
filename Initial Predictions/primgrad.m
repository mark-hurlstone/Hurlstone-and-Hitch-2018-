function pgrad = primgrad(parms) 
pgrad = zeros(1,parms.ll);
for i=1:parms.ll	
	pgrad(i) = parms.GradStart * parms.GradDecrease.^(i-1);
end