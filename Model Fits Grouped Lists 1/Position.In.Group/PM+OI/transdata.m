function indTrnsMtx = transdata(data,parms)
% TRANSDATA Calculates transposition matrix for each trial
%
% DESCRIPTION
% This function takes the data or model predictions and spits out an LDF
% for each experimental or simulation trial. This is necessary so that the
% RTs for each transposition displacement can be sorted into bins.
%..........................................................................

% Initialise data storage
%--------------------------------------------------------------------------
n            = length(data); 
responses    = data(:,1:parms.ll);
latencies    = data(:,parms.ll+1:parms.ll*2);
indTrnsMtx   = zeros(n,parms.ll*parms.ll);
displacement = -parms.ll+1:1:parms.ll-1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% REMOVE REPETITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

repetitions = zeros(n,parms.ll); 
for i = 1:parms.ll-1
    for j = i+1:parms.ll
        reptmp1 = responses(:,i) == responses(:,j);
        reptmp2 = responses(:,j) == responses(:,i);
        repetitions(:,j) = repetitions(:,j) + reptmp1;
        repetitions(:,i) = repetitions(:,i) + reptmp2;       
    end
end
repetitions = repetitions > 0;
responses(repetitions) = NaN; 
latencies(repetitions) = NaN;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCT TRANSPOSITION MATRIX FOR EACH TRIAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

counter = 0;
for i = 1:parms.ll
    for j = 1:parms.ll
        transTmp = responses(:,i) == j;
        transTmp = transTmp>0;
        indTrnsMtx(transTmp,counter+j) = latencies(transTmp,i);    
    end
    counter = counter + parms.ll;
end
