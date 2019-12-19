function scoring(parms,allRes,allRts) 

global score

% Remove omissions
omissions = allRes == 0;
allRes(omissions) = NaN;
allRts(omissions) = NaN;

% ACC SPC
acctmp = allRes == repmat(1:parms.ll, parms.nTrials,1);
score.accspc = sum(acctmp); 

% RT SPC
allrtstmp = allRts;
incorrect = acctmp == 0;
allrtstmp(incorrect) = 0; 
score.crtspc = sum(allrtstmp);
score.rtspc = nansum(allRts);

% Transposition Gradient
displacement = (-parms.ll+1:1:parms.ll-1); 
trans = zeros(parms.nTrials,parms.ll);
for i=1:parms.ll*2-1
    transtmp = (repmat(1:parms.ll,parms.nTrials,1) - allRes) == displacement(i); 
	trans(:,i) = sum(transtmp');
end
score.trans = sum(trans); 

% Identify repetitions
repetitions = zeros(parms.nTrials,parms.ll);
for row=1:parms.ll-1
    for col=row+1:parms.ll
        reptmp1 = allRes(:,row) == allRes(:,col);
        reptmp2 = allRes(:,col) == allRes(:,row);
        repetitions(:,col) = repetitions(:,col)+reptmp1;
        repetitions(:,row) = repetitions(:,row)+reptmp2;    
    end
end

% Exclude repetitions
repetitions = repetitions > 0;
allRes(repetitions) = NaN;
allRts(repetitions) = NaN; 

% Unfiltered transposition latencies
trans = zeros(parms.nTrials,parms.ll);
transrt = zeros(parms.nTrials,parms.ll);
for i=1:parms.ll*2-1
	allrtstmp = allRts;
    transtmp = (repmat(1:parms.ll,parms.nTrials,1) - allRes) == displacement(i); 
	trans(:,i) = sum(transtmp');
    transrttmp = transtmp == 0; 
    allrtstmp(transrttmp) = 0;
    transrt(:,i) = sum(allrtstmp');
end
score.transreps = sum(trans);
score.transrt = sum(transrt);

% Filtered transposition latencies 
trans = zeros(parms.nTrials,parms.ll);
allRts = allRts - repmat(nanmean(allRts),parms.nTrials,1);
for i=1:parms.ll*2-1
    allrtstmp = allRts;
    transtmp = (repmat(1:parms.ll,parms.nTrials,1) - allRes) == displacement(i); 
	trans(:,i) = sum(transtmp');
    transrttmp = transtmp == 0; 
    allrtstmp(transrttmp) = 0;
    transrt(:,i) = sum(allrtstmp');
end
score.fltrdtransrt = sum(transrt);