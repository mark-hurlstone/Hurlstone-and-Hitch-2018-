function prd = cq(parmarray) 
% CQ Main function for competitive queuing model
%..........................................................................

global parms score

% Retrieve parameter values from parmarray
%--------------------------------------------------------------------------
fminparms.GradStart    = parmarray(1);
fminparms.GradDecrease = parmarray(2);
fminparms.ResSupp      = parmarray(3);
fminparms.Scaling      = parmarray(4);

% Retrieve lateral inhibition weight matrix
%--------------------------------------------------------------------------
W = weightMatrix(parms);

% Initialize scoring structure
%--------------------------------------------------------------------------
score.accspc       = zeros(1,parms.ll); 
score.rtspc        = zeros(1,parms.ll);
score.crtspc       = zeros(1,parms.ll);
score.trans        = zeros(1,parms.ll*2-1);
score.transreps    = zeros(1,parms.ll*2-1);
score.transrt      = zeros(1,parms.ll*2-1);
score.fltrdtransrt = zeros(1,parms.ll*2-1);

% Responses + latencies
%--------------------------------------------------------------------------
allRes = zeros(parms.nTrials,parms.ll); % All responses
allRts = zeros(parms.nTrials,parms.ll); % All latencies

% Randomization stuff
%--------------------------------------------------------------------------
randn('state',parms.seed);
rand('state',parms.seed);

% RUN THE MODEL 
%==========================================================================

gmarkers = combine(parms,fminparms);
for trial=1:parms.nTrials	
    markers = gmarkers;
    for pos=1:parms.ll 
        Vin = markers(pos,:)+(parms.NoiseSD.*randn(1,parms.ll)*parms.OutInt*pos);
        negVals = Vin < 0;
        Vin(negVals) = 0;
        Vout = W*Vin'+(parms.NoiseMean+parms.NoiseSD.*randn(1,parms.ll))';
        for cycle=1:parms.MaxIters
            if (max(Vout)>parms.CQThresh)
                [~, recall] = max(Vout);
                markers(:,recall) = markers(:,recall)*(1-fminparms.ResSupp);
                allRes(trial,pos) = recall;
                allRts(trial,pos) = cycle;
            break
            else
                Vin = Vout; 
                negVals = Vin<0;
                Vin(negVals) = 0; 
                Vout = W*Vin+(parms.NoiseMean+parms.NoiseSD.*randn(1,parms.ll))';
            continue 
            end
        end
    end
end

% PREPARE PREDICTIONS
%==========================================================================

% Add constant to RTs at first output position
%--------------------------------------------------------------------------
%allRts(:,1) = allRts(:,1)+40;

% Apply scaling to RTs
%--------------------------------------------------------------------------
allRts = allRts * fminparms.Scaling;

% Remove the effects of output position
%--------------------------------------------------------------------------
allRtsF = allRts - repmat(mean(allRts),parms.nTrials,1);

% Store responses and RTs in prd
%--------------------------------------------------------------------------
prd = [allRes,allRtsF];

% Score the data
%--------------------------------------------------------------------------
scoring(parms,allRes,allRts);
