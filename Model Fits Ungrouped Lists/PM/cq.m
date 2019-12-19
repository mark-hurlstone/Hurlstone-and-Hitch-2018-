function prd = cq(parmarray) 
% CQ Main function for competitive queuing model
%..........................................................................

global parms score

% Retrieve parameter values from parmarray
%--------------------------------------------------------------------------
fminparms.ItemWeight   = parmarray(1);
fminparms.ItemDistinct = parmarray(2);
fminparms.Scaling      = parmarray(3);

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

% 1. RUN THE MODEL 
%==========================================================================

% CHANGE THIS CODE IN THE OTHER MODELS AS WELL!!!
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
                markers(:,recall) = markers(:,recall)*(1-parms.ResSupp);
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

% 2. Prepare predictions
%==========================================================================

% 2.1 Add constant to RTs at first output position
%allRts(:,1) = allRts(:,1)+40;

% 2.2 Apply scaling to RTs
allRts = allRts * fminparms.Scaling;

% 2.3 Remove the effects of output position
allRtsF = allRts - repmat(mean(allRts),parms.nTrials,1);

% 2.4 Store responses and RTs in prd
prd = [allRes,allRtsF];

% 2.5 Score the data
scoring(parms,allRes,allRts);
