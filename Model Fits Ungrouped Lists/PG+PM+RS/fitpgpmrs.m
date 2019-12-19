
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                            %
% COMPETITIVE QUEUING MODEL OF SERIAL RECALL %
% BASED ON FARRELL & LEWANDOWSKY (2004)      %
%                                            %                          
% THIS IS USED TO FIT THE PG+PM+RS MODEL     %
%                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mark Hurlstone                  %
% School of Psychology            %
% University of Western Australia %
% mark.hurlstone@uwa.edu.au       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
global parms score
parms.seed = 111211;

% Import data for fitting
obsrts = dlmread('E2Six.txt');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

parms.ll               = 6;    % List length
parms.nTrials          = 5000; % N simulation trials

% CQ MODEL PARAMETERS
parms.ExciteWeight     = 1.1;  % Excitatory weight
parms.InhibitWeight    = -0.1; % Inhibitory weight
parms.CQThresh         = 1;    % Threshold for response
parms.MaxIters         = 200;  % Max iterations for each response
parms.NoiseMean        = 0;    % Mean noise
parms.NoiseSD          = .04;  % Std.Dev of noise

% POSITION MARKING 
fminparms.ItemWeight   = 1;    % Activation of target item 
fminparms.ItemDistinct = .65;  % Distinctiveness of position markers

% PRIMACY GRADIENT 
fminparms.GradStart    = .6;   % Start value for primacy gradient
fminparms.GradDecrease = .85;  % Decrease in primacy gradient
parms.Mix              = .5;   % Weighting of primacy gradient and position
                               % markers

% RESPONSE SUPPRESSION
fminparms.ResSupp      = .95;  % Extent of response suppression

% OUTPUT INTERFERENCE
parms.OutInt           = 0;   % Weighting of output interference

% SCALING OF RT
fminparms.Scaling      = 50;   % Iteration to ms scaling parameter


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% FMINSEARCH FOR PG+PM+RS MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Convert parameters for fminsearch
%--------------------------------------------------------------------------
temp = struct2cell(fminparms);
for i=1:length(temp)
	parmarray(i)=temp{i};
end

% Create parameter boundaries 
%--------------------------------------------------------------------------
parms.LB = zeros(1,length(parmarray)); 
parms.UB = ones(1,length(parmarray)); parms.UB(end) = 200;

% Data prep: Remove the effects of output position for each participant
%--------------------------------------------------------------------------
zeros = obsrts == 0; obsrts(zeros) = NaN; 
for i = 1:max(obsrts(:,1))    
    index   = obsrts(:,1) == i;
    nTrials = sum(index);
    obsrts(index,parms.ll+2:end) = obsrts(index,parms.ll+2:end) ...
        - repmat(nanmean(obsrts(index,parms.ll+2:end)),nTrials,1);
end

% Set up fminsearch
%--------------------------------------------------------------------------
x          = parmarray; 
nfuncevals = 500;
tolerance  = .1;
defopts    = optimset('fminsearch');
options    = optimset(defopts,'Display','iter','TolFun',tolerance, ...
    'MaxFunEvals',nfuncevals);

% Initialize data storage
%--------------------------------------------------------------------------
fits.accspc      = zeros(1,parms.ll);
fits.crtspc      = zeros(1,parms.ll);
fits.trans       = zeros(1,parms.ll*2-1);
fits.transrt     = zeros(1,parms.ll*2-1);
fits.finalstate  = zeros(1,length(x)+3);

% Fit the model 
%--------------------------------------------------------------------------
binEdges              = vincentAveraging(parms,obsrts);
[x,fval,dummy,output] = mywrapperLoopfmin(parmarray,obsrts,binEdges,parms);
    
% Generate predictions for final parameters
%--------------------------------------------------------------------------
prdrts            = cq(x);
[lnL,AIC,BIC]     = maximumLogLikelihood(parms,binEdges,obsrts,prdrts);
finalstate        = [lnL AIC BIC x];

% Store the predicitions
%--------------------------------------------------------------------------
fits.accspc       = score.accspc./parms.nTrials; 
fits.crtspc       = score.crtspc./score.accspc;
fits.trans        = score.trans./sum(score.trans);
fits.transrt      = sum(score.transrt)./score.transreps;
fits.fltrdtransrt = score.fltrdtransrt./score.transreps;
fits.finalstate   = finalstate;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% PLOT PREDICTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Accuracy SPC
%-------------
subplot(2,2,1)
plot(fits.accspc)  
title('Accuracy SPC')
xlabel('Serial Position')
ylabel('Proportion Correct')

% Latency SPC
%------------
subplot(2,2,2)
plot(fits.crtspc)
title('Latency SPC')
xlabel('Serial Position')
ylabel('Latency (Iterations)')

% Transposition Gradients
%------------------------
subplot(2,2,3)
plot(fits.trans)
title('Transposition Gradient')
xlabel('Transposition Displacement')
ylabel('Proportion Responses')

% Latency-Displacement Functions
%-------------------------------
subplot(2,2,4)
plot(fits.fltrdtransrt)
title('Displacement Latencies')
xlabel('Transposition Displacement')
ylabel('Latency (Iterations)')

% PREDICTIONS

% FOUR-ITEMS
%--------------------------------------------------------------------------
%       accspc: [0.7520 0.5876 0.5568 0.6780]
%1.5305e+03
%       crtspc: [4.6303e+03 1.5638e+03 1.4350e+03 1.0338e+03]
%        trans: [0.0060 0.0285 0.1503 0.6436 0.1227 0.0395 0.0094]
%      transrt: [5.8428e+03 3.4800e+03 2.7744e+03 2.2854e+03 1.3839e+03 1.2545e+03 1.0870e+03]
%   finalstate: [-1.7214e+04 3.4439e+04 3.4457e+04 0.4334 0.5593 0.6006 0.9127 0.9060 77.4960]
% fltrdtransrt: [1.0676e+03 570.4163 245.5888 -93.6969 6.9332 -9.2106 29.3912]

% FIVE-ITEMS
%--------------------------------------------------------------------------
%       accspc: [0.6284 0.4442 0.4080 0.4194 0.5426]
%1.8265e+03
%       crtspc: [5.0757e+03 1.8655e+03 1.8222e+03 1.6657e+03 1.2003e+03]
%        trans: [0.0070 0.0209 0.0586 0.1753 0.4885 0.1474 0.0644 0.0296 0.0083]
%      transrt: [6.5409e+03 4.1013e+03 3.2686e+03 2.8238e+03 2.4731e+03 1.6563e+03 1.5754e+03 1.4541e+03 1.2494e+03]
%   finalstate: [-2.5515e+04 5.1041e+04 5.1060e+04 0.4223 0.6143 0.5529 0.9333 0.9623 81.2296]
% fltrdtransrt: [1.2356e+03 675.0287 381.1611 135.8117 -128.8318 -36.3136 -8.1688 18.7680 10.7423]
    
% SIX-ITEMS
%--------------------------------------------------------------------------
%       accspc: [0.5784 0.4172 0.3580 0.3544 0.3700 0.5048]
%2.1021e+03
%       crtspc: [5.8100e+03 2.1464e+03 2.1223e+03 1.9591e+03 1.6744e+03 1.0064e+03]
%        trans: [0.0050 0.0140 0.0325 0.0676 0.1659 0.4305 0.1495 0.0782 0.0375 0.0158 0.0036]
%      transrt: [1x11 double]
%   finalstate: [-3.3484e+04 6.6979e+04 6.6999e+04 0.2305 0.4501 0.9055 0.9619 0.8666 92.6983]
% fltrdtransrt: [1.4012e+03 926.8494 594.8348 366.6463 120.2997 -167.3993 -57.2729 -17.9499 -18.4183 -36.6058 -29.7390]


