
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                            %
% COMPETITIVE QUEUING MODEL OF SERIAL RECALL %
% BASED ON FARRELL & LEWANDOWSKY (2004)      %
%                                            %                          
% THIS IS USED TO FIT THE PG+RS MODEL        %
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
obsrts = dlmread('E2Four.txt');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

parms.ll               = 4;    % List length
parms.nTrials          = 5000; % N simulation trials

% CQ MODEL PARAMETERS
parms.ExciteWeight     = 1.1;  % Excitatory weight
parms.InhibitWeight    = -0.1; % Inhibitory weight
parms.CQThresh         = 1;    % Threshold for response
parms.MaxIters         = 200;  % Max iterations for each response
parms.NoiseMean        = 0;    % Mean noise
parms.NoiseSD          = .04;  % Std.Dev of noise

% POSITION MARKING 
parms.ItemWeight       = 1;    % Activation of target item 
parms.ItemDistinct     = .65;  % Distinctiveness of position markers

% PRIMACY GRADIENT 
fminparms.GradStart    = .6;   % Start value for primacy gradient
fminparms.GradDecrease = .85;  % Decrease in primacy gradient
parms.Mix              = 1;    % Weighting of primacy gradient and position
                               % markers

% RESPONSE SUPPRESSION
fminparms.ResSupp      = .95;  % Extent of response suppression

% OUTPUT INTERFERENCE
parms.OutInt           = 0;    % Weighting of output interference

% SCALING OF RT
fminparms.Scaling      = 50;   % Iteration to ms scaling parameter


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% FMINSEARCH FOR PG+RS MODEL
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
%       accspc: [0.6744 0.4702 0.4700 0.6656]
%       crtspc: [1.7488e+03 1.7926e+03 1.6345e+03 1.2740e+03]
%        trans: [0.0057 0.0377 0.1607 0.5700 0.1928 0.0316 0.0014]
%      transrt: [6.5072e+03 3.9629e+03 3.2577e+03 2.6002e+03 1.3582e+03 1.1382e+03 833.4448]
%   finalstate: [-1.7735e+04 3.5479e+04 3.5490e+04 0.4311 0.8738 0.9742 85.4815]
% fltrdtransrt: [1.1660e+03 474.6091 254.6599 -36.2029 -184.2110 -289.7088 -388.2200]

% FIVE-ITEMS
%--------------------------------------------------------------------------
%       accspc: [0.5282 0.3496 0.3142 0.3352 0.5492]
%       crtspc: [1.7601e+03 1.7647e+03 1.6915e+03 1.4499e+03 980.3474]
%        trans: [0.0068 0.0222 0.0676 0.1780 0.4153 0.2205 0.0737 0.0150 0.0010]
%      transrt: [5.9526e+03 4.2389e+03 3.2006e+03 2.7175e+03 2.2653e+03 1.3360e+03 1.1569e+03 966.8921 715.1221]
%   finalstate: [-2.6217e+04 5.2443e+04 5.2455e+04 0.4895 0.9298 0.8738 77.0606]
% fltrdtransrt: [902.1036 547.5929 296.4278 127.8058 -42.7911 -115.0477 -196.0864 -218.8930 -219.4823]

% SIX-ITEMS
%--------------------------------------------------------------------------
%         accspc: [0.4330 0.2840 0.2370 0.2390 0.2850 0.4702]
%       crtspc: [2.8438e+03 2.8637e+03 2.8381e+03 2.5432e+03 2.1469e+03 1.3381e+03]
%        trans: [0.0064 0.0175 0.0390 0.0810 0.1681 0.3247 0.2136 0.1051 0.0348 0.0087 0.0011]
%      transrt: [4.2118e+03 4.1447e+03 3.5807e+03 3.1474e+03 2.8090e+03 2.3239e+03 2.1460e+03 1.9769e+03 1.7626e+03 1.4859e+03 1.1171e+03]
%   finalstate: [-3.4552e+04 6.9112e+04 6.9126e+04 0.5389 0.9543 0.8197 114.5741]
% fltrdtransrt: [1.0340e+03 1.0484e+03 587.0059 275.2593 110.3823 -68.0810 -161.8856 -195.7037 -234.4962 -256.9135 -163.5659]