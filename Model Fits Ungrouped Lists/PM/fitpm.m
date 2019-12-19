
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                            %
% COMPETITIVE QUEUING MODEL OF SERIAL RECALL %
% BASED ON FARRELL & LEWANDOWSKY (2004)      %
%                                            %                          
% THIS IS USED TO FIT THE PM MODEL           %
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
fminparms.ItemWeight   = 1;    % Activation of target item 
fminparms.ItemDistinct = .65;  % Distinctiveness of position markers

% PRIMACY GRADIENT 
parms.GradStart        = .6;   % Start value for primacy gradient
parms.GradDecrease     = .85;  % Decrease in primacy gradient
parms.Mix              = 0;    % Weighting of primacy gradient and position
                               % markers

% RESPONSE SUPPRESSION
parms.ResSupp          = 0;    % Extent of response suppression

% OUTPUT INTERFERENCE
parms.OutInt           = 0;   % Weighting of output interference

% SCALING OF RT
fminparms.Scaling      = 50;   % Iteration to ms scaling parameter


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% FMINSEARCH FOR PM MODEL
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
%       accspc: [0.7184 0.5698 0.5796 0.7150]
%       crtspc: [1.6038e+03 1.7564e+03 1.7776e+03 1.6067e+03]
%        trans: [0.0118 0.0375 0.1260 0.6457 0.1288 0.0389 0.0114]
%      transrt: [5.5943e+03 3.8800e+03 3.1517e+03 2.5079e+03 2.1373e+03 2.4935e+03 2.5878e+03]
%   finalstate: [-2.4419e+04 4.8845e+04 4.8854e+04 0.4848 0.6484 70.9243]
% fltrdtransrt: [1.0677e+03 740.6915 296.0919 -96.5425 318.2714 714.7362 898.3500]

% FIVE-ITEMS
%--------------------------------------------------------------------------
%       accspc: [0.6282 0.4824 0.4456 0.4812 0.6326]
%       crtspc: [2.2058e+03 2.4259e+03 2.4701e+03 2.4027e+03 2.1852e+03]
%        trans: [0.0115 0.0299 0.0550 0.1382 0.5340 0.1370 0.0564 0.0280 0.0102]
%      transrt: [7.4555e+03 5.0850e+03 4.2926e+03 3.8650e+03 3.2424e+03 2.8926e+03 3.3363e+03 3.3178e+03 3.7565e+03]
%   finalstate: [-3.5356e+04 7.0718e+04 7.0728e+04 0.4728 0.6588 91.6214]
% fltrdtransrt: [1.4024e+03 1.0064e+03 668.8957 265.9090 -228.4165 275.4364 741.7773 782.4250 1.3612e+03]

% SIX-ITEMS
%--------------------------------------------------------------------------
%       accspc: [0.5640 0.4252 0.3780 0.3686 0.4306 0.5674]
%       crtspc: [2.0232e+03 2.2075e+03 2.2478e+03 2.2572e+03 2.1954e+03 2.0281e+03]
%        trans: [0.0101 0.0240 0.0415 0.0655 0.1318 0.4556 0.1328 0.0648 0.0416 0.0218 0.0106]
%      transrt: [1x11 double]
%   finalstate: [-4.6034e+04 9.2075e+04 9.2085e+04 0.3974 0.6157 80.6719]
% fltrdtransrt: [998.6443 693.4141 785.0163 515.8330 149.7888 -304.0034 227.6351 432.7833 700.7613 842.1230 1.1300e+03]