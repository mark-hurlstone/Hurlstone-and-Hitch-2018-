
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                            %
% COMPETITIVE QUEUING MODEL OF SERIAL RECALL %
% BASED ON FARRELL & LEWANDOWSKY (2004)      %
%                                            %                          
% THIS IS USED TO FIT THE PM+OI MODEL        %
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
fminparms.OutInt       = .5;   % Weighting of output interference

% SCALING OF RT
fminparms.Scaling      = 50;   % Iteration to ms scaling parameter


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% FMINSEARCH FOR PM+OI MODEL
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

% FOUR-ITEMS
%--------------------------------------------------------------------------
%       accspc: [0.7548 0.5764 0.5408 0.6240]
%       crtspc: [1.5120e+03 1.6585e+03 1.6003e+03 1.4323e+03]
%        trans: [0.0082 0.0316 0.1300 0.6240 0.1436 0.0491 0.0135]
%      transrt: [5.6042e+03 3.4077e+03 2.8620e+03 2.4488e+03 1.8612e+03 2.1781e+03 2.3207e+03]
%   finalstate: [-2.4292e+04 4.8593e+04 4.8604e+04 0.5742 0.6611 0.4354 71.0798]
% fltrdtransrt: [1.1761e+03 639.3518 278.5316 -97.3490 219.3428 574.5770 788.2510]

% FIVE-ITEMS
% %--------------------------------------------------------------------------
%       accspc: [0.6496 0.4910 0.4294 0.4288 0.5316]
%       crtspc: [1.9594e+03 2.0999e+03 2.0959e+03 1.9852e+03 1.7876e+03]
%        trans: [0.0093 0.0276 0.0576 0.1406 0.5061 0.1475 0.0673 0.0315 0.0125]
%      transrt: [6.5126e+03 4.6894e+03 3.6334e+03 3.1840e+03 2.8878e+03 2.4058e+03 2.5640e+03 2.6933e+03 2.8316e+03]
%   finalstate: [-3.5125e+04 7.0259e+04 7.0272e+04 0.4804 0.6411 0.3856 82.6199]
% fltrdtransrt: [1.0853e+03 817.9073 514.9429 191.8531 -203.4084 233.5581 421.8851 638.6968 841.5277]

% SIX-ITEMS
% %--------------------------------------------------------------------------
%       accspc: [0.6106 0.4186 0.3518 0.3370 0.3280 0.3940]
%       crtspc: [1.9785e+03 2.1194e+03 2.1249e+03 1.9909e+03 1.8612e+03 1.6404e+03]
%        trans: [0.0081 0.0210 0.0382 0.0659 0.1397 0.4067 0.1474 0.0821 0.0489 0.0288 0.0130]
%      transrt: [1x11 double]
%   finalstate: [-4.6001e+04 9.2009e+04 9.2023e+04 0.4410 0.6122 0.5334 82.1571]
% fltrdtransrt: [1.3216e+03 721.0820 663.9834 383.7089 71.4117 -251.9306 31.2153 290.5450 361.1807 342.6167 632.8365]