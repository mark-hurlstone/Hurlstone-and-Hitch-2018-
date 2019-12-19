
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
obsrts = dlmread('E3Grouped.txt');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

parms.ll               = 6;     % List length
parms.Grouping         = [3 3]; % Grouping pattern
parms.nTrials          = 5000;  % N simulation trials

% CQ MODEL PARAMETERS
parms.ExciteWeight     = 1.1;   % Excitatory weight
parms.InhibitWeight    = -0.1;  % Inhibitory weight
parms.CQThresh         = 1;     % Threshold for response
parms.MaxIters         = 200;   % Max iterations for each response
parms.NoiseMean        = 0;     % Mean noise
parms.NoiseSD          = .04;   % Std.Dev of noise

% POSITION MARKING 
fminparms.ItemWeight   = 1;     % Activation of target item 
fminparms.ItemDistinct = .65;   % Distinctiveness of position markers
parms.Sg               = .5;    % Weighting of two sets of markers

% PRIMACY GRADIENT 
fminparms.GradStart    = .6;    % Start value for primacy gradient
fminparms.GradDecrease = .85;   % Decrease in primacy gradient
parms.Mix              = .5;    % Weighting of primacy gradient and position
                                % markers

% RESPONSE SUPPRESSION
fminparms.ResSupp      = .95;   % Extent of response suppression

% OUTPUT INTERFERENCE
parms.OutInt           = 0;     % Weighting of output interference

% SCALING OF RT
fminparms.Scaling      = 50;    % Iteration to ms scaling parameter


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

%       accspc: [0.5584 0.4222 0.4394 0.4680 0.3800 0.4586]
%       crtspc: [2.6529e+03 2.6295e+03 2.4702e+03 2.3003e+03 2.0317e+03 1.2326e+03]
%        trans: [0.0045 0.0127 0.0540 0.0600 0.1445 0.4544 0.1225 0.0727 0.0535 0.0165 0.0047]
%      transrt: [1x11 double]
%   finalstate: [-2.1389e+04 4.2790e+04 4.2809e+04 0.0399 0.4546 0.6902 0.9547 0.9996 113.6254]
% fltrdtransrt: [1.9402e+03 1.3495e+03 429.5573 569.7584 185.9598 -229.9779 -20.0732 53.1192 -108.4093 -22.0978 56.1889]