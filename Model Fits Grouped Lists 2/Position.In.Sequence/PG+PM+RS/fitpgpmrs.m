
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
%       accspc: [0.5300 0.4066 0.3844 0.4042 0.3780 0.4418]
%       crtspc: [2.8189e+03 2.8265e+03 2.6740e+03 2.5484e+03 2.1848e+03 1.4006e+03]
%        trans: [0.0039 0.0113 0.0264 0.0745 0.1912 0.4242 0.1392 0.0665 0.0340 0.0191 0.0097]
%      transrt: [1x11 double]
%   finalstate: [-2.1654e+04 4.3321e+04 4.3339e+04 0.0307 0.6149 0.7454 0.9649 0.7478 121.5719]
% fltrdtransrt: [1.8053e+03 1.4905e+03 987.8048 467.4563 76.5203 -226.9971 -67.3141 6.8833 123.9058 91.2689 56.6022]
