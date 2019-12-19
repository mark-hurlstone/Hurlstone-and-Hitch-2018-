
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                            %
% COMPETITIVE QUEUING MODEL OF SERIAL RECALL %
% BASED ON FARRELL & LEWANDOWSKY (2004)      %
%                                            %                          
% THIS IS USED TO FIT THE PM+RS MODEL        %
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
parms.GradStart        = .6;   % Start value for primacy gradient
parms.GradDecrease     = .85;  % Decrease in primacy gradient
parms.Mix              = 0;    % Weighting of primacy gradient and position
                               % markers

% RESPONSE SUPPRESSION
fminparms.ResSupp      = 0;    % Extent of response suppression

% OUTPUT INTERFERENCE
parms.OutInt           = 0;    % Weighting of output interference

% SCALING OF RT
fminparms.Scaling      = 50;   % Iteration to ms scaling parameter


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% FMINSEARCH FOR PM+RS MODEL
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
%       accspc: [0.7050 0.5484 0.5180 0.5770]
%       crtspc: [2.2989e+03 2.3385e+03 2.1569e+03 1.5900e+03]
%        trans: [0.0075 0.0403 0.1870 0.5871 0.1032 0.0485 0.0266]
%      transrt: [8.3493e+03 4.9023e+03 3.7742e+03 3.3625e+03 2.2529e+03 2.3038e+03 2.5157e+03]
%   finalstate: [-1.8576e+04 3.7160e+04 3.7171e+04 0.9142 0.8000 0.9069 107.7731]
% fltrdtransrt: [1.5178e+03 702.7017 213.5311 -214.2553 131.0257 370.6421 741.7149]

% FIVE-ITEMS
% %--------------------------------------------------------------------------
%           accspc: [0.5546 0.4172 0.3866 0.3992 0.4454]
%           crtspc: [1.7293e+03 1.7998e+03 1.7681e+03 1.6303e+03 1.4140e+03]
%            trans: [0.0095 0.0293 0.0737 0.2001 0.4406 0.1226 0.0660 0.0366 0.0216]
%          transrt: [5.6607e+03 3.8077e+03 3.0123e+03 2.5301e+03 2.3121e+03 1.7684e+03 1.8372e+03 1.8587e+03 1.8912e+03]
%       finalstate: [-2.8318e+04 5.6643e+04 5.6656e+04 0.8184 0.8172 1.0000 71.3107]
%     fltrdtransrt: [849.8828 595.9391 289.4743 45.4944 -189.2162 21.6174 195.3417 306.5683 404.5901]

% SIX-ITEMS
% %--------------------------------------------------------------------------
%           accspc: [0.5076 0.3766 0.3342 0.3284 0.3396 0.3758]
%           crtspc: [2.0557e+03 2.1636e+03 2.1948e+03 2.0625e+03 1.9835e+03 1.7368e+03]
%            trans: [0.0074 0.0188 0.0382 0.0872 0.2016 0.3770 0.1199 0.0641 0.0414 0.0276 0.0168]
%          transrt: [1x11 double]
%       finalstate: [-3.7707e+04 7.5421e+04 7.5435e+04 0.9768 0.8354 0.9828 83.7230]
%     fltrdtransrt: [1.0539e+03 928.6335 577.5130 312.2407 3.4821 -270.1891 -18.3457 171.8553 264.6780 416.6325 404.0496]