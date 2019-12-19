
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                            %
% COMPETITIVE QUEUING MODEL OF SERIAL RECALL %
% BASED ON FARRELL & LEWANDOWSKY (2004)      %
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
global score
parms.seed = 12;


%%%%%%%%%%%%%%
% PARAMETERS %
%%%%%%%%%%%%%%

parms.ll            = 6;      % List length
parms.nTrials       = 10000;  % N simulation trials

% CQ MODEL PARAMETERS
parms.ExciteWeight  = 1.1;    % Excitatory weight
parms.InhibitWeight = -0.1;   % Inhibitory weight
parms.CQThresh      = 1;      % Threshold for response
parms.MaxIters      = 200;    % Max iterations for each response
parms.NoiseMean     = 0;      % Mean noise
parms.NoiseSD       = .04;    % Std.Dev of noise
parms.Iters2MS      = 50;     % Iteration-to-ms scaling 

% POSITION MARKING 
parms.ItemWeight    = 1;      % Activation of target item 
parms.ItemDistinct  = .65;    % Distinctiveness of position markers 

% PRIMACY GRADIENT 
parms.GradStart     = .65;    % Start value for primacy gradient 
parms.GradDecrease  = .90;    % Decrease in primacy gradient     
parms.Mix           = .5;     % Weighting of primacy gradient and position markers 

% RESPONSE SUPPRESSION
parms.ResSupp       = .95;    % Extent of response suppression

% OUTPUT INTERFERENCE
parms.OutInt        = 0;      % Amount of output interference


% Show predictions
cq(parms)
score.accspc./parms.nTrials
score.crtspc./score.accspc
score.trans./sum(score.trans)
score.transrt./score.transreps 
score.fltrdtransrt./score.transreps 

% Accuracy SPC
subplot(2,2,1)
plot(score.accspc./parms.nTrials)  
title('Accuracy SPC')
xlabel('Serial Position')
ylabel('Proportion Correct')

% Latency SPC
subplot(2,2,2)
plot(score.crtspc./score.accspc)
title('Latency SPC')
xlabel('Serial Position')
ylabel('Latency (Iterations)')

% Transposition Gradients
subplot(2,2,3)
plot(log2(score.trans./sum(score.trans)))
title('Transposition Gradient')
xlabel('Transposition Displacement')
ylabel('Proportion Responses')

% Latency-Displacement Function
subplot(2,2,4)
plot(score.fltrdtransrt./score.transreps)
title('Displacement Latencies')
xlabel('Transposition Displacement')
ylabel('Latency (Iterations)')