function [bestx,bestFval,bestdummy,bestoutput] = mywrapperLoopfmin(parmarray,obsrts,binEdges,parms)
% MYWRAPPERLOOPFMIN Wrapper function for fminsearch
%..........................................................................

bestFval   = realmax;
bestx      = parmarray.*0;
bestdummy  = 0;
bestoutput = 0;

% Generate starting points for fminsearch
%--------------------------------------------------------------------------
for p1 = .5             % fminparms.ItemWeight
    for p2 = .65        % fminparms.ItemDistinct
        for p3 = .95    % fminparms.ResSupp
            for p4 = 50 % fminparms.Scaling            
                tic
                [x,fval,dummy,output] = fminsearchbnd(@bof,[p1 p2 p3 p4],...
                    parms.LB,parms.UB);
                toc
                if fval < bestFval
                    bestFval = fval;
                    bestx = x;
                    bestdummy = dummy;
                    bestoutput = output;
                end
            end
        end
    end
end

% Nested function inherits data from wrapper function
%--------------------------------------------------------------------------
    function lnL = bof(fminparms)
        prdrts = cq(fminparms);
        [lnL,~,~] = maximumLogLikelihood(parms,binEdges,obsrts,prdrts);
        lnL = -lnL; % make lnL negative for fitting      
    end
end