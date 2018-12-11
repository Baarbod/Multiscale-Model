function hd = computehematocritPries(G,hd,D,Q,edgenod,boundaryInd)

nnod = numnodes(G);

for inod = 1:nnod
    
    if ismember(inod,boundaryInd)
        continue
    end
    
    segInInd = find(edgenod(:,2) == inod);
    segOutInd = find(edgenod(:,1) == inod);
    
    if isempty(segInInd) || isempty(segOutInd)
        % boundary nodes useless
        continue
    end
    
    if length(segInInd) > 1
        % weighted average
        numer = 0;
        denom = 0;
        Qf = 0;
        for i = 1:length(segInInd)
            iIn = segInInd(i);
            numer = numer + hd(iIn)*Q(iIn);
            denom = denom + Q(iIn);
            Qf = Qf + Q(iIn);
            Df = mean(D(segInInd));
        end
        hf = numer/denom;
    else
        hf = hd(segInInd);
        Qf = Q(segInInd);
        Df = D(segInInd);
    end
    
    if Qf == 0
        continue
    end
    
    if length(segOutInd) == 2
        % phase separation
        alpha = segOutInd(1);
        beta = segOutInd(2);
        Qa = Q(alpha);
        Qb = Q(beta);
        Da = D(alpha);
        Db = D(beta);
        A = -6.96*log(Da/Db)/Df;
        B = 1 + 6.98*(1-hf)/Df;
        X0 = 0.4/Df;
        FQB = Qa/Qf; % fractional blood flow into alpha daughter
        if FQB <= X0
            FQE = 0;
        elseif (FQB > X0) && (FQB < 1-X0)
            x = (FQB - X0)/(1 - 2*X0);
            y = A + B*log(x/(1 - x));
            FQE = exp(y)/(1 + exp(y));
        elseif FQB >= 1-X0
            FQE = 1;
        end
        hd(alpha) = FQE*hf*Qf/Qa;
        hd(beta) = (1 - FQE)*hf*Qf/Qb;
        
    elseif length(segOutInd) > 2
        hd(segOutInd) = hf;
        
    elseif length(segOutInd) == 1
        % conservation law
        hd(segOutInd) = hf;
    end
    if isnan(hd(segOutInd))
        fprintf('Node: %i, segIn: %i, segOut: %i\n',inod,segInInd,segOutInd)
    end
end

