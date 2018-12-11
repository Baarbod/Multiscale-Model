function h = computehematocritGould(G,h,edgenod,d,M,q)

nnod = numnodes(G);

for inod = 1:nnod
    
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
        qf = 0;
        for i = 1:length(segInInd)
            iIn = segInInd(i);
            numer = numer + h(iIn)*q(iIn);
            denom = denom + q(iIn);
            qf = qf + q(iIn);
            df = mean(d(segInInd));
        end
        hf = numer/denom;
    else
        hf = h(segInInd);
        qf = q(segInInd);
        df = d(segInInd);
    end
    
    theta = zeros(length(segOutInd),1);
    sumqtheta = 0;
    for iseg = 1:length(segOutInd)
        seg = segOutInd(iseg);
        theta(iseg) = (d(seg)/df)^(1/M);
        sumqtheta = sumqtheta + theta(iseg)*q(seg);
    end
    
    adjustedHematocrit = qf*hf/sumqtheta;
    
    for iseg = 1:length(segOutInd)
        seg = segOutInd(iseg);
        h(seg) = theta(iseg)*adjustedHematocrit;
    end
            
end

        
        
        
        
        
        