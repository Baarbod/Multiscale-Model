function edgenod = setflowdirection(edgenod,negFlow)

% reverse pressure nodes based on negative flow
for ineg = 1:length(negFlow)
    edgeneg = negFlow(ineg);
    pneg1 = edgenod(edgeneg,1);
    pneg2 = edgenod(edgeneg,2);
    edgenod(edgeneg,1) = pneg2;
    edgenod(edgeneg,2) = pneg1;
end
