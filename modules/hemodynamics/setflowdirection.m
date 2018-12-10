function edgenod = setflowdirection(edgenod,negFlow)

% reverse pressure nodes based on negative flow
f = waitbar(0,"Reversing pressure nodes for negative flows...");
for ineg = 1:length(negFlow)
    edgeneg = negFlow(ineg);
    pneg1 = edgenod(edgeneg,1);
    pneg2 = edgenod(edgeneg,2);
    edgenod(edgeneg,1) = pneg2;
    edgenod(edgeneg,2) = pneg1;
    waitbar(ineg/length(negFlow),f)
end
close(f)
