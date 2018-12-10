function checkconservation(G,nodnod,J,P,boundaryInd)

nnod = numnodes(G);

% conservation of flow
for inod = 1:nnod
    nodes = nodnod(inod,:);
    nodes(nodes==0) = [];
    if ismember(inod, boundaryInd)
        continue
    end
    qSum = 0;
    for j = 1:length(nodes)
        jedge = findedge(G,inod,nodes(j));
        qtemp = (P(inod) - P(nodes(j)))*J(jedge);
        qSum = qSum + qtemp;
    end
    if qSum > 0.01
        error('Flow Convervation Violation on node: %i\n',inod)
    end
end



