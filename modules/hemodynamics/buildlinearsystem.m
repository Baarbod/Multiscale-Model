function [A, b] = buildlinearsystem(G,J,nodnod,bcnodind,...
    inflow_ind,outflow_ind,inflow_val,outflow_val)

dbstop if error


nnod = numnodes(G);
A = zeros(nnod,nnod); % linear system matrix A
b = zeros(nnod,1); % linear system matrix b

% setup [A] matrix of coefficients
f = waitbar(0,'Setting up [A] matrix...');
for inod = 1:nnod
    nodnear = nonzeros(nodnod(inod,:)); % nodes connected to inod
    edgenear = findedge(G,inod,nodnear); % edges connected to inod
    A(inod,nodnear) = J(edgenear);
    A(inod,inod) = -sum(A(inod,nodnear));
    waitbar(inod/nnod,f)
end
close(f)

% setup [b] for boundary nodes
f = waitbar(0,'Setting up [b] for boundary nodes...');
for ibc = 1:length(bcnodind)
    bcnod = bcnodind(ibc);
    if ismember(bcnod,inflow_ind)
        bcval = inflow_val;
    elseif ismember(bcnod,outflow_ind)
        bcval = outflow_val;
    end
    col = find(A(:,bcnod)~=0);
    for i = 1:length(col)
        row = col(i);
        b(row) = b(row) - bcval*A(row,bcnod);
        A(row,bcnod) = 0;
    end
    waitbar(ibc/length(inflow_ind),f)
end
close(f)