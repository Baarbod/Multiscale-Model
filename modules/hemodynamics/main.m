clear, clc, close all

%% Load network
[file, path] = uigetfile('C:\Users\Baarbod\Desktop\multiscale_sandbox\data\processed');
network_path = [path file];
load(network_path)

xyz = o.xyz;
Adj = o.Adj;
d = o.d;
lseg = o.lseg;
edgenod = o.edgenod;
nodnod = o.nodnod;
% inflowInd = o.inflowInd;
% outflowInd = o.outflowInd;
G = graph(Adj);

% network au_sect_150_125_200 only
deg1 = find(degree(G)==1);
deg1zcoord = xyz(degree(G)==1,3);
inIndex = find(deg1zcoord < 500);
inflowInd = deg1(inIndex);
outIndex = find(deg1zcoord >= 500);
outflowInd = deg1(outIndex);

nseg = length(lseg); 
nnod = numnodes(G); 

%% Show inflow and outflow of network
figure
H = plot(G,'XData',xyz(:,1),'YData',xyz(:,2),'ZData',xyz(:,3),...
    'LineWidth',0.5,'MarkerSize',2,'EdgeAlpha',1,'MarkerSize',0.001);
highlight(H,inflowInd,'MarkerSize',5,'Marker','s','NodeColor','g')
highlight(H,outflowInd,'MarkerSize',5,'Marker','o','NodeColor','r')
ax = gca;
ax.FontWeight = 'bold';
ax.FontName = 'Times New Roman';
ax.Title.String = 'Network with boundary nodes labeled (green-in,red-out)';
ax.Title.FontSize = 14;
ax.AmbientLightColor = 'magenta';
ax.XGrid = 'off';
ax.YGrid = 'off';
ax.ZGrid = 'off';
ax.XAxis.Label.String = 'x (\mum)';
ax.YAxis.Label.String = 'y (\mum)';
ax.ZAxis.Label.String = 'z (\mum)';
axis equal

%% Display node degrees
degreeList = degree(G);
fprintf("Degree 1 nodes: %i\n",length(find(degreeList == 1)))
fprintf("Degree 2 nodes: %i\n",length(find(degreeList == 2)))
fprintf("Degree 3 nodes: %i\n",length(find(degreeList == 3)))
fprintf("Degree 4 nodes: %i\n",length(find(degreeList == 4)))
fprintf("Degree 5 nodes: %i\n",length(find(degreeList == 5)))
fprintf("Degree 6 nodes: %i\n",length(find(degreeList == 6)))

%% Parameters
pIn = 50; % inflow pressure [mmHg]
pOut = 5; % outflow pressure [mmHg]
maxiter = 20000;
tol = 1e-2;

%% Initial values
q0 = 1e5; % initial volumetric flow [um^3/s]
u0 = 4e-3; % initial viscosity [Pa*s]
h0 = 0.45; % intial vessel hematocrit discharge
errh = 1;

%% Preallocation
J = zeros(nseg,1); % conductance (1/resistance)

%% Initial Vectors
q = ones(nseg,1)*q0;
u = ones(nseg,1)*u0/133.322; % [Pa*s] to [mmHg*s]
h = ones(nseg,1)*h0;

%% Fix zero diameter vessels
d(d==0) = mean(d);

%% Effects of Endothelial Surface Layer
[Deff, Dph] = ESLeffect(d,h);

%% P vector to fill with pressures
P = ones(nnod,1);
inflowInd = reshape(inflowInd,length(inflowInd),1);
outflowInd = reshape(outflowInd,length(outflowInd),1);
boundaryInd = [inflowInd;outflowInd];   
pUnknown = 1:nnod;
pUnknown(boundaryInd) = [];

%% Alternate between linear system and nonlinear system until convergence

iter = 0;

startmain = tic;

while (errh > tol) && (iter < maxiter)
    
    iter = iter + 1;
    
    hprev = h;
    
    % Define conductance of edges, units: [um^3/mmHg/s]
    J = pi.*d.^4./(128*u.*lseg);

    % setup matricies for Ax = b
    [A, b] = buildlinearsystem(G,J,nodnod,boundaryInd,...
        inflowInd,outflowInd,pIn,pOut);
    
    % remove rows and columns corresponding to boundary nodes
    A(boundaryInd,:) = [];
    A(:,boundaryInd) = [];
    b(boundaryInd) = [];
    
    % convert to sparse
    A = sparse(A);
    
%     solve linear system
    try
        X = A\b;
    catch
        X = pseudoinverse(A)*b;
    end

    % populate pressure vector
    P(pUnknown) = X;
    P(inflowInd) = pIn;
    P(outflowInd) = pOut;
    
    % solve volume flow rate, units: [um^3/s] or [fL/s]
    q = (P(edgenod(:,1)) - P(edgenod(:,2))).*J(1:nseg);
    
    % get negative flows on first iteration
    negFlow = find(q < 0);
        
    % reverse edge directions as needed
    if ~isempty(negFlow)
        
        edgenod = setflowdirection(edgenod,negFlow);
        
        % recalcuate flow
        q = (P(edgenod(:,1)) - P(edgenod(:,2))).*J(1:nseg);
        
    end
        
    % solve hematocrit
%     h = computehematocritPries(G,h,d,q,edgenod,boundaryInd);
    h = computehematocritGould(G,h,edgenod,d,5,q);
   
    % avoid extreme hematocrit
    hHigh = h >= 1;
    hLow = h <= 0;
    h(hHigh) = 0.95;
    h(hLow) = 0.1;
    
    % solve viscosity, units % [mPa*s] to [mmHg*s]
    u = computeviscosity(d,Deff,Dph,h)*1e-3/133.322;
    
    errh = abs(max(hprev - h));
    
    fprintf('\nIteration: %i\n',iter)
    fprintf('Hematocrit Error: %f\n',errh)
    fprintf('High hematocrit: %i\n',nnz(hHigh))
    fprintf('Low hematocrit: %i\n',nnz(hLow))
    
    % check conservation laws
    checkconservation(G,nodnod,J,P,boundaryInd)

    errhVector(iter) = errh;
end

endmain = toc(startmain);

%% Make Convergence plot
if iter > 1
    figure
    hplot = plot(1:length(errhVector),errhVector);
    hplot.LineWidth = 1.5;
    hplot.Color = 'k';
    ax = gca;
    ax.Title.String = 'Convergence Plot';
    ax.Title.FontSize = 14;
    ax.YLabel.String = 'Error';
    ax.YLabel.FontSize = 12;
    ax.XLabel.String = 'Iterations';
    ax.XLabel.FontSize = 12;
end

%% Make histogram for hematocrit
figure
hhist = histogram(h,'Normalization','probability');
hhist.FaceColor = 'k';
ax = gca;
ax.Title.String = 'Hematocrit Distribution';
ax.Title.FontSize = 14;
ax.YLabel.String = 'Proportion';
ax.YLabel.FontSize = 12;
ax.XLabel.String = 'Hematocrit';
ax.XLabel.FontSize = 12;
ax.XLim = [0 1];

%% Find velocity
crossSecArea = pi*(d/2).^2;
velocity = q./crossSecArea; % [um/s]
velocity = velocity*0.001; % [um/s] to [mm/s]

%% Total inflow and outflow (must be equal via conservation law)
qtotIn = 0;
for iIn = 1:length(inflowInd)
    bcnod = inflowInd(iIn);
    nodnear = nonzeros(nodnod(bcnod,:));
    edgenear = findedge(G,bcnod,nodnear);
    qtotIn = qtotIn + sum(q(edgenear));
end

qtotOut = 0;
for iOut = 1:length(outflowInd)
    bcnod = outflowInd(iOut);
    nodnear = nonzeros(nodnod(bcnod,:));
    edgenear = findedge(G,bcnod,nodnear);
    qtotOut = qtotOut + sum(q(edgenear));
end

%% Simulation Information
fprintf('\nTotal inflow: %f nL/min\n',qtotIn*1e-6*60);
fprintf('Total outflow: %f nL/min\n',qtotOut*1e-6*60);
fprintf('Minimum flow: %f nL/min\n',min(q)*1e-6*60)
fprintf('Maximum flow: %f nL/min\n',max(q)*1e-6*60)
fprintf('Mean flow: %f nL/min\n',mean(q)*1e-6*60)
fprintf('High hematocrit: %i\n',nnz(hHigh))
fprintf('Low hematocrit: %i\n',nnz(hLow))
fprintf('Elapsed Time: %f seconds\n',endmain)
fprintf('Iterations: %i\n' ,iter)
fprintf('Number of nodes: %i\n',nnod)
fprintf('Number of edges: %i\n',nseg)

q = q*1e-6*60; % [fL/s] to [nL/min]

%% Visualization
variable = q;
% variable = velocity;
numfig =  findobj('type','figure');
nfig = length(numfig);
climits = [min(variable) max(variable)];
visualize(G,variable,xyz,nfig+1,'Blood Flow Distribution','Flow',climits,2)

variable = h;

numfig =  findobj('type','figure');
nfig = length(numfig);
climits = [0.1 0.95];
visualize(G,variable,xyz,nfig+1,'Hematocrit Distribution','Hematocrit',climits,2)
