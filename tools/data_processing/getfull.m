function getfull(data)
% GETFULL Extracts information from a full kleinfeld network.
%
% Input:
%   data - Data structure of a kleinfeld network.
%
% Output:
%   matfile - Data structure 'o' with network information as fields.
%
% Example of usage:  
%   getfull(au)
%
% Example of output:
%   'getfull(au)' yields 'au_full' matfile
%
% Note: input data structure should be stored as a variable that contains
%       a field called "vectorizedStructure".
%
% Note: Save path is hard coded. Make sure to change the path if needed.

%% Load data
name = inputname(1);
network = data;
xyz = double(network.vectorizedStructure.Vertices.AllVerts);
rnod = network.vectorizedStructure.Vertices.AllRadii;
artTip = network.vectorizedStructure.plungingArterioleTipList';
venTip = network.vectorizedStructure.plungingVenuleTipList';
artFullStrandPath = network.vectorizedStructure.allPlungingArterioleList';
venFullStrandPath = network.vectorizedStructure.allPlungingVenuleList';

starttime = tic;
    
%% Adjacency Matrix
if strcmp(name,'au')
    structcol = 3;
else
    structcol = 4;
end     
   
temp = struct2cell(network.vectorizedStructure.Strands.').';
strandList = struct('StartToEndIndices', temp(:,structcol));
clear temp

nnod = length(xyz);
row = zeros(nnod,1);
col = zeros(nnod,1);
val = zeros(nnod,1);

f = waitbar(0,'Generating adjacency matrix...');

nstrand = length(strandList);
nodcount = 0;
for istrand = 1:nstrand
    nodes = cell2mat(struct2cell(strandList(istrand)))';
    for ii = 1:length(nodes)-1
        if isempty(nodes) == 1
            break
        end
        nodcount = nodcount + 1;
        if nodcount > nnod
            break
        end
        inod = nodes(ii);
        next = nodes(ii+1);
        row(nodcount) = inod;
        col(nodcount) = next;
        val(nodcount) = norm(xyz(inod,:)-xyz(next,:));
    end
end

close(f)

Adj = sparse(row,col,val,nnod,nnod);
Adj = (Adj + Adj')/2;
G = graph(Adj);

%% Get nodes for penetrating arterioles and venuoles
artFullNodePath = [];
for i = 1:length(artFullStrandPath)
    istrand = artFullStrandPath(i);
    nodes = cell2mat(struct2cell(strandList(istrand)))';
    artFullNodePath = [artFullNodePath;nodes]; 
end

venFullNodePath = [];
for i = 1:length(venFullStrandPath)
    istrand = venFullStrandPath(i);
    nodes = cell2mat(struct2cell(strandList(istrand)))';
    venFullNodePath = [venFullNodePath;nodes];
end

%% Specify boundary nodes
inflowInd = zeros(length(artTip),1);
outflowInd = zeros(length(venTip),1);
nInflow = length(artTip);
nOutflow = length(venTip);

for i = 1:nInflow
    istrand = artTip(i);
    nodes = cell2mat(struct2cell(strandList(istrand)))';
    inflowInd(i) = nodes(1);
end

for i = 1:nOutflow
    istrand = venTip(i);
    nodes = cell2mat(struct2cell(strandList(istrand)))';
    outflowInd(i) = nodes(1);
end

%% Compute segment diameter
nseg = numedges(G);
d = zeros(nseg,1);
edgenod = G.Edges{:,1};
nodsta = edgenod(:,1);
nodend = edgenod(:,2);

for iseg = 1:nseg
    d(iseg) = (rnod(nodsta(iseg)) + rnod(nodend(iseg))); % average
end

%% Get length from graph object
lseg = 2*G.Edges{:,2};

%% Node connectivity
f = waitbar(0,'Setting up node connectivity...');
nnod = numnodes(G);
nodnod = zeros(nnod,3); % nodes connected to nodes
for inod = 1:nnod
    nodind = find(Adj(:,inod) ~= 0);
    for jnod = 1:length(nodind)
        nod = nodind(jnod);
        nodnod(inod,jnod) = nod;
    end
    waitbar(inod/nnod,f);
end
close(f)

%% Time
endtime = toc(starttime);
fprintf('Elapsed time: %f seconds\n',endtime);

%% Plot
figure
h = plot(G,'XData',xyz(:,1),'YData',xyz(:,2),'ZData',xyz(:,3),...
    'LineWidth',0.5,'EdgeAlpha',0.2,'MarkerSize',0.001);
highlight(h,artFullNodePath,'MarkerSize',1,'Marker','s','NodeColor','g')
highlight(h,venFullNodePath,'MarkerSize',1,'Marker','o','NodeColor','r')
ax = gca;
ax.FontWeight = 'bold';
ax.FontName = 'Times New Roman';
ax.Title.String = 'Arteriole nodes: green || Venule nodes: red';
ax.Title.FontSize = 14;
ax.AmbientLightColor = 'magenta';
ax.XGrid = 'off';
ax.YGrid = 'off';
ax.ZGrid = 'off';
ax.XAxis.Label.String = 'x (\mum)';
ax.YAxis.Label.String = 'y (\mum)';
ax.ZAxis.Label.String = 'z (\mum)';
axis equal

%% Store
o.xyz = xyz;
o.Adj = Adj;
o.d = d;
o.lseg = lseg;
o.name = name;
o.edgenod = edgenod;
o.nodnod = nodnod;
o.inflowInd = inflowInd;
o.outflowInd = outflowInd;
o.artFullNodePath = artFullNodePath;
o.venFullNodePath = venFullNodePath;
o.strandList = strandList;

%% Save
matpath = 'G:\My Drive\matlab\research\kleinfeld_processing\output\matfiles\full\';
figpath = 'G:\My Drive\matlab\research\kleinfeld_processing\output\figures\vanilla\full\';
figname = sprintf([name '_full']);
matname = sprintf([name '_full']);
saveas(gcf,[figpath figname])
save([matpath matname],'o')

%#ok<*AGROW>
%#ok<*STRNU>