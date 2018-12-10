function getfullreduced(data,segperstrand)
% GETFULLREDUCED Extracts information from a reduced full kleinfeld network.
%
% Input:
%   data - Data structure of a kleinfeld network.
%   segperstrand - Desired number of segments per strand.
%
% Output:
%   matfile - Data structure 'o' with network information as fields.
%
% Example of usage:  
%   getfullreduced(au,2)
%   getfullreduced(au,10)   
%
% Example of output:
%   getfullreduced(au,2) yields 'au_reduced_nps_3' matfile
%
% Note: 'nps' in the name of the saved matfile means 'nodes per strand'
%
% Note: input data structure should be stored as a variable that contains
%       a field called "vectorizedStructure".

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

%% Get strand list
if strcmp(name,'au')
    structcol = 3;
else
    structcol = 4;
end     
   
temp = struct2cell(network.vectorizedStructure.Strands.').';
strandList = struct('StartToEndIndices', temp(:,structcol));
clear temp

%% Loop strand segemnts to remove nodes
f = waitbar(0,'Reducing strand nodes...');

nstrand = length(strandList);
newStrand = cell(nstrand,1);
xyz_new = zeros(length(xyz),3);
r_new = zeros(length(xyz),1);

for istrand = 1:nstrand
    strandNodList = uint32(cell2mat(struct2cell(strandList(istrand))))';
    if length(strandNodList) ~= (segperstrand + 1) % strand reduction needed
        indJumpSize = ceil(length(strandNodList)/segperstrand);
        nodToKeep = strandNodList(1:indJumpSize:end);
        if ~isequal(segperstrand,length(nodToKeep))
            error('New strand length mismatch on strand: %i',istrand)
        end
        strandNodList = [nodToKeep;strandNodList(end)];
    end
    xyz_new(strandNodList,:) = xyz(strandNodList,:);
    rnodList = rnod(strandNodList);
    rnodList(rnodList==0) = mean(rnod);
    r_new(strandNodList) = rnodList;
    newStrand{istrand} = strandNodList;
    waitbar(istrand/nstrand,f)
end
close(f)

xyz_new(any(xyz_new==0,2),:) = [];
r_new(any(r_new==0,2),:) = [];

%% Adjacency Matrix
f = waitbar(0,'Generating Sparse Adj. Matrix... Please wait...');

nnod = length(xyz_new);
row = zeros(nnod,1);
col = zeros(nnod,1);
val = zeros(nnod,1);

nstrand = length(newStrand);
nodcount = 0;
for istrand = 1:nstrand
    nodes = cell2mat(newStrand(istrand));
    for ii = 1:length(nodes)-1
        if isempty(nodes) == 1
            break
        end
        nodcount = nodcount + 1;
        inod = nodes(ii);
        next = nodes(ii+1);
        row(nodcount) = inod;
        col(nodcount) = next;
        val(nodcount) = norm(xyz(inod,:)-xyz(next,:));
    end
end
close(f)

%% Rename nodes and grab correct xyz
rowlen = length(row);
[~,~,rename] = unique([row;col]);
row_rename = rename(1:rowlen);
col_rename = rename(rowlen+1:end);

Adj = sparse(row_rename,col_rename,val);
Adj = (Adj + Adj')/2;
G = graph(Adj);

%% Get nodes for penetrating arterioles and venuoles
artFullNodePath = [];
for i = 1:length(artFullStrandPath)
    istrand = artFullStrandPath(i);
    nodes = cell2mat(newStrand(istrand));
    artFullNodePath = [artFullNodePath;nodes]; 
end

venFullNodePath = [];
for i = 1:length(venFullStrandPath)
    istrand = venFullStrandPath(i);
    nodes = cell2mat(newStrand(istrand));
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

%% Specify Boundary Nodes
% inflowInd = [];
% outflowInd = [];
% 
% for istrand = 1:nstrand
%     nodes = cell2mat(newStrand(istrand));    
%     if ismember(istrand,artTip)
%         inflowInd = [inflowInd;nodes(1)];
%     end
%     if ismember(istrand,venTip)
%         outflowInd = [outflowInd;nodes(1)];
%     end
% end

%% Rename nodes numbers in artioles, venules, and strands
list1 = [sort([row;col]),sort([row_rename;col_rename])];
list2 = unique(list1,'rows');

for in = 1:length(inflowInd)
    ind = list2(:,1) == inflowInd(in);
    inflowInd_rename(in,1) = list2(ind,2);
end
    
for out = 1:length(outflowInd)
    ind = list2(:,1) == outflowInd(out);
    outflowInd_rename(out,1) = list2(ind,2);
end

for inod = 1:length(artFullNodePath)
    ind = list2(:,1) == artFullNodePath(inod);
    artFullNodePath_rename(inod,1) = list2(ind,2);
end

for inod = 1:length(venFullNodePath)
    ind = list2(:,1) == venFullNodePath(inod);
    venFullNodePath_rename(inod,1) = list2(ind,2);
end

newStrandRenamed = zeros(nstrand,3);
for istrand = 1:nstrand
    nodes = cell2mat(newStrand(istrand));
    for i = 1:length(nodes)
        ind = list2(:,1) == nodes(i);
        nodes(i) = list2(ind,2);
    end
    newStrandRenamed(istrand,:) = nodes;
end

%% Compute segment diameter
nseg = numedges(G);

if segperstrand == 1
    d = network.vectorizedStructure.effectiveStrandRadiusList*2;
else
    d = zeros(nseg,1);
    edgenod = G.Edges{:,1};
    nodsta = edgenod(:,1);
    nodend = edgenod(:,2);
    for iseg = 1:nseg
        d(iseg) = (r_new(nodsta(iseg)) + r_new(nodend(iseg))); % average
    end
end

%% Grab segment length
lseg = G.Edges{:,2}*2;

% lseg = zeros(nseg,1);
% f = waitbar(0,'Calculating segment lengths...');
% for iseg = 1:nseg
%     [nod1, nod2] = findedge(G,iseg);
%     lseg(iseg) = norm(xyz(nod2,:)-xyz(nod1,:));
%     if lseg(iseg) == 0
%         error('Length of segment %i is zero',iseg)
%     end
%     waitbar(iseg/nseg,f)
% end
% close(f)

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
h = plot(G,'XData',xyz_new(:,1),'YData',xyz_new(:,2),'ZData',xyz_new(:,3),...
    'LineWidth',0.5,'EdgeAlpha',0.2,'MarkerSize',0.001);
highlight(h,artFullNodePath_rename,'MarkerSize',3,'Marker','s','NodeColor','g')
highlight(h,venFullNodePath_rename,'MarkerSize',3,'Marker','o','NodeColor','r')
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

%% Output data structure
o.xyz = xyz_new;
o.Adj = Adj;
o.d = d;
o.lseg = lseg;
o.name = name;
o.edgenod = edgenod;
o.nodnod = nodnod;
o.inflowInd = inflowInd_rename;
o.outflowInd = outflowInd_rename;
o.artFullNodePath = artFullNodePath_rename;
o.venFullNodePath = venFullNodePath_rename;
o.strandList = newStrandRenamed;

%% Save
nps = length(strandNodList);
matpath = 'G:\My Drive\matlab\research\kleinfeld_processing\output\matfiles\reduced\';
figpath = 'G:\My Drive\matlab\research\kleinfeld_processing\output\figures\vanilla\reduced\';
figname = sprintf([name '_reduced_nps_' num2str(nps)]);
matname = sprintf([name '_reduced_nps_' num2str(nps)]);
saveas(gcf,[figpath figname])
save([matpath matname],'o')

%#ok<*AGROW>
%#ok<*STRNU>