function [sG,se,sn,sv,sparams,Gred] = subgraph_params(G,pts,e,n,params)
%SUBGRAPH_PARAMS  One-line summary of what the function does.
%
%   [sG,se,sn,snd,sparams,Gred] = SUBGRAPH_PARAMS(G,pts,e,n,params)
%
%   DESCRIPTION: Takes the original graph and parameters and uses the pts
%   grouping to divide the existing problem information into cells. Also
%   finds the interconnections between the subgraphs and produces a reduced
%   graph where the subgraphs are nodes.
%
%   INPUTS:
%       G   - Graph of network.
%       pts - Cell of the two node groups.
%       e   - Structure of edge information.
%       n   - Structure of sizes.
%       params  - Structure of network component parameters.
%
%   OUTPUTS:
%       sG  - Cell of graphs of partitioned network.
%       se  - Cell of structures of partitioned edge information.
%       sn  - Cell of structures of partitioned sizes.
%       sv  - Cell of structures of partitioned node information.
%       sparams - Cell of structures of partitioned network componenet
%       parameters.
%       Gred    - Graph of subgrah interconnections.

%% Preallocate storage
sG = cell(n.sg,1);
se = cell(n.sg,1);
sn = cell(n.sg,1);
sv = cell(n.sg,1);
sparams = cell(n.sg,1);

%% Divide Graph
for idx_sg = 1:n.sg
    sG{idx_sg} = subgraph(G,pts{idx_sg});
    se{idx_sg}.e = sG{idx_sg}.Edges.Names';
    se{idx_sg}.u = intersect(se{idx_sg}.e,e.u);
    se{idx_sg}.by = intersect(se{idx_sg}.e,e.by);
    se{idx_sg}.hot = intersect(se{idx_sg}.e,e.hot);
    se{idx_sg}.cold = intersect(se{idx_sg}.e,e.cold);
    
    % index values instead of names
    [~,se{idx_sg}.u_idx] = find(se{idx_sg}.u'==se{idx_sg}.e);
    [~,se{idx_sg}.by_idx] = find(se{idx_sg}.by'==se{idx_sg}.e);
    se{idx_sg}.u_idx = se{idx_sg}.u_idx';
    se{idx_sg}.by_idx = se{idx_sg}.by_idx';
end

%% Subgraph elements
for idx_sg = 1:n.sg
    sn{idx_sg} = n;
    % number of elements
    sn{idx_sg}.e = numedges(sG{idx_sg});
    sn{idx_sg}.v = numnodes(sG{idx_sg});
    sn{idx_sg}.u = numel(se{idx_sg}.u);
    
    % root node
    sv{idx_sg}.root_idx = find(indegree(sG{idx_sg})==0);
    sv{idx_sg}.root = sG{idx_sg}.Nodes.Names(sv{idx_sg}.root_idx);
    
    % terminal node
    sv{idx_sg}.term_idx = find(outdegree(sG{idx_sg})==0);
    sv{idx_sg}.term = sG{idx_sg}.Nodes.Names(sv{idx_sg}.term_idx);
    
    % non-user edges
    se{idx_sg}.nu = setdiff(se{idx_sg}.e, se{idx_sg}.u);
    [~,se{idx_sg}.nu_idx] = find(se{idx_sg}.nu'==se{idx_sg}.e);
    se{idx_sg}.nu_idx = se{idx_sg}.nu_idx';
    sn{idx_sg}.nu = numel(se{idx_sg}.nu);
end

%% In and out edges of users
for idx_sg = 1:n.sg
    if ~isempty(se{idx_sg}.u)
        % inedges
        nodes = sG{idx_sg}.Edges.EndNodes(se{idx_sg}.u_idx,1)';
        se{idx_sg}.ui_idx = zeros(1,numel(nodes));
        for j = 1:numel(nodes)
            if ~isempty(inedges(sG{idx_sg},nodes(j)))
                se{idx_sg}.ui_idx(j) = inedges(sG{idx_sg},nodes(j));
            else
                se{idx_sg}.ui_idx(j) = NaN;
            end
        end
        
        % outedges
        nodes = sG{idx_sg}.Edges.EndNodes(se{idx_sg}.u_idx,2)';
        se{idx_sg}.uo_idx = zeros(1,numel(nodes));
        for j = 1:numel(nodes)
            if ~isempty(outedges(sG{idx_sg},nodes(j)))
                se{idx_sg}.uo_idx(j) = outedges(sG{idx_sg},nodes(j));
            else
                se{idx_sg}.uo_idx(j) = NaN;
            end
        end
        % other naming conventions
        se{idx_sg}.ui = sG{idx_sg}.Edges.Names(se{idx_sg}.ui_idx)';
        se{idx_sg}.uo = sG{idx_sg}.Edges.Names(se{idx_sg}.uo_idx)';
        [se{idx_sg}.ui_idxnu,~] = find(se{idx_sg}.ui==se{idx_sg}.nu');
        [se{idx_sg}.uo_idxnu,~] = find(se{idx_sg}.uo==se{idx_sg}.nu');
    end
end

%% Subgraph parameters
for idx_sg = 1:n.sg
    % Graph information
    sparams{idx_sg}.I = -incidence(sG{idx_sg});
    % Pipe information
    sparams{idx_sg}.pipes = params.pipes(se{idx_sg}.e,:);
    % User information
    u_num = se{idx_sg}.u==e.u';
    sparams{idx_sg}.Qb = zeros(sn{idx_sg}.u,size(params.Qb,2));
    sparams{idx_sg}.Cap_l = zeros(sn{idx_sg}.u,size(params.Cap_l,2));
    sparams{idx_sg}.Cap_u = zeros(sn{idx_sg}.u,size(params.Cap_u,2));
    sparams{idx_sg}.w_flex = zeros(sn{idx_sg}.u);
    for j = 1:sn{idx_sg}.u
        sparams{idx_sg}.Qb(j,:) = params.Qb(u_num(:,j),:);
        sparams{idx_sg}.Cap_l(j,:) = params.Cap_l(u_num(:,j),:);
        sparams{idx_sg}.Cap_u(j,:) = params.Cap_u(u_num(:,j),:);
        sparams{idx_sg}.w_flex(j,j) = params.w_flex(u_num(:,j),u_num(:,j));
    end
end

%% Create reduced graph

edg_rm = cell2mat(cellfun(@(x) x.Edges.Names',sG,'UniformOutput',false)');
Gred = rmedge(G,edg_rm);
Gred.Edges.Names = NaN(numedges(Gred),1);
for i = 1:n.sg
    Gred = addedge(Gred,sv{i}.root,sv{i}.term);
    edg_idx = all(Gred.Edges.EndNodes == [sv{i}.root,sv{i}.term],2)&Gred.Edges.Names==0;
    Gred.Edges.Names(edg_idx) = i;
end

Gred = rmnode(Gred,find(indegree(Gred)==0&outdegree(Gred)==0));
