function [sG,se,sn,snd,sparams,Gred] = subgraph_params(G,pts,e,n,params)
%GENERATE_SUBGRAPH Summary of this function goes here
%   Detailed explanation goes here

%% Preallocate storage
sG = cell(n.sg,1);
se = cell(n.sg,1);
sn = cell(n.sg,1);
snd = cell(n.sg,1);
sparams = cell(n.sg,1);
%% Divide Graph
for i = 1:n.sg
    sG{i} = subgraph(G,pts{i});
    se{i}.e = sG{i}.Edges.Names';
    se{i}.u = intersect(se{i}.e,e.u);
    se{i}.by = intersect(se{i}.e,e.by);
    se{i}.hot = intersect(se{i}.e,e.hot);
    se{i}.cold = intersect(se{i}.e,e.cold);
    
    % index values instead of names
    [~,se{i}.u_idx] = find(se{i}.u'==se{i}.e);
    [~,se{i}.by_idx] = find(se{i}.by'==se{i}.e);
    se{i}.u_idx = se{i}.u_idx';
    se{i}.by_idx = se{i}.by_idx';
    
    %% Subgraph elements
    
    sn{i} = n;
    % number of elements
    sn{i}.e = numedges(sG{i});
    sn{i}.v = numnodes(sG{i});
    sn{i}.u = numel(se{i}.u);
    
    % root node
    snd{i}.root_idx = find(indegree(sG{i})==0);
    snd{i}.root = sG{i}.Nodes.Names(snd{i}.root_idx);
    
    % terminal node
    snd{i}.term_idx = find(outdegree(sG{i})==0);
    snd{i}.term = sG{i}.Nodes.Names(snd{i}.term_idx);
    
    % non-user edges
    se{i}.nu = setdiff(se{i}.e, se{i}.u);
    [~,se{i}.nu_idx] = find(se{i}.nu'==se{i}.e);
    se{i}.nu_idx = se{i}.nu_idx';
    sn{i}.nu = numel(se{i}.nu);
    
    %% In and out edges of users
    if ~isempty(se{i}.u)
        % inedges
        nodes = sG{i}.Edges.EndNodes(se{i}.u_idx,1)';
        se{i}.ui_idx = zeros(1,numel(nodes));
        for j = 1:numel(nodes)
            if ~isempty(inedges(sG{i},nodes(j)))
                se{i}.ui_idx(j) = inedges(sG{i},nodes(j));
            else
                se{i}.ui_idx(j) = NaN;
            end
        end
        
        % outedges
        nodes = sG{i}.Edges.EndNodes(se{i}.u_idx,2)';
        se{i}.uo_idx = zeros(1,numel(nodes));
        for j = 1:numel(nodes)
            if ~isempty(outedges(sG{i},nodes(j)))
                se{i}.uo_idx(j) = outedges(sG{i},nodes(j));
            else
                se{i}.uo_idx(j) = NaN;
            end
        end
        
        se{i}.ui = sG{i}.Edges.Names(se{i}.ui_idx)';
        se{i}.uo = sG{i}.Edges.Names(se{i}.uo_idx)';
        [se{i}.ui_idxnu,~] = find(se{i}.ui==se{i}.nu');
        [se{i}.uo_idxnu,~] = find(se{i}.uo==se{i}.nu');
        % se{i}.ui_idxnu = find(any(se{i}.nu==se{i}.ui',1));
        % se{i}.uo_idxnu = find(any(se{i}.nu==se{i}.uo',1));
    end
    %% Subgraph parameters
    
    sparams{i}.I = -incidence(sG{i});
    sparams{i}.pipes = params.pipes(se{i}.e,:);
    u_num = se{i}.u==e.u';
    sparams{i}.Qb = zeros(sn{i}.u,size(params.Qb,2));
    sparams{i}.Cap_l = zeros(sn{i}.u,size(params.Cap_l,2));
    sparams{i}.Cap_u = zeros(sn{i}.u,size(params.Cap_u,2));
    sparams{i}.w_flex = zeros(sn{i}.u);
    for j = 1:sn{i}.u
        sparams{i}.Qb(j,:) = params.Qb(u_num(:,j),:);
        sparams{i}.Cap_l(j,:) = params.Cap_l(u_num(:,j),:);
        sparams{i}.Cap_u(j,:) = params.Cap_u(u_num(:,j),:);
        sparams{i}.w_flex(j,j) = params.w_flex(u_num(:,j),u_num(:,j));
    end
end

%% Create reduced graph

edg_rm = cell2mat(cellfun(@(x) x.Edges.Names',sG,'UniformOutput',false)');
Gred = rmedge(G,edg_rm);
Gred.Edges.Names = NaN(numedges(Gred),1);
for i = 1:n.sg
    Gred = addedge(Gred,snd{i}.root,snd{i}.term);
    edg_idx = all(Gred.Edges.EndNodes == [snd{i}.root,snd{i}.term],2)&Gred.Edges.Names==0;
    Gred.Edges.Names(edg_idx) = i;
end

Gred = rmnode(Gred,find(indegree(Gred)==0&outdegree(Gred)==0));
