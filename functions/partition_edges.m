function [pts] = partition_edges(G,w,e,n,flag)
%PARTITION_EDGES  Partitions edges based on weighted adjacency matrix to
%minimize modularity. 
%
%   [pts] = PARTITION_EDGES(G,w,e,n,flag)
%
%   DESCRIPTION: Uses the Fielder vector to partition the graph G into two
%   subgraphs based on the weights w. Because the pipes are the edges, this
%   function also adds in the removed edges to one of the subgraphs
%   empirically based on the existing node grouping and the missing edge
%   type.
%
%   INPUTS:
%       G       - Graph to be partitionined.
%       w       - Vector of edge weights for partitioning.
%       e       - Structure of edge information.
%       n       - Structure of sizes.
%       flag    - Binary indicating plotting of partition.
%
%   OUTPUTS:
%       pts - Cell of the two node groups.

%% Convert to Weighted Undirected Graph
Grm = rmedge(G,e.u_idx);
w(e.u_idx) = [];
Gu = graph(Grm.Edges);

%% Calculate Partition
A = adjacency(Gu,w);
D = diag(sum(A));
rD = sqrtm(inv(D));
L = eye(n.v)-rD*A*rD;
[V,~] = eigs(L,2,'smallestabs');
pts_idx{1} = find(V(:,2)>=0)';
pts_idx{2} = find(V(:,2)<0)';

% Plotting
if flag
    figure
    h = plot(G); 
    highlight(h,pts_idx{1},'NodeColor','r')
    highlight(h,pts_idx{2},'NodeColor','g')
end

%% Add in missing edges
G1 = subgraph(Gu,pts_idx{1});
G2 = subgraph(Gu,pts_idx{2});
e_inc = [G1.Edges.Names; G2.Edges.Names];
% identify missing edges
e_miss = setdiff(Grm.Edges.Names, e_inc);
[e_miss_idx,~] = find(e_miss'==G.Edges.Names);
% for all missing edges
for idx_e = 1:numel(e_miss)
    % identify start and end nodes
    n_in = G.Edges.EndNodes(e_miss_idx(idx_e),1);
    n_out = G.Edges.EndNodes(e_miss_idx(idx_e),2);
    % for the hot edges
    if any(e_miss(idx_e)==e.hot)
        % if the out node is in group one
        if any(n_out==pts_idx{1})
            % add the in node to group one
            pts_idx{1}=[pts_idx{1} n_in];
        else
            % add the in node to group two
            pts_idx{2}=[pts_idx{2} n_in];
        end
    else
        % if the in node is in group one
        if any(n_in==pts_idx{1})
            % add the out node to group one
            pts_idx{1}=[pts_idx{1} n_out];
        else
            % add the out node to group two
            pts_idx{2}=[pts_idx{2} n_out];
        end
    end
end
% remove duplicates
pts_idx = cellfun(@unique,pts_idx,'UniformOutput',false);
% replace index with node name
for idx_pts = 1:2
    pts{idx_pts} = G.Nodes.Names(pts_idx{idx_pts})';
end


end
