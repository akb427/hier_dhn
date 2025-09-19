function [pts] = partition_edges(G,w,e,n,flag)
%PARTITION_EDGES  Partitions edges based on weighted adjacency matrix to
%minimize modularity. 
%
%   [pts] = PARTITION_EDGES(G,w,e,n,flag)
%
%   DESCRIPTION:
%   
%
%   INPUTS:
%       in1  - Description of input 1 (type, format, units if applicable)
%       in2  - Description of input 2
%
%   OUTPUTS:
%       out1 - Description of output 1 (what it represents)
%       out2 - Description of output 2
%
%   DEPENDENCIES:
%
%   SEE ALSO:

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
e_miss = setdiff(Grm.Edges.Names, e_inc);
[e_miss_idx,~] = find(e_miss'==G.Edges.Names);
for i = 1:numel(e_miss)
    n_in = G.Edges.EndNodes(e_miss_idx(i),1);
    n_out = G.Edges.EndNodes(e_miss_idx(i),2);
    if any(e_miss(i)==e.hot)
        if any(n_out==pts_idx{1})
            pts_idx{1}=[pts_idx{1} n_in];
        else
            pts_idx{2}=[pts_idx{2} n_in];
        end
    else
        if any(n_in==pts_idx{1})
            pts_idx{1}=[pts_idx{1} n_out];
        else
            pts_idx{2}=[pts_idx{2} n_out];
        end
    end
end
pts_idx = cellfun(@unique,pts_idx,'UniformOutput',false);
for i = 1:2
    pts{i} = G.Nodes.Names(pts_idx{i})';
end


end
