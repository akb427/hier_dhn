function [A,B,E] = graph2ss(G,params,sparams,n,v,e,mdot_e,ismat)
%GRAPH2SS Convert graph, parameters and mass flow rates to state space
%matrices
%
%   [A,B,E] = GRAPH2SS(G,params,sparams,n,nd,e,mdot_e,ismat)
%
%   DESCRIPTION: Converts the network parameter and current edge mass flow
%   rates into the variables needed for the continuous time state space
%   calculation. Uses either matrices with known mdot_e or CasADi variables
%   depending on ismat value. 
%
%   INPUTS:
%       G       - Graph of network.
%       params  - Structure of network parameters.
%       n       - Strucutre of sizes.
%       mdot_e  - Vector of mass flow rate in edges.
%       ismat   - Binary indicator of numeric (true) or CasAdi (false).
%
%   OUTPUTS:
%       A   - State transition matrix.
%       B   - Supply temperature input matrix. 
%       E   - Return and ambient input matrix.

import casadi.*

%% Coefficients
if ismat
    c = zeros(n.e,3);
else
    c = MX(n.e,3);
end
c(e.nu_idx,1) = mdot_e(e.nu_idx).*sparams.pipes(e.nu_idx,4);
c(e.nu_idx,2) = (params.h*pi*sparams.pipes(e.nu_idx,1).*sparams.pipes(e.nu_idx,2)).*sparams.pipes(e.nu_idx,4)./params.cp;
c(e.nu_idx,3) = -(c(e.nu_idx,1)+c(e.nu_idx,2));

%% A matrix
edges = table2array(G.Edges);

if ismat
    A = zeros(n.nu);
else
    A = MX(n.nu,n.nu);
end

for i = 1:n.nu
    i2 = e.nu_idx(1,i);
    A(i,i) = c(i2,3);
    in_node = edges(i2,1);
    in_edge_all = find(edges(:,2)==in_node);
    in_edge = in_edge_all(all(in_edge_all~=e.u_idx,2))';
    if numel(in_edge)>0
        for j = in_edge
            j2 = find(j==e.nu_idx);
            if isscalar(in_edge_all)
                A(i,j2) = c(i2,1);
            else
                A(i,j2) = mdot_e(j)/mdot_e(i2)*c(i2,1);
            end
        end
    end
end

%% B
% Preallocate B
if ismat
    B = zeros(n.nu,2);
else
    B = MX(n.nu,2);
end

% Populate B
e.supply = find(edges(:,1)==v.root_idx);
for i = 1:n.nu
    i2 = e.nu_idx(1,i);
    if any(i2==e.uo_idx)
        in_node = edges(i2,1);
        in_edge_all = find(edges(:,2)==in_node);
        in_edge = in_edge_all(any(in_edge_all==e.u_idx,2))';
        B(i,2) = mdot_e(in_edge)/mdot_e(i2)*c(i2,1);
    elseif any(i2==e.supply)
        B(i,1) = c(i2,1);
    end     
end

%% E
% Preallocate E
if ismat
    E = zeros(n.nu,1);
else
    E = MX(n.nu,1);
end

% Populate E
for i = 1:n.nu
    i2 = e.nu_idx(1,i);
    E(i,1) = c(i2,2);
end

end