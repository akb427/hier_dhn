function [G,n,e,v,params,temp_prof] = generate_params(G,n,e,params)
%GENERATE_PARAMS Generate parameters for network components.
%
%   [G,n,e,v,params,temp_prof] = GENERATE_PARAMS(G,n,e,params)
%
%   DESCRIPTION: Classifies the elements in G. Contains hard-coded values
%   for the network parameters, pipe parameters, temperature profiles
%   followed by each building type. Additionally loads building heat
%   demands and ambient conditions for 1/28/2018-2/4/2018. Finally loads
%   the buildings heat capacity. 
%
%   INPUTS:
%       G       - Graph of network.
%       n       - Structure of sizes.
%       params  - Structure of network components parameter, including
%       pipes in the form [L D zetaB zetaF mdot Lbypass Dbypass mdotbypass]
%       and users in the form [mu Q Ls1 Ls2 Ls3 Ds1 Ds2 Ds3].
%
%   OUTPUTS:
%       G           - Graph of updated network.
%       n           - Structure of updated sizes.
%       e           - Structure of edge information.
%       v           - Structure of node information.
%       params      - Structure of update network component parameters.
%       temp_prof   - Structure of allowable temperature diviations.
%
%   DEPENDENCIES: heat profiles.mat, USA_IL_Chicago-Midway.AP_Tamb.xlsx,
%       Cbuild.mat 

%% Data loading

pth_data = pwd+filesep+"data"+filesep; 

%% Timing

time = 0:params.dt:params.tf;
time_T = 0:params.dt_T:params.tf;

date_start = datetime(2018,1,28);
date_end = datetime(2018,2,4);
delta_date = seconds(date_end-date_start);

%% Graph 

% Graph elements
params.Inc = -incidence(G);
n.v = numnodes(G);
n.e = numedges(G);
G.Edges.Names = (1:n.e)';
G.Nodes.Names = (1:n.v)';

% Edge set
edges = [G.Edges.Names G.Edges.EndNodes];

% Bypass edges
[~,ia,~] = unique(edges(:,2:3),'rows');
e.by = setdiff(1:n.e,ia);
e.by_idx = e.by;
% Users connected to bypass
for i = e.by
    u = find(all(edges(i,[2,3])==edges(:,[2 3]),2));
    e.u = [e.u u(u~=i)];
end

% Overall user information
e.u_idx = e.u;
e.nu = setdiff(1:n.e,e.u);
e.nu_idx = e.nu;
n.u = numel(e.u);
n.nu = numel(e.nu);

% Plant information
v.root=find(indegree(G)==0);
v.term = find(outdegree(G)==0);
v.root_idx = v.root;
v.term_idx = v.term;

% Hot and cold elements
vhot = find(outdegree(G)>1)';
vcold = find(indegree(G)>1)';

G_cold = flipedge(subgraph(G,vcold));
G_hot = subgraph(G, vhot);
P = isomorphism(G_hot,G_cold);
G_cold = reordernodes(G_cold,P);
e.hot = G_hot.Edges.Names';
e.cold = G_cold.Edges.Names';

%% Inedges and outedges of users
% Inedges
nodes = G.Edges.EndNodes(e.u_idx,1)';
e.ui_idx = zeros(1,numel(nodes));
for i = 1:numel(nodes)
    if ~isempty(inedges(G,nodes(i)))
        e.ui_idx(i) = inedges(G,nodes(i));
    else
        e.ui_idx(i) = NaN;
    end
end

% Outedges
nodes = G.Edges.EndNodes(e.u_idx,2)';
e.uo_idx = zeros(1,numel(nodes));
for i = 1:numel(nodes)
    if ~isempty(outedges(G,nodes(i)))
        e.uo_idx(i) = outedges(G,nodes(i));
    else
        e.uo_idx(i) = NaN;
    end
end

% Index to name
e.ui = G.Edges.Names(e.ui_idx)';
e.uo = G.Edges.Names(e.uo_idx)';
% Index in nonuser edge list
[e.ui_idxnu,~] = find(e.ui==e.nu');
[e.uo_idxnu,~] = find(e.uo==e.nu');

%% Network parameters

params.p = 971;
params.cp = 4179;
params.h = 1.5;
params.mI = n.u*5;
params.Ts = 80*ones(1,numel(time_T));
params.TsetR = 30;

%% Building parameters

% Flexibible hours
idx_day = floor(time_T/(24*60*60))+1;
idx_hr = mod(floor(time_T/3600),24);

t_all = false(n.u,length(time_T));

t_res = false(size(time_T));
t_res(idx_day~=1&idx_day~=7&idx_hr<18&idx_hr>9)=1;
nr = find(ismember(e.u,[8 10 12 13 20 22 24 25 57 59 60]))';
t_all(nr,:) = repmat(t_res,numel(nr),1);

t_com = false(size(time_T));
t_com(~(idx_day~=1&idx_day~=7&idx_hr<18&idx_hr>6)) = 1;
nr = find(ismember(e.u,[49 51]))';
t_all(nr,:) = repmat(t_com,numel(nr),1);

t_retail = false(size(time_T));
t_retail(idx_hr>22|idx_hr<6)=1;
nr = find(ismember(e.u,[38 41 44 54]))';
t_all(nr,:) = repmat(t_retail,numel(nr),1);

% Store Profiles
temp_prof.res = t_res;
temp_prof.com = t_com;
temp_prof.retail = t_retail;
temp_prof.med = false(size(t_res));

% Temperature Limits
T_all = zeros(size(t_all));
T_all(t_all) = 2;
T_all(~t_all) = 2;

% Get Capacity
load(pth_data+"Cbuild.mat",'pcpV');

params.Cap_l = pcpV.*-T_all/10^6; %MJ
params.Cap_u = pcpV.*T_all/10^6;  %MJ

params.w_flex = diag(1./pcpV*10^3);

%% Pipe Parameters
lambda = 0.01;

% pipes = [L(1) D(2) zeta(3) 1/pV(4)]
params.pipes = zeros(n.e, 3);
params.pipes(e.hot,1) = [60 100 50 30 15 20 10 10 10 10 10 10 10 15 20 30 45 15 20 10 10];
params.pipes(e.by,1) = 3;
params.pipes(e.cold,1) = params.pipes(e.hot,1);
params.pipes(:,2) = .40*ones(n.e,1);
params.pipes(e.by,2) = .15;
params.pipes(e.u,2) = 0;
params.pipes(:,3) = lambda*params.pipes(:,1)./params.pipes(:,2)./(2*params.p*(pi/4*params.pipes(:,2).^2).^2);
params.pipes(:,4) = 1./params.p./(pi/4.*params.pipes(:,2).^2.*params.pipes(:,1));

%% Heat Demand

load(pth_data+"demand profiles",'Qb','t');
Qb = Qb(:,t>=date_start & t<=date_end);
time_Qb_orig = 0:15*60:delta_date;
params.Qb = interp1(time_Qb_orig, Qb', time)'; %kW

%% Ambient Temperature

Tamb = readmatrix(pth_data+"USA_IL_Chicago-Midway.AP_Tamb.xlsx")';
time_amb_orig = datetime(2018,1,1,0,0,0)+seconds(Tamb(1,:));
Tamb= Tamb(2,time_amb_orig>=date_start & time_amb_orig<=date_end);
time_amb_orig = 0:60*60:delta_date;
params.Tamb = interp1(time_amb_orig, Tamb, time_T);

end