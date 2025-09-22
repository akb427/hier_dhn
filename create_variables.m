%CREATE_VARIABLES  Generates the hierarchical control problem.
%
%   DESCRIPTION: Takes the prescribed graph and partitions it into
%   subgraphs. Additionally stores the timing and network parameters for
%   these subgraphs. Uses these to generate the CasADi functions to be
%   solved in the hierarchical structure for the subgraphs.
%
%   OUTPUTS:
%       opti_variables.m - File storing control problem
%
%   DEPENDENCIES: generate_params, partition_edges, subgraph_params,
%   opt_low_tfn
%
%   SEE ALSO: solve_flex

%% Setup

clc, clear, close all

pth = pwd;
addpath(pth+append(filesep+"functions"))

%% Problem setup

G = digraph([1 1 1 1 2 2 3 3 4 4 5 5 6 6 7 8 9 10 11 11 12 12 13 13 14 14 15 16 17 18 19 21 21 22 23 23 23 24 24 25 26 26 27 28 28 29 30 31 31 32 32 33 33 34 34 35 35 36 36 37 37 38 39 40 41 42 43 44],...
    [2 21 23 31 3 11 4 10 5 9 6 8 7 7 8 9 10 19 12 18 13 17 14 16 15 15 16 17 18 19 20 22 22 20 24 26 28 25 25 30 27 27 30 29 29 30 20 32 44 33 43 34 35 41 41 36 40 37 39 38 38 39 40 42 42 43 44 20]);

e.u = [8 10 12 20 22 24 49 51 57 59]; % users without direct bypass connections

params.tf = (7*24+1)*60*60;

% Mass flow timing parameters
params.tf_opt = 1*60*60;       % [s]
params.dt = 10*60;
n.seg = params.tf_opt/params.dt;
n.step = n.seg+1;
% Temperature timing parameters
params.dt_T = 30;
n.seg_T = params.tf_opt/params.dt_T;
n.step_T = n.seg_T+1;

% Simulation timing parameters
params.tf_sim = 10*60;
n.seg_sim = params.tf_sim/params.dt;
n.step_sim = n.seg_sim+1;
n.seg_T_sim = params.tf_sim/params.dt_T;
n.step_T_sim = n.seg_T_sim+1;

[G,n,e,v,params,temp_prof] = generate_params(G,n,e, params);

% Cost weights
params.w_T = 10^-4;

% Pressure Drops
n_step = (params.tf-params.tf_opt)/params.tf_sim;
vec_dP = [.1:.1:1 1.5:.5:4.5 5:5:10];
params.dP_max = max(vec_dP);
n.dP = numel(vec_dP);

%% Subgraphs with manual changes

[pts] = partition_edges(G,params.pipes(:,4),e,n,0);
n.sg = 2;
[sG,se,sn,~,sparams] = subgraph_params(G,pts,e,n,params);

% Partition 2
[pts_n] = partition_edges(sG{2},sparams{2}.pipes(:,4),se{2},sn{2},0);
pts = [pts_n, pts(1)];
n.sg = 3;
[sG,se,sn,~,sparams] = subgraph_params(G,pts,e,n,params);

% Partition 2
[pts_n] = partition_edges(sG{2},sparams{2}.pipes(:,4),se{2},sn{2},0);
pts = [pts_n, pts([1 3])];
n.sg = 4;
[sG,se,sn,~,sparams] = subgraph_params(G,pts,e,n,params);

% Partition 1
[pts_n] = partition_edges(sG{1},sparams{1}.pipes(:,4),se{1},sn{1},0);
pts = [pts_n, pts([2 3 4])];
pts{1}(pts{1}==1) = [];
pts{1}(pts{1}==20) = [];
n.sg = 5;
[sG,se,sn,snd,sparams, Gred] = subgraph_params(G,pts,e,n,params);

%% Clean Up params
params_all = rmfield(params,{'tf','tf_opt','dt','dt_T','tf_sim','Inc','p','cp','h','mI','pipes','Ts','TsetR','Tamb','w_T','dP_max'});
params = rmfield(params,{'Cap_l','Cap_u','w_flex','Qb'});

%% Setup Optimization

vec_par = [repelem(1:n.sg,n.dP);repmat(vec_dP,1,n.sg)];
n.par = size(vec_par,2);
vLL_i = cell(1,n.par);
t_elapsed_i = zeros(1,n.par);

%% Create functions
M = cell(1,n.sg);
Ms = cell(1,n.sg);
e.fin = zeros(1,n.sg);
for i = 1:n.sg
    [M{i},e.fin(i)] = opt_low_tfn(sG{i},sn{i},snd{i},se{i},params,sparams{i},false);
    Ms{i} = opt_low_tfn(sG{i},sn{i},snd{i},se{i},params,sparams{i},true);
end
M = repelem(M,n.dP);
Ms = repelem(Ms,n.dP);

%% Store Variables

save('opti_variables')