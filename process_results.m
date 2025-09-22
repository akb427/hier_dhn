%PROCESS_RESULTS  Processes the hierarchical optimization results.
%
%   DESCRIPTION: Loads the step level solutions and combines it into a
%   single structure. Solves the flow and pressure to meet nominal demand
%   while minimizing mass flow rate in the bypass segments. Calculates cost
%   metrics and plots results. 
%   
%   INPUTS:
%       opti_variables.mat  - File of problem variables.
%       vLL_step  - Folder of stepwise solutions.
%
%   DEPENDENCIES: nominal_flow, nominal_flow_slack, fig_cost, fig_demand,
%   fig_graph, fig_intQ, fig_mI, fig_mdot_u, fig_part, fig_tempprofiles,
%   fig_graphred
%
%   SEE ALSO: create_variables, solve_flex

%% Load problem variables
clc, clear, close all

pth = string(pwd);
load(pth+filesep+"opti_variables.mat","G","e","v","n","params","params_all","vec_dP","sG","Gred","temp_prof")

%% Load and Combine simulation steps

pth_save = pth+filesep+"vLL_step"+filesep;
flname = pth_save+"vLLi_"+num2str(1);
load(flname,"vsim_i");
vsim_all = vsim_i;

% for all the files in the folder
for idx_solve = 2:length(dir(pth_save+"*.mat"))
    flname = pth_save+"vLLi"+num2str(idx_solve);
    load(flname, "vsim_i");
    vsim_all.dP_e = [vsim_all.dP_e vsim_i.dP_e];
    vsim_all.P_nd = [vsim_all.P_nd vsim_i.P_nd];
    vsim_all.mdot_e = [vsim_all.mdot_e vsim_i.mdot_e];
    vsim_all.T = [vsim_all.T vsim_i.T];
    vsim_all.Qp = [vsim_all.Qp vsim_i.Qp];
    vsim_all.intQ = [vsim_all.intQ vsim_i.intQ];
    vsim_all.mI = [vsim_all.mI vsim_i.mI];
    vsim_all.cost = vsim_all.cost+vsim_i.cost;
end

%% Nominal Full Network Cost

% Timing parameters
n_nom = n;
n_nom.seg = 4;
n_nom.step = n_nom.seg+1;
n_nom.seg_T = n_nom.seg*(params.dt/params.dt_T);
n_nom.step_T = n_nom.seg_T+1;

n_cycle = floor(params.tf/(n_nom.seg*params.dt));
vnom = cell(1,n_cycle);

% Nominal flow function
Mnom = nominal_flow_tfn(G,n,v,e,params,0);
Mnom_slack = nominal_flow_tfn(G,n,v,e,params,1);

% First cycle
params_nom.T0 = params.Ts(1)*ones(n.nu,1);
params_nom.T0(e.uo_idxnu,1) = params.TsetR;
idx_solve=1;
idx_i = (idx_solve-1)*n_nom.seg+1:idx_solve*n_nom.seg;
idx_i_T = (idx_solve-1)*n_nom.seg_T+1:idx_solve*n_nom.seg_T;
params_nom.Qb = params_all.Qb(:,idx_i);
params_nom.Tamb = params.Tamb(:,idx_i_T);
params_nom.Ts = params.Ts(:,idx_i_T);

% Solve nominal flow problem
vnom{1} = Mnom.call(params_nom);
if vnom{1}.vld
    params_nom.i_mdot_e = vnom{1}.mdot_e;
    params_nom.i_mI = vnom{1}.mI;
    params_nom.i_zeta_u = vnom{1}.zeta_u;
    params_nom.i_dP_e = vnom{1}.dP_e;
    params_nom.i_P_n = vnom{1}.P_n;
    vnom{1} = Mnom_slack.call(params_nom);
else
    error('Could not solve nominal flow')
end

%% Rest
for idx_solve = 2:n_cycle
    disp(idx_solve)
    idx_i = (idx_solve-1)*n_nom.seg+1:idx_solve*n_nom.seg;
    idx_i_T = (idx_solve-1)*n_nom.seg_T+1:idx_solve*n_nom.seg_T;
    params_nom.T0 = vnom{idx_solve-1}.T(:,end);
    params_nom.Qb = params_all.Qb(:,idx_i);
    params_nom.Tamb = params.Tamb(:,idx_i_T);
    params_nom.Ts = params.Ts(:,idx_i_T);
    vnom{idx_solve} = Mnom_slack.call(params_nom);
    disp(sum(vnom{idx_solve}.Qp_slack.^2,'all'))
    if ~vnom{idx_solve}.vld
        disp(idx_solve);
    end
end

%% Store Nominal in single vector

vnom_all = vnom{1};
for idx_solve = 2:n_cycle
    vnom_all.dP_e = [vnom_all.dP_e vnom{idx_solve}.dP_e];
    vnom_all.P_nd = [vnom_all.P_nd vnom{idx_solve}.P_nd];
    vnom_all.mdot_e = [vnom_all.mdot_e vnom{idx_solve}.mdot_e];
    vnom_all.T = [vnom_all.T vnom{idx_solve}.T];
    vnom_all.Qp = [vnom_all.Qp vnom{idx_solve}.Qp];
    vnom_all.mI = [vnom_all.mI vnom{idx_solve}.mI];
    vnom_all.cost = vnom_all.cost+vnom{idx_solve}.cost;
end

%% Cost Calculations
e_fin = inedges(G,v.term_idx);
e_fin_nu = find(any(e.nu_idx==e_fin));
cost_T_nom = 0;
cost_T_sim = 0;
idx = 0;
for idx_solve = 1:size(vsim_all.T,2)
    if mod((idx_solve-1)*params.dt_T,params.dt)==0      % when mdot can change
        idx = idx+1;
    end
    cost_T_nom = cost_T_nom+(vnom_all.mdot_e(e_fin,idx).*vnom_all.T(e_fin_nu,idx_solve))'*(params.w_T*eye(numel(e_fin)))*(vnom_all.mdot_e(e_fin,idx).*vnom_all.T(e_fin_nu,idx_solve));
    cost_T_sim = cost_T_sim+(vsim_all.mdot_e(e_fin,idx).*vsim_all.T(e_fin_nu,idx_solve))'*(params.w_T*eye(numel(e_fin)))*(vsim_all.mdot_e(e_fin,idx).*vsim_all.T(e_fin_nu,idx_solve));
end


%% Reduction Calc
mI_sim = params.dt*sum(vsim_all.mI,'all');
mI_nom = params_nom.dt*sum(vnom.mI,'all');
(mI_sim-mI_nom)/mI_nom

mby_sim = params.dt*sum(vsim_all.mdot_e(e.by,:),'all');
mby_nom = params_nom.dt*sum(vnom.mdot_e(e.by,:),'all');
(mby_sim-mby_nom)/mby_nom

%% Timing

tp.time = 0:params.dt:params.tf;
tp.time_T = 0:params.dt_T:params.tf;

tp.date_start = datetime(2018,1,28);
tp.date_end = datetime(2018,2,4);
tp.delta_date = seconds(tp.date_end-tp.date_start);

tp.time_day = interp1([0, params.tf], [tp.date_start tp.date_end],tp.time);
tp.time_T_day = interp1([0, params.tf], [tp.date_start tp.date_end],tp.time_T);

%% Figures

addpath(pth+filesep+"figures")

fig_graph(G,n,e,v)
fig_part(G,Gred,n,e,sG)
%fig_graphrd(Gred,n)
fig_demand(params,params_all,e,n,tp,1003)
fig_intQ(vsim_all,params,params_all,e,n)
fig_mdot_u(vsim_all,params,e,n)
fig_cost(vLL{12*6},vec_dP,n)
fig_mI(vsim_all, vnom,params,e)
fig_tempprofiles(tp,temp_prof)