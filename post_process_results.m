%% Load Results

%% Load and Combine simulation steps
flname = pwd+"\vLL notime1\"+"vLL_nt_"+num2str(1)+".mat";
load(flname);
vsim_all = vsim_i;
for i = 2:37
    flname = pwd+"\vLL notime1\"+"vLL_nt_"+num2str(i);
    load(flname);
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

% First cycle
params_nom.T0 = params.Ts(1)*ones(n.nu,1);
params_nom.T0(e.uo_idxnu,1) = params.TsetR;
i=1;
idx_i = (i-1)*n_nom.seg+1:i*n_nom.seg;
idx_i_T = (i-1)*n_nom.seg_T+1:i*n_nom.seg_T;
params_nom.Qb = params_all.Qb(:,idx_i);
params_nom.Tamb = params.Tamb(:,idx_i_T);
params_nom.Ts = params.Ts(:,idx_i_T);

vnom{1} = nominal_flow(G,n_nom,v,e,params,params_nom,1);
if vnom{1}.vld
    params_nom.mdot_e = vnom{1}.mdot_e;
    params_nom.mI = vnom{1}.mI;
    params_nom.zeta_u = vnom{1}.zeta_u;
    params_nom.dP_e = vnom{1}.dP_e;
    params_nom.P_nd = vnom{1}.P_nd;
    vnom{1} = nominal_flow(G,n_nom,v,e,params,params_nom,0);
end

%% Rest
for i = 157:n_cycle
    disp(i)
    idx_i = (i-1)*n_nom.seg+1:i*n_nom.seg;
    idx_i_T = (i-1)*n_nom.seg_T+1:i*n_nom.seg_T;
    params_nom.T0 = vnom{i-1}.T(:,end);
    params_nom.Qb = params_all.Qb(:,idx_i);
    params_nom.Tamb = params.Tamb(:,idx_i_T);
    params_nom.Ts = params.Ts(:,idx_i_T);
    vnom{i} = nominal_flow_slack(G,n_nom,v,e,params,params_nom,0);
    disp(sum(vnom{i}.Qp_slack.^2,'all'))
    if ~vnom{i}.vld
        disp(i);
    end
end

%% Store Nominal in single vector

vnom_all = vnom{1};
for i = 2:252
    vnom_all.dP_e = [vnom_all.dP_e vnom{i}.dP_e];
    vnom_all.P_nd = [vnom_all.P_nd vnom{i}.P_nd];
    vnom_all.mdot_e = [vnom_all.mdot_e vnom{i}.mdot_e];
    vnom_all.T = [vnom_all.T vnom{i}.T];
    vnom_all.Qp = [vnom_all.Qp vnom{i}.Qp];
    vnom_all.mI = [vnom_all.mI vnom{i}.mI];
    vnom_all.cost = vnom_all.cost+vnom{i}.cost;
end

%% Cost Calculations
e_fin = inedges(G,v.term_idx);
e_fin_nu = find(any(e.nu_idx==e_fin));
cost_T_nom = 0;
cost_T_sim = 0;
idx = 0;
for i = 1:size(vsim_all.T,2)
    if mod((i-1)*params.dt_T,params.dt)==0      % when mdot can change
        idx = idx+1;
    end
        cost_T_nom = cost_T_nom+(vnom_all.mdot_e(e_fin,idx).*vnom_all.T(e_fin_nu,i))'*(params.w_T*eye(numel(e_fin)))*(vnom_all.mdot_e(e_fin,idx).*vnom_all.T(e_fin_nu,i));
        cost_T_sim = cost_T_sim+(vsim_all.mdot_e(e_fin,idx).*vsim_all.T(e_fin_nu,i))'*(params.w_T*eye(numel(e_fin)))*(vsim_all.mdot_e(e_fin,idx).*vsim_all.T(e_fin_nu,i));
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

cd 'C:\Users\akb42\OneDrive - The Ohio State University\DistrictHeatingNetwork\Project Codes\Time Flexible Optimization\figures'
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

fig_graph(G,n,e,v)
fig_part(G,Gred,n,e,sG)
%fig_graphrd(Gred,n)
fig_demand(params,params_all,e,n,tp,1003)
fig_intQ(vsim_all,params,params_all,e,n)
fig_mdot_u(vsim_all,params,e,n)
fig_cost(vLL{12*6},vec_dP,n)
fig_mI(vsim_all, vnom,params,e)
fig_tempprofiles(tp,temp_prof)

set(groot,'defaultAxesTickLabelInterpreter','tex');  
set(groot,'defaulttextinterpreter','tex');
set(groot,'defaultLegendInterpreter','tex');
cd 'C:\Users\akb42\OneDrive - The Ohio State University\DistrictHeatingNetwork\Project Codes\Time Flexible Optimization'