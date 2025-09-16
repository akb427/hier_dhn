%% Load variables
clc, clear, close all
load('opti_variables_notime')

%% Parallel Pool Setup

%parpool(19);

%% Load Initial Guesses
load('init_v');
vLL_ig = vLL_ig1;

%% Set values for i=1
params_HL.T0 = params.Ts(1)*ones(n.nu,1);
params_HL.T0(e.uo_idxnu,1) = params.TsetR;
params_HL.intQ = zeros(n.u,1);
flag_rmlam = false(1,n.par);
flag_slack = false(1,n.par);
flag_slack_new = false(1,n.par);
M_i = M;

%% Run Optimization
%load('start_for136');
tic
for i = 38:n_step
    disp(i)
    % Timing indices
    idx_i = i:(i-1+n.seg);
    idx_i_sim = i:(i-1+n.seg_sim);
    idx_i_T = (n.seg_T_sim*(i-1)+1):(n.seg_T_sim*(i-1)+n.seg_T);
    idx_i_T_sim = (n.seg_T_sim*(i-1)+1):(n.seg_T_sim*(i-1)+n.seg_T_sim);
    % Deal variables for low level controller
    params_i = cell(1,n.par);
    for j = 1:n.sg
        idx_sg = vec_par(1,:) == j;
        st.Tamb = params.Tamb(:,idx_i_T);
        st.Ts = params.Ts(:,idx_i_T);
        if i==1
            % Initial conditions
            st.T0 = params.Ts(1)*ones(sn{j}.nu,1);
            st.T0(se{j}.uo_idxnu,1) = params.TsetR;
            st.intQ0 = zeros(sn{j}.u,1);
        else
            % Initial conditions
            st.T0 = init.T(se{j}.nu,:);
            st.intQ0 = init.intQ(se{j}.u,:);
        end
        st.Qb = sparams{j}.Qb(:,idx_i);
        st.Cap_u = sparams{j}.Cap_u(:,idx_i_T);
        st.Cap_l = sparams{j}.Cap_l(:,idx_i_T);
        params_i(1,idx_sg) = repmat({st},1,n.dP);
    end
    % Low Level Controller
    %ticBytes(gcp)
    parfor idx = 1:n.par
        tic
        % Initial guesses
        params_idx = params_i{1,idx};
        params_idx.i_mdot_e = vLL_ig{1,idx}.mdot_e;
        params_idx.i_zeta_u = vLL_ig{1,idx}.zeta_u;
        params_idx.i_dP_e = vLL_ig{1,idx}.dP_e;
        params_idx.i_mI = vLL_ig{1,idx}.mI;
        params_idx.i_P_n = vLL_ig{1,idx}.P_n;
        params_idx.i_lam_g = vLL_ig{1,idx}.lam_g;
        params_idx.dP = vec_par(2,idx);
        if flag_rmlam(idx)
            params_idx = rmfield(params_idx,'i_lam_g');
        end
        % Call function
        vt = M_i{idx}.call(params_idx);
        vt.valid = M_i{idx}.stats.success;
        vLL_i{1,idx} = structfun(@full,vt,'UniformOutput',false);
        t_elapsed_i(1,idx) = toc;
    end
    %tocBytes(gcp)
    % High level
    [~, ~, idx_HL] = opt_high_par(vLL_i,n,vec_dP,vec_par,params.pipes);
    params_HL.zeta_u = zeros(n.u,n.seg);
    params_HL.mdot_e = zeros(n.e,n.seg);
    mI_partion = zeros(n.sg,n.seg);
    for j = 1:2
        for k = 1:sn{j}.u
            params_HL.zeta_u(se{j}.u(k)==e.u_idx,:) = vLL_i{1,(j-1)*n.dP+idx_HL(1)}.zeta_u(k,:);
        end
        params_HL.mdot_e(se{j}.e,:) = vLL_i{1,(j-1)*n.dP+idx_HL(1)}.mdot_e;
        params_HL.dP_e(se{j}.e,:) = vLL_i{1,(j-1)*n.dP+idx_HL(1)}.dP_e;
        mI_partion(j,:) = vLL_i{1,(j-1)*n.dP+idx_HL(1)}.mI;
    end
    for j = 3:5
        for k = 1:sn{j}.u
            params_HL.zeta_u(se{j}.u(k)==e.u_idx,:) = vLL_i{1,(j-1)*n.dP+idx_HL(2)}.zeta_u(k,:);
        end
        mI_partion(j,:) = vLL_i{1,(j-1)*n.dP+idx_HL(2)}.mI;
        params_HL.mdot_e(se{j}.e,:) = vLL_i{1,(j-1)*n.dP+idx_HL(2)}.mdot_e;
        params_HL.dP_e(se{j}.e,:) = vLL_i{1,(j-1)*n.dP+idx_HL(2)}.dP_e;
    end

    % Simulation
    params_HL.mI = sum(mI_partion);
    params_HL.mI = params_HL.mI(:,1:n.seg_sim);
    params_HL.zeta_u = params_HL.zeta_u(:,1:n.seg_sim);
    params_HL.Tamb = params.Tamb(:,idx_i_T_sim);
    params_HL.Ts = params.Ts(:,idx_i_T_sim);
    params_HL.Cap_u = params_all.Cap_u(:,idx_i_T_sim);
    params_HL.Cap_l = params_all.Cap_l(:,idx_i_T_sim);
    params_HL.Qb = params_all.Qb(:,idx_i_sim);
    
    vsim_i = simulate_flow(G,n,v,e,params,params_HL);
    if ~vsim_i.valid
        vsim_i = simulate_flow2(G,n,v,e,params,params_HL);
    end
    
    % initial values for next optimization
    init.T = zeros(n.e,1);
    init.T(e.nu,1) = vsim_i.T(:,end);
    init.intQ = zeros(n.e,1);
    init.intQ(e.u,1) = vsim_i.intQ(:,end);
    % inital values for next simulation
    params_HL.T0 = vsim_i.T(:,end);
    params_HL.intQ = vsim_i.intQ(:,end);
    % Save
    flname = pwd+"\vLL notime1\"+"vLL_nt_"+num2str(i);
    save(flname,'vLL_i','vsim_i','t_elapsed_i','idx_HL','init','params_HL');
    
    % Update Initial guess & create flags for next run
    vLL_ig =vLL_i;
    cost = cellfun(@(x)x.cost, vLL_i);
    vLL_ig(cost>10^4) = vLL_ig1(cost>10^4);
    for j = 1:n.sg
        flag_slack_new(vec_par(1,:)==j) = any(abs(init.intQ(se{j}.u,:))>sparams{j}.Cap_u(:,idx_i_T(1)));
    end
    flag_rmlam = flag_slack_new~=flag_slack;
    flag_rmlam(cost>10^4) = 1;
    flag_slack = flag_slack_new;
    M_i = M;
    M_i(flag_slack) = Ms(flag_slack);
end
toc