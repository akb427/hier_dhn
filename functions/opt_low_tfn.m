function [M,e_fin] = opt_low_tfn(G,n,nd,e,params,sparams)
%OPTIMIZE_FLOW Calculate mass flow rate breakdown given a set dP
%   G: graph of network
%   params: parameters of network
%   n: structure of sizes
%   given mI, dP 


%% Setup Problem 
import casadi.*
opti_flow = casadi.Opti();

%% Flow variables

% Mass flow
mdot_e = opti_flow.variable(n.e, n.seg);
opti_flow.subject_to(mdot_e(:)>=0);

mI = opti_flow.variable(1, n.seg);
opti_flow.subject_to(mI(:)>=0);

mdot_n = MX.zeros(n.v-1,n.seg);
mdot_n(nd.root_idx,:) = mI;

% Zeta
zeta_u = opti_flow.variable(n.u, n.seg);
opti_flow.subject_to(zeta_u(:)>=0);

% Pressures
P_n = opti_flow.variable(n.v, n.seg);
dP_e = opti_flow.variable(n.e, n.seg);
opti_flow.subject_to(dP_e(:)>=0);
opti_flow.subject_to(P_n(:)>=0);

dP = opti_flow.parameter(1);

% Remove overconstraint for mdot
Ired = sparams.I;
Ired(nd.term_idx,:) = [];

%% Constraints

for i = 1:n.seg
    % KCL
    opti_flow.subject_to(Ired*mdot_e(:,i)==mdot_n(:,i));
    
    % KVL
    opti_flow.subject_to(dP_e(:,i)==sparams.I'*P_n(:,i));
    % Set reference pressue
    opti_flow.subject_to(P_n(nd.term_idx,i)==0);
    % Enforce required pressure drop
    if (i-1)*params.dt<=params.tf_sim
        opti_flow.subject_to(P_n(nd.root_idx,i)==dP);
    else
        opti_flow.subject_to(P_n(nd.root_idx,i)<=params.dP_max);
    end
    % Edge equations
    opti_flow.subject_to(dP_e(e.nu_idx,i) == sparams.pipes(e.nu_idx,3).*mdot_e(e.nu_idx,i).^2);
    opti_flow.subject_to(dP_e(e.u_idx,i) == zeta_u(:,i).*mdot_e(e.u_idx,i).^2);
end

%% Temperature Variables

T0 = opti_flow.parameter(n.nu,1);
Ts = opti_flow.parameter(1,n.seg_T);
Tamb = opti_flow.parameter(1,n.seg_T);
intQ0 = opti_flow.parameter(n.u,1);
Qb = opti_flow.parameter(n.u,n.seg);
Cap_u = opti_flow.parameter(n.u,n.seg_T);
Cap_l = opti_flow.parameter(n.u,n.seg_T);

T = MX(n.nu,n.seg_T);
Qp = MX(n.u, n.seg_T);
intQ = MX(n.u,n.seg_T);

idx = 1;
i = 1;
T(:,i) = T0;
Qp(:,i) = params.cp*mdot_e(e.u_idx,idx).*(T(e.ui_idxnu,i)-params.TsetR)/10^3; % kW
intQ(:,i) = intQ0+((Qp(:,i)-Qb(:,idx))*params.dt_T)/10^3; % MJ
[A,B,E] = graph2ss(G,params,sparams,n,nd,e,mdot_e(:,idx),0);

e_fin = find(e.nu_idx==inedges(G,nd.term_idx));
cost_T = MX(1,1);
cost_Q = MX(1,1);

cost_T = cost_T+(mI(1,idx)*T(e_fin,i))*params.w_T*(mI(1,idx)*T(e_fin,i));
cost_Q = cost_Q+intQ(:,i)'*sparams.w_flex*intQ(:,i);

%% Calculate Dynamics

for i = 2:n.seg_T
    if mod((i-1)*params.dt_T,params.dt)==0
        idx = idx+1;
        % Update SS
        [A,B,E] = graph2ss(G,params,sparams,n,nd,e,mdot_e(:,idx),0);
        T(:,i)=(MX.eye(size(A,1))+A*params.dt_T)*T(:,i-1)+params.dt_T*[B,E]*[Ts(1,i-1);params.TsetR;Tamb(1,i-1)];
        % Calculate provided heat
        Qp(:,i) = params.cp*mdot_e(e.u_idx,idx).*(T(e.ui_idxnu,i)-params.TsetR)/10^3;   % kW
        % Ensure envelope is never exceeded
        intQ(:,i) = intQ(:,i-1)+((Qp(:,i)-Qb(:,idx))*params.dt_T)/10^3; % MJ
        opti_flow.subject_to(Cap_l(:,i)<=intQ(:,i)<=Cap_u(:,i));
        % Update Cost
        cost_T = cost_T+(mI(1,idx)*T(e_fin,i))*params.w_T*(mI(1,idx)*T(e_fin,i));
        cost_Q = cost_Q+intQ(:,i)'*sparams.w_flex*intQ(:,i);
    else
        T(:,i)=(MX.eye(size(A,1))+A*params.dt_T)*T(:,i-1)+params.dt_T*[B,E]*[Ts(1,i-1);params.TsetR;Tamb(1,i-1)];
        % Calculate provided heat
        Qp(:,i) = params.cp*mdot_e(e.u_idx,idx).*(T(e.ui_idxnu,i)-params.TsetR)/10^3;   % kW
        % Ensure envelope is never exceeded
        intQ(:,i) = intQ(:,i-1)+((Qp(:,i)-Qb(:,idx))*params.dt_T)/10^3; % MJ
        opti_flow.subject_to(Cap_l(:,i)<=intQ(:,i)<=Cap_u(:,i));
        % Update Cost
        cost_T = cost_T+(mI(1,idx)*T(e_fin,i))*params.w_T*(mI(1,idx)*T(e_fin,i));
        cost_Q = cost_Q+intQ(:,i)'*sparams.w_flex*intQ(:,i);
    end
end

opti_flow.subject_to(-30<T(:)<100);
cost = cost_T;%+cost_Q;

%% Solver
mdot_by = mdot_e(e.by_idx,:);

opti_flow.minimize(cost);

opti_flow.solver('ipopt',struct('print_time',0),struct('print_level',0,'tol', 1e-2,'max_iter',900))

M = opti_flow.to_function('M',{dP,T0,intQ0,Qb,Tamb,Ts,Cap_u,Cap_l,mdot_e,zeta_u,dP_e,mI,P_n,opti_flow.lam_g},...
    {dP_e,P_n,mdot_e,zeta_u,T,Qp,intQ,mI,opti_flow.lam_g,mdot_by,cost_T,cost_Q,cost},...
    {'dP','T0','intQ0','Qb','Tamb','Ts','Cap_u','Cap_l','i_mdot_e','i_zeta_u','i_dP_e','i_mI','i_P_n','i_lam_g'},...
    {'dP_e','P_n','mdot_e','zeta_u','T','Qp','intQ','mI','lam_g','mdot_by','cost_T','cost_Q','cost'});

end







