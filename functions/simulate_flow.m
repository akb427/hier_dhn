function [sol] = simulate_flow(G,n,v,e,params,params_sim,changeIG)
%SIMULATE_FLOW Calculate mass flow rate and pressure losses given zeta.
%
%   [sol] = simulate_flow(G,n,v,e,params,params_sim,changeIG)
%
%   DESCRIPTION: Simulates the network flow given the valve positions,
%   plant conditions, and ambient variables. Resolves pressure drop and
%   mass flow in all network edges. It outputs the simulation results in a
%   structure. Has two different intial guess magnitudes for if the problem
%   does not solve. Uses CasADi.
%
%   INPUTS:
%       G       - Graph of network.
%       n       - Structure of sizes.
%       v       - Structure of node infomration.
%       e       - Structure of edge information.
%       params  - Structure of network components parameter.
%       params_sim  - Structure of set variables for the simulation.
%       changeIG    - Binary flag indicating which intial guess to use.
%
%   OUTPUTS:
%       sol - Strucuture of solution variables.
%
%   DEPENDENCIES: graph2ss

%% Setup Problem 
import casadi.*
opti_flow = casadi.Opti();

%% Flow variables

% Parameters
mI = opti_flow.parameter(1,n.seg_sim);
zeta_u = opti_flow.parameter(n.u,n.seg_sim);

% Mass flow
mdot_e = opti_flow.variable(n.e, n.seg_sim);
opti_flow.subject_to(mdot_e(:)>=0);

mdot_n = MX.zeros(n.v-1,n.seg_sim);
mdot_n(v.root_idx,:) = mI;

% Pressures
P_n = opti_flow.variable(n.v, n.seg_sim);
opti_flow.subject_to(P_n(:)>=0);

dP_e = opti_flow.variable(n.e, n.seg_sim);
opti_flow.subject_to(dP_e(:)>=0);

% Remove overconstraint for mdot
Ired = params.Inc;
Ired(v.term_idx,:) = [];

%% Constraints
for i = 1:n.seg_sim
    % KCL
    opti_flow.subject_to(Ired*mdot_e(:,i)==mdot_n(:,i));
    
    % KVL
    opti_flow.subject_to(dP_e(:,i)==params.Inc'*P_n(:,i));
    
    % Set reference pressue
    opti_flow.subject_to(P_n(v.term_idx,i)==0);
    
    % Edge equations
    opti_flow.subject_to(dP_e(e.nu_idx,i) == params.pipes(e.nu_idx,3).*mdot_e(e.nu_idx,i).^2);
    opti_flow.subject_to(dP_e(e.u_idx,i) == zeta_u(:,i).*mdot_e(e.u_idx,i).^2);
end
%% Temperature Dynamics

T0 = opti_flow.parameter(n.nu,1);
Ts = opti_flow.parameter(1,n.seg_T_sim);
Tamb = opti_flow.parameter(1,n.seg_T_sim);
intQ0 = opti_flow.parameter(n.u,1);
Qb = opti_flow.parameter(n.u,n.seg_sim);

T = MX(n.nu,n.seg_T_sim);
Qp = MX(n.u, n.seg_T_sim);
intQ = MX(n.u,n.seg_T_sim);

idx = 1;
i = 1;
T(:,i) = T0;
Qp(:,i) = params.cp*mdot_e(e.u_idx,idx).*(T(e.ui_idxnu,i)-params.TsetR)/10^3; % kW
intQ(:,i) = intQ0+(Qp(:,i)-Qb(:,idx))*params.dt_T/10^3; %MJ
[A,B,E] = graph2ss(G,params,params,n,v,e,mdot_e(:,idx),0);

for i = 2:n.seg_T_sim
    if mod((i-1)*params.dt_T,params.dt)==0      % when mdot can change
        idx = idx+1;
        % Update SS
        [A,B,E] = graph2ss(G,params,params,n,v,e,mdot_e(:,idx),0);
        T(:,i)=(MX.eye(size(A,1))+A*params.dt_T)*T(:,i-1)+params.dt_T*[B,E]*[Ts(1,i-1);params.TsetR;Tamb(1,i-1)];
        % Calculate provided heat
        Qp(:,i) = params.cp*mdot_e(e.u_idx,idx).*(T(e.ui_idxnu,i)-params.TsetR)/10^3; %kW
        intQ(:,i) = intQ(:,i-1)+((Qp(:,i)-Qb(:,idx))*params.dt_T)/10^3; %MJ
    else
        T(:,i)=(MX.eye(size(A,1))+A*params.dt_T)*T(:,i-1)+params.dt_T*[B,E]*[Ts(1,i-1);params.TsetR;Tamb(1,i-1)];
        Qp(:,i) = params.cp*mdot_e(e.u_idx,idx).*(T(e.ui_idxnu,i)-params.TsetR)/10^3;   %kW
        intQ(:,i) = intQ(:,i-1)+((Qp(:,i)-Qb(:,idx))*params.dt_T)/10^3; %MJ
    end
end

%% Set initial conditions and parameters

% Parameter values
opti_flow.set_value(mI, params_sim.mI);
opti_flow.set_value(zeta_u, params_sim.zeta_u);
opti_flow.set_value(Qb, params_sim.Qb);
opti_flow.set_value(Tamb,params_sim.Tamb);
opti_flow.set_value(Ts, params_sim.Ts);

opti_flow.set_value(T0, params_sim.T0);
opti_flow.set_value(intQ0, params_sim.intQ);

% Initial conditions
if changeIG
    opti_flow.set_initial(mdot_e,10*ones(n.e,n.seg_sim));
else
    opti_flow.set_initial(mdot_e,ones(n.e,n.seg_sim));
end
opti_flow.set_initial(dP_e,params_sim.dP_e(:,1:n.seg_sim));


%% Solve optimization problem

mdot_by = mdot_e(e.by_idx,:);
opti_flow.minimize(1);

opti_flow.solver('ipopt',struct('print_time',0),struct('print_level',0,'tol', 1e-2))

try
    sol = opti_flow.solve;
    % output values
    sol.dP_e = sol.value(dP_e);
    sol.P_nd = sol.value(P_n);
    sol.mdot_e = sol.value(mdot_e);
    sol.T = sol.value(T);
    sol.Qp = sol.value(Qp);
    sol.intQ = sol.value(intQ);
    sol.mI = params_sim.mI;
    sol.cost = sum(sol.value(mdot_by),'all');
    sol.valid = 1;
catch
    opti_flow.debug.show_infeasibilities(10^-2)
    % output values if solver doesnt converge
    sol.dP_e = opti_flow.debug.value(dP_e);
    sol.P_nd = opti_flow.debug.value(P_n);
    sol.mdot_e = opti_flow.debug.value(mdot_e);
    sol.T = opti_flow.debug.value(T);
    sol.Qp = opti_flow.debug.value(Qp);
    sol.mI = params_sim.mI;
    sol.cost = sum(opti_flow.debug.value(mdot_by),'all');
    sol.valid = 0;
end





