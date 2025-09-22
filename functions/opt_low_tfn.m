function [M,e_fin] = opt_low_tfn(G,n,v,e,params,sparams,isSlack)
%OPT_LOW_TFN  Generate CasADi function to calculate losses given a set
%pressure drop. 
%
%   [M,e_fin] = OPT_LOW_TFN(G,n,nd,e,params,sparams,isSlack)
%
%   DESCRIPTION: Creates the CasADi function that optimizes the network
%   behavior subject to a set pressure drop in the first time step. The
%   resulting function takes the set pressure drop, supply temperature,
%   initial conditions, ambient conditions and allowable temperature
%   deviations, and initial guesses. It returns pressure varaiables, flow
%   variables, friction coefficients, temperatures, building SOE, various
%   costs, and the gradients. The slack flag allows slack in the building
%   SOE and penalizes it in the cost function. 
%
%   INPUTS:
%       G   - Graph of network.
%       n   - Structure of sizes.
%       v   - Structure of node information.
%       e   - Structure of edge information.
%       params  - Structure of network-wide parameters.
%       sparams - Structure of parameters in current network.
%       isSlack - Binary indicating slack in SOE limits.
%
%   OUTPUTS:
%       M   - CasADi function
%       e_fin   - Vector of edges leading to terminal nodes.

%#ok<*CHAIN>   % This is okay in CASADI
%#ok<*FNDSB>   % This is nessecary in CASADI

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
mdot_n(v.root_idx,:) = mI;

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
Ired(v.term_idx,:) = [];

%% Constraints

for i = 1:n.seg
    % KCL
    opti_flow.subject_to(Ired*mdot_e(:,i)==mdot_n(:,i));
    
    % KVL
    opti_flow.subject_to(dP_e(:,i)==sparams.I'*P_n(:,i));
    % Set reference pressue
    opti_flow.subject_to(P_n(v.term_idx,i)==0);
    % Enforce required pressure drop
    if (i-1)*params.dt<=params.tf_sim
        opti_flow.subject_to(P_n(v.root_idx,i)==dP);
    else
        opti_flow.subject_to(P_n(v.root_idx,i)<=params.dP_max);
    end
    % Edge equations
    opti_flow.subject_to(dP_e(e.nu_idx,i) == sparams.pipes(e.nu_idx,3).*mdot_e(e.nu_idx,i).^2);
    opti_flow.subject_to(dP_e(e.u_idx,i) == zeta_u(:,i).*mdot_e(e.u_idx,i).^2);
end

%% Temperature Variables

% Set parameters
T0 = opti_flow.parameter(n.nu,1);
Ts = opti_flow.parameter(1,n.seg_T);
Tamb = opti_flow.parameter(1,n.seg_T);
intQ0 = opti_flow.parameter(n.u,1);
Qb = opti_flow.parameter(n.u,n.seg);
Cap_u = opti_flow.parameter(n.u,n.seg_T);
Cap_l = opti_flow.parameter(n.u,n.seg_T);

% Derived variables
T = MX(n.nu,n.seg_T);
Qp = MX(n.u, n.seg_T);
intQ = MX(n.u,n.seg_T);
cost_T = MX(1,1);
cost_Q = MX(1,1);

% Optional slack in heat capacity
if isSlack
    slack_Cap = opti_flow.variable(n.u,n.seg);
    opti_flow.subject_to(slack_Cap(:)>=0);
end

% Initial conditions
idx = 1;
i = 1;
T(:,i) = T0;
Qp(:,i) = params.cp*mdot_e(e.u_idx,idx).*(T(e.ui_idxnu,i)-params.TsetR)/10^3; % kW
intQ(:,i) = intQ0+((Qp(:,i)-Qb(:,idx))*params.dt_T)/10^3; % MJ
[A,B,E] = graph2ss(G,params,sparams,n,v,e,mdot_e(:,idx),0);

% Edges leaving subsystem
e_fin = find(e.nu_idx==inedges(G,v.term_idx));

cost_T = cost_T+(mI(1,idx)*T(e_fin,i))*params.w_T*(mI(1,idx)*T(e_fin,i));
cost_Q = cost_Q+intQ(:,i)'*sparams.w_flex*intQ(:,i);

%% Calculate Dynamics

% Discretized over dt_T, controllable at dt
for i = 2:n.seg_T
    % If step allows control change
    if mod((i-1)*params.dt_T,params.dt)==0
        idx = idx+1;
        % Update SS
        [A,B,E] = graph2ss(G,params,sparams,n,v,e,mdot_e(:,idx),0);
        T(:,i)=(MX.eye(size(A,1))+A*params.dt_T)*T(:,i-1)+params.dt_T*[B,E]*[Ts(1,i-1);params.TsetR;Tamb(1,i-1)];
        % Calculate provided heat
        Qp(:,i) = params.cp*mdot_e(e.u_idx,idx).*(T(e.ui_idxnu,i)-params.TsetR)/10^3;   % kW
        % Ensure envelope is never exceeded
        intQ(:,i) = intQ(:,i-1)+((Qp(:,i)-Qb(:,idx))*params.dt_T)/10^3; % MJ
        if isSlack
            opti_flow.subject_to((Cap_l(:,i)-slack_Cap(:,idx))<=intQ(:,i)<=Cap_u(:,i)+slack_Cap(:,idx));
        else
            opti_flow.subject_to(Cap_l(:,i)<=intQ(:,i)<=Cap_u(:,i));
        end
        % Update Cost
        cost_T = cost_T+(mI(1,idx)*T(e_fin,i))*params.w_T*(mI(1,idx)*T(e_fin,i));
        cost_Q = cost_Q+intQ(:,i)'*sparams.w_flex*intQ(:,i);
    else % If valve positions cant change
        % State space remains the same
        T(:,i)=(MX.eye(size(A,1))+A*params.dt_T)*T(:,i-1)+params.dt_T*[B,E]*[Ts(1,i-1);params.TsetR;Tamb(1,i-1)];
        % Calculate provided heat
        Qp(:,i) = params.cp*mdot_e(e.u_idx,idx).*(T(e.ui_idxnu,i)-params.TsetR)/10^3;   % kW
        % Ensure envelope is never exceeded
        intQ(:,i) = intQ(:,i-1)+((Qp(:,i)-Qb(:,idx))*params.dt_T)/10^3; % MJ
        if isSlack
            opti_flow.subject_to((Cap_l(:,i)-slack_Cap(:,idx))<=intQ(:,i)<=(Cap_u(:,i)+slack_Cap(:,idx)));
        else
            opti_flow.subject_to(Cap_l(:,i)<=intQ(:,i)<=Cap_u(:,i));
        end
        % Update Cost
        cost_T = cost_T+(mI(1,idx)*T(e_fin,i))*params.w_T*(mI(1,idx)*T(e_fin,i));
        cost_Q = cost_Q+intQ(:,i)'*sparams.w_flex*intQ(:,i);
    end
end

% Temperature limits
opti_flow.subject_to(-30<T(:)<100);

%% Solver

% Cost variables
cost = cost_T;
mdot_by = mdot_e(e.by_idx,:);

% Include slack if needed
if isSlack
    opti_flow.minimize(cost+sum(slack_Cap(:)));
else
    opti_flow.minimize(cost);
end

% Solver setup
if isSlack
    n_iter = 1500;
else
    n_iter = 900;
end
opti_flow.solver('ipopt',struct('print_time',0),struct('print_level',0,'tol', 1e-2,'max_iter',n_iter))

% Create function
inpt = {dP,T0,intQ0,Qb,Tamb,Ts,Cap_u,Cap_l,mdot_e,zeta_u,dP_e,mI,P_n,opti_flow.lam_g};
outpt = {dP_e,P_n,mdot_e,zeta_u,T,Qp,intQ,mI,opti_flow.lam_g,mdot_by,cost_T,cost_Q,cost};
inpt_name = {'dP','T0','intQ0','Qb','Tamb','Ts','Cap_u','Cap_l','i_mdot_e','i_zeta_u','i_dP_e','i_mI','i_P_n','i_lam_g'};
outpt_name = {'dP_e','P_n','mdot_e','zeta_u','T','Qp','intQ','mI','lam_g','mdot_by','cost_T','cost_Q','cost'};

% Add slack variables
if isSlack
    outpt{end+1} = slack_Cap;
    outpt_name{end+1} = 'Cap_slack';
end

% Create function
M = opti_flow.to_function('M',inpt,outpt,inpt_name,outpt_name);


end







