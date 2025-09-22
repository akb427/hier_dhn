function [Mnom] = nominal_flow_tfn(G,n,v,e,params,isSlack)
%NOMINAL_FLOW_TFN  Create function to solve for network with nominal demands.
%
%   [Mnom] = NOMINAL_FLOW_TFN(G,n,nd,e,params,isSlack)
%
%   DESCRIPTION: Creates the CasADi function that optimizes the network
%   behavior subject to set heat demands in the buildings. The resulting 
%   function takes the supply temperature, initial conditions, ambient 
%   conditions building conditions, and initial guesses. It returns 
%   pressure variables, flow variables, friction coefficients, 
%   temperatures, various costs, and the gradients. The slack flag allows 
%   slack in the equality constraint for meeting heat demand and penalizes 
%   it in the cost function. 
%   
%   INPUTS:
%       G   - Graph of network.
%       n   - Structure of sizes.
%       v   - Structure of node information.
%       e   - Structure of edge information.
%       params  - Structure of network-wide parameters.
%       isSlack - Binary indicating slack in demand equality constraint.
%
%   OUTPUTS:
%       Mnom - CasADI function.
%
%   DEPENDENCIES: graph2ss

%#ok<*CHAIN>   % This is okay in CASADI
%#ok<*FNDSB>   % This is nessecary in CASADI

%% Setup Problem 
import casadi.*
opti_flow = casadi.Opti();

%% Flow variables

% Mass flow
mdot_e = opti_flow.variable(n.e, n.seg);
opti_flow.subject_to(mdot_e(:)>=0);

mI = opti_flow.variable(1,n.seg);
opti_flow.subject_to(mI(:)>0);

mdot_n = MX.zeros(n.v-1,n.seg);
mdot_n(v.root_idx,:) = mI;

% Zeta
zeta_u = opti_flow.variable(n.u, n.seg);
opti_flow.subject_to(zeta_u(:)>=0);

% Pressures
P_n = opti_flow.variable(n.v, n.seg);
opti_flow.subject_to(P_n(:)>=0);

dP_e = opti_flow.variable(n.e, n.seg);
opti_flow.subject_to(dP_e(:)>=0);

% Remove overconstraint for mdot
Ired = params.Inc;
Ired(v.term_idx,:) = [];

%% Constraints
for i = 1:n.seg
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
Ts = opti_flow.parameter(1,n.seg_T);
Tamb = opti_flow.parameter(1,n.seg_T);
Qb = opti_flow.parameter(n.u,n.seg);

T = MX(n.nu,n.seg_T);
Qp = MX(n.u, n.seg_T);
if isSlack
    Qp_slack = opti_flow.variable(n.u,n.seg);
end

i = 1;
idx = 1;
[A,B,E] = graph2ss(G,params,params,n,v,e,mdot_e(:,idx),0);
T(:,i)=(MX.eye(size(A,1))+A*params.dt_T)*T0+params.dt_T*[B,E]*[Ts(1,i-1);params.TsetR;Tamb(1,i-1)];
Qp(:,i) = params.cp*mdot_e(e.u_idx,idx).*(T(e.ui_idxnu,i)-params.TsetR)/1000;
if isSlack
    opti_flow.subject_to(Qp(:,i)==Qb(:,idx)+Qp_slack(:,idx));
end

e_fin = inedges(G,v.term_idx);
e_fin_nu = find(any(e.nu_idx==e_fin));
cost_T = MX(1,1);
cost_T = cost_T+(mdot_e(e_fin,idx).*T(e_fin_nu,i))'*(params.w_T*eye(numel(e_fin)))*(mdot_e(e_fin,idx).*T(e_fin_nu,i));

%% Calculate Dynamics

for i = 2:n.seg_T
    if mod((i-1)*params.dt_T,params.dt)==0      % when mdot can change
        idx = idx+1;
        % Update SS
        [A,B,E] = graph2ss(G,params,params,n,v,e,mdot_e(:,idx),0);
        T(:,i)=(MX.eye(size(A,1))+A*params.dt_T)*T(:,i-1)+params.dt_T*[B,E]*[Ts(1,i-1);params.TsetR;Tamb(1,i-1)];
        % Calculate provided heat
        Qp(:,i) = params.cp*mdot_e(e.u_idx,idx).*(T(e.ui_idxnu,i)-params.TsetR)/1000;
        % Ensure provided heat is equal to demanded heat 
        if isSlack
            opti_flow.subject_to(Qp(:,i)==Qb(:,idx)+Qp_slack(:,idx));
        else
            opti_flow.subject_to(Qp(:,i)==params.Qb(:,idx));
        end
        cost_T = cost_T+(mdot_e(e_fin,idx).*T(e_fin_nu,i))'*(params.w_T*eye(numel(e_fin)))*(mdot_e(e_fin,idx).*T(e_fin_nu,i));
    else
        T(:,i)=(MX.eye(size(A,1))+A*params.dt_T)*T(:,i-1)+params.dt_T*[B,E]*[Ts(1,i-1);params.TsetR;Tamb(1,i-1)];
        % Calculate provided heat
        Qp(:,i) = params.cp*mdot_e(e.u_idx,idx).*(T(e.ui_idxnu,i)-params.TsetR)/1000;
        cost_T = cost_T+(mdot_e(e_fin,idx).*T(e_fin_nu,i))'*(params.w_T*eye(numel(e_fin)))*(mdot_e(e_fin,idx).*T(e_fin_nu,i));
    end
end

opti_flow.subject_to(-30<T(:)<100);

%% Solve zeta for minimum
mdot_by = mdot_e(e.by_idx,:);

if isSlack
    opti_flow.minimize(cost_T+10^8*sum(Qp_slack(:).^2));
else
    opti_flow.minimize(sum(mdot_by(:)));
end

%% Create Function

opti_flow.solver('ipopt',struct('print_time',0),struct('print_level',0,'tol', 1e-2,'max_iter', 100000))

inpt = {T0,Qb,Tamb,Ts,mdot_e,zeta_u,dP_e,mI,P_n,opti_flow.lam_g};
outpt = {dP_e,P_n,mdot_e,zeta_u,T,Qp,mI,opti_flow.lam_g,mdot_by,cost_T};
inpt_name = {'T0','Qb','Tamb','Ts','i_mdot_e','i_zeta_u','i_dP_e','i_mI','i_P_n','i_lam_g'};
outpt_name = {'dP_e','P_n','mdot_e','zeta_u','T','Qp','intQ','mI','lam_g','mdot_by','cost'};

% Add slack output
if isSlack
    outpt{end+1} = Qp_slack;
    outpt_name{end+1} = "Qp_slack";
end

% Create function
Mnom = opti_flow.to_function('M',inpt,outpt,inpt_name,outpt_name);

end





