function [cfin, Pfin, idxfin] = opt_high_par(v_list,n,vec_dP,vec_par,pipes)
%OPT_HIGH_PAR  One-line summary of what the function does.
%
%   [cfin, Pfin, idxfin] = OPT_HIGH_PAR(v_list,n,vec_dP,vec_par,pipes)
%
%   DESCRIPTION: Problem specific calculation of high level cost associated
%   with each pressure drop pairs. Finds the minimum overall cost where all
%   subsystems have a valid solution. Outputs the cost, final pressure
%   drop, and problem index.
%
%   INPUTS:
%       v_list  - Cell of solutions corresponding to vec_par
%       n       - Structure of sizes.
%       vec_dP  - Vector of pressure losses.
%       vec_par - Vector of subgraph pressure loss pairs considered.
%       pipes   - Matrix of pipe parameters.
%
%   OUTPUTS:
%       cfin    - Cost of optimal solution
%       Pfin    - Vector of pressure drops for optimal solution.
%       idxfin  - Vector of indices of optimal solution.

%% Parse information into matrices by subgraph

vld = zeros(n.sg, n.dP);
cost = zeros(n.sg, n.dP);
v = cell(n.sg, 1);
for idx_sg = 1:n.sg
    idx_par = vec_par(1,:)==idx_sg;
    v{idx_sg,:} = v_list(1,idx_par)';
    vld(idx_sg,:) = cellfun(@(x)x.valid, v_list(1,idx_par));
    cost(idx_sg,:) = cellfun(@(x)x.cost, v_list(1,idx_par))';   
end

vld(isnan(vld))=0;
vld = logical(vld);
cost(~vld)=NaN;

%% Find valid costs

% Preallocate
P_c = zeros(3,n.dP);
c = zeros(5,n.dP);
idx_c = zeros(2,n.dP);

idx_solution = 1;
% for any pressure in partition 1 with a valid solution
for dPg1_idx = find(vld(1,:))
    % if partition 2 has a valid solution
    if vld(2,dPg1_idx)    
        % calculate pressure loss in feeding+return edge
        dP_feeding = 2*pipes(1,3)*(v{1}{dPg1_idx}.mI(1)+v{2}{dPg1_idx}.mI(1))^2;
        % Pressure loss value closest to pressure loss p1+p2
        [~,dPg2_idx] = min(abs(vec_dP-(vec_dP(dPg1_idx)+dP_feeding)));
        % if all other partitions have solutions for this pressure
        if vld(3,dPg2_idx) && vld(4,dPg2_idx) && vld(5,dPg2_idx)
            % store pressures of 2 subgraphs
            P_c(1,idx_solution) = vec_dP(dPg1_idx);
            P_c(2,idx_solution) = vec_dP(dPg2_idx);
            P_c(3,idx_solution) = vec_dP(dPg1_idx)+dP_feeding;
            % store costs
            c(1,idx_solution) = cost(1,dPg1_idx);
            c(2,idx_solution) = cost(2,dPg1_idx);
            c(3,idx_solution) = cost(3,dPg2_idx);
            c(4,idx_solution) = cost(4,dPg2_idx);
            c(5,idx_solution) = cost(5,dPg2_idx);
            % store results index
            idx_c(1,idx_solution) = dPg1_idx;
            idx_c(2,idx_solution) = dPg2_idx;
            idx_solution = idx_solution+1;
        end
    end
end

% Trim unused space
P_c = P_c(:,1:idx_solution-1);
c = c(:,1:idx_solution-1);
idx_c = idx_c(:,1:idx_solution-1);

%% Find true minimum
if idx_solution>1
    [cfin, idx_solution] = min(sum(c));
    Pfin = P_c(:,idx_solution);
    idxfin = idx_c(:,idx_solution);
else
    error('No valid high level solutions found')
end

end