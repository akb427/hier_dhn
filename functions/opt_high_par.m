function [cfin, Pfin, idxfin] = opt_high_par(v_list,n,vec_dP,vec_par,pipes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Parse information into matrices

vld = zeros(n.sg, n.dP);
cost = zeros(n.sg, n.dP);
v = cell(n.sg, 1);
for sg_idx = 1:n.sg
    idx = vec_par(1,:)==sg_idx;
    v{sg_idx,:} = v_list(1,idx)';
    vld(sg_idx,:) = cellfun(@(x)x.valid, v_list(1,idx));
    cost(sg_idx,:) = cellfun(@(x)x.cost, v_list(1,idx))';   
end

vld(isnan(vld))=0;
vld = logical(vld);
cost(~vld)=NaN;

%% Find valid costs


solution_idx = 1;
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
            P_c(1,solution_idx) = vec_dP(dPg1_idx);
            P_c(2,solution_idx) = vec_dP(dPg2_idx);
            P_c(3,solution_idx) = vec_dP(dPg1_idx)+dP_feeding;
            % store costs
            c(1,solution_idx) = cost(1,dPg1_idx);
            c(2,solution_idx) = cost(2,dPg1_idx);
            c(3,solution_idx) = cost(3,dPg2_idx);
            c(4,solution_idx) = cost(4,dPg2_idx);
            c(5,solution_idx) = cost(5,dPg2_idx);
            % store results index
            idx_c(1,solution_idx) = dPg1_idx;
            idx_c(2,solution_idx) = dPg2_idx;
            solution_idx = solution_idx+1;
        end
    end
end


%% Find true minimum
if solution_idx>1
    [cfin, solution_idx] = min(sum(c));
    Pfin = P_c(:,solution_idx);
    idxfin = idx_c(:,solution_idx);
else
    cfin = NaN;
    Pfin = NaN;
    idxfin = NaN;
    error('No valid high level solutions found')
end

%% Sub functions
function out = get_value(x,value)
    if isa(x,'struct')
        out = x.(value);
    else
        out = NaN;
    end
end
end