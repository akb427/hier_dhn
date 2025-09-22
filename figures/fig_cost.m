function fig_cost(v,vec_dP,n)
%FIG_COST  Plots the low level costs vs pressure drops.
%
%   FIG_COST(v,vec_dP,n)
%
%   DESCRIPTION: Plots the low level costs vs pressure drops in a single
%   optimization time step. 
%   
%
%   INPUTS:
%       v       - Cell of low level results at a single timestep.
%       vec_dP  - Vector of pressure losses
%       n       - Structure of sizes.

%% Parse Data

vld = zeros(n.sg, n.dP);
cost = zeros(n.sg, n.dP);
for idx_sg = 1:n.sg
    isSolved = cellfun(@(x) isstruct(x) && isfield(x,'valid'), v{idx_sg});
    vld(idx_sg,:) = NaN(size(isSolved));
    vld(idx_sg,isSolved) = cellfun(@(x)x.valid,v{idx_sg});

    cost(idx_sg,:) = NaN(size(isSolved));
    cost(idx_sg,isSolved) = cellfun(@(x)x.cost,v{idx_sg});
end
cost = cost.*vld;

%% Plot Results

figure('Name','costs')
tiledlayout(3,4,'TileSpacing','compact')
pos = [1 3 5 7 10];
for idx_sg = 1:n.sg
    nexttile(pos(idx_sg),[1 2])
    hold on
    plot(vec_dP,cost(idx_sg,:),'.','MarkerSize',12);
    plot(vec_dP(isnan(cost(idx_sg,:))),zeros(1, sum(isnan(cost(idx_sg,:)))),'x','Color',[0.6350 0.0780 0.1840],'MarkerSize',8,'Linewidth',1.5)
    ax = gca;
    ax.FontSize = 14;
    xlabel('Pressure drop [Pa]','FontSize',14)
    xlim([0, vec_dP(end)])
    xticks([0, 25, 50])
    ylabel('Cost [kg]','FontSize',14)
    title(strcat('Partition',32,num2str(idx_sg)),'FontSize',14);
    if idx_sg==n.sg
        if sum(isnan(cost(idx_sg,:)))==0
            plot(nan,nan,'x','Color',[0.6350 0.0780 0.1840],'MarkerSize',8,'Linewidth',1.5)
        end
        L = legend('Valid','Invalid','FontSize',12);
        L.ItemTokenSize(1) = 15;
    end
    box on; grid on;hold off
end

end

