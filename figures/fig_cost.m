function fig_cost(v,vec_dP,n)
%FIG_COST Summary of this function goes here
%   Detailed explanation goes here

%% Parse Data

vld = zeros(n.sg, n.dP);
cost = zeros(n.sg, n.dP);
for i = 1:n.sg
    vld(i,:) = cellfun(@(x)get_value(x,'valid'),v{i});
    vld(vld==0) = NaN;
    cost(i,:) = cellfun(@(x)get_value(x,'cost'),v{i});
end
cost = cost.*vld;
%% Plot Results

figure('Name','costs')
tiledlayout(3,4,'TileSpacing','compact')
pos = [1 3 5 7 10];
for i = 1:n.sg
    nexttile(pos(i),[1 2])
    hold on
    plot(vec_dP,cost(i,:),'.','MarkerSize',12);
    plot(vec_dP(isnan(cost(i,:))),zeros(1, sum(isnan(cost(i,:)))),'x','Color',[0.6350 0.0780 0.1840],'MarkerSize',8,'Linewidth',1.5)
    %ylim([min(cost,[],'all') max(cost,[],'all')])
    ax = gca;
    ax.FontSize = 14;
    xlabel('Pressure drop [Pa]','FontSize',14)
    xlim([0, vec_dP(end)])
    xticks([0, 25, 50])
    ylabel('Cost [kg]','FontSize',14)
    title(strcat('Partition',32,num2str(i)),'FontSize',14);
    if i==n.sg
        if sum(isnan(cost(i,:)))==0
            plot(nan,nan,'x','Color',[0.6350 0.0780 0.1840],'MarkerSize',8,'Linewidth',1.5)
        end
        L = legend('Valid','Invalid','FontSize',12);
        L.ItemTokenSize(1) = 15;
    end
    box on; grid on;hold off
end

%% Subfunction
function out = get_value(x,value)
    if isa(x,'struct')
        out = x.(value);
    else
        out = NaN;
    end
end
end

