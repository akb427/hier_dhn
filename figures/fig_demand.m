function fig_demand(params_all,e,n,tp,idx_end)
%FIG_DEMAND Plot building demands.
%
%   FIG_DEMAND(params_all,e,n,tp,idx_end)
%
%   DESCRIPTION: Plots the building demands over time, seperated by
%   building type and labeled according to edge number.
%   
%   INPUTS:
%       params_all  - Structure of network-wide parameters.
%       e   - Structure of edge information.
%       n   - Structure of sizes.
%       tp  - Structure of timings.
%       idx_end - Upper limit for time plotting

%% Legend Vector
% Residential
bld.res = [3561,3801,4017,4090,4177,20041,28770,22719,80368,80372,80387];
e.res = [8,12,22,24,25,57,59,60,20,10,13];
[e.res_sort, idx_res] = sort(e.res);
n.res = numel(bld.res);
leg.res = cell(1,n.res);
for i = 1:n.res
    leg.res{i} = strcat('$e_{',num2str(e.res_sort(i)),'}$');
end
% Commercial
bld.com = [1700,18740,95364,123604,177428,232839,343832];
e.com = [38,49,54,51,32,41,44];
[e.com_sort, idx_com] = sort(e.com);
n.com = numel(bld.com);
leg.com = cell(1,n.com);
for i = 1:n.com
    leg.com{i} = strcat('$e_{',num2str(e.com_sort(i)),'}$');
end

for i = 1:n.res
    Qres(i,:) = params_all.Qb(e.res(i)==e.u,1:idx_end);
end
Qres=Qres(idx_res,:);
for i = 1:n.com
    Qcom(i,:) = params_all.Qb(e.com(i)==e.u,1:idx_end);
end
Qcom=Qcom(idx_com,:);

%% Plot
figure('Name','Qb')
tiledlayout(2,1,'TileSpacing','tight');
set(gcf,'Position',[275,153,827.3333333333333,420])
nexttile
hold on
plot(tp.time_day(1:idx_end),Qres(1:7,:),'Linewidth',1.5);
plot(tp.time_day(1:idx_end),Qres(8:end,:),'--','Linewidth',1.5);
L = legend(leg.res);
L.Location = 'Northwest';
L.NumColumns = 2;
L.FontSize = 14;
title('Residential','FontSize',16)
ax = gca;
xlabel('Time','FontSize',14)
ax.XAxis.SecondaryLabel.Visible='off';
ax.FontSize = 14;
ylabel('Nominal Heat Demand [$kW$]','FontSize',14)
ylim([0 100]);
box on; grid on; hold off

nexttile
hold on
plot(tp.time_day(1:idx_end),Qcom,'Linewidth',1.5);
legend(leg.com,'Location','Northwest','FontSize',14)
hold on
ax = gca;
ax.FontSize = 14;
title('Commercial','FontSize',16)
xlabel('Time','FontSize',14)
ax.XAxis.SecondaryLabel.Visible='off';
ylabel('Nominal Heat Demand [$kW$]','FontSize',14)
box on; grid on; hold off

end

