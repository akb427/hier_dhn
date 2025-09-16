
%%
G = digraph([1 2 2 3 3 4 4 5 5 6 6 7 8 8 9 9 10 11 12 13 14 15],[2 3 15 4 14 5 13 6 8 7 7 12 9 11 10 10 11 12 13 14 15 16])
e.u = [3 5 7 14]; % users without direct bypass connections

%%

cd 'C:\Users\akb42\OneDrive - The Ohio State University\DistrictHeatingNetwork\Optimal Control\rom\figures'
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

fig_graph(G,n,e,nd)
fig_part(G,Gred,n,e,sG)
%fig_graphrd(Gred,n)
fig_demand(params,e,n)
fig_cost(vLL{3},vec_dP,n)
fig_intQ(vsim_all,params,e,n)



set(groot,'defaultAxesTickLabelInterpreter','tex');  
set(groot,'defaulttextinterpreter','tex');
set(groot,'defaultLegendInterpreter','tex');
cd 'C:\Users\akb42\OneDrive - The Ohio State University\DistrictHeatingNetwork\Optimal Control\rom'

%% Calculate Costs

vsim_all = vsim{1};

for i = 2:57
    vsim_all.dP_e = [vsim_all.dP_e vsim{i}.dP_e];
    vsim_all.P_nd = [vsim_all.P_nd vsim{i}.P_nd];
    vsim_all.mdot_e = [vsim_all.mdot_e vsim{i}.mdot_e];
    vsim_all.mI = [vsim_all.mI vsim{i}.mI];
    vsim_all.T = [vsim_all.T vsim{i}.T];
    vsim_all.Qp = [vsim_all.Qp vsim{i}.Qp];
    vsim_all.intQ = [vsim_all.intQ vsim{i}.intQ];
    vsim_all.cost = vsim_all.cost+vsim{i}.cost;
end

% Nominal cost
cnom_mI = sum(vnom.mI(1,1:57));
cflex_mI = sum(vsim_all.mI);

cnom_mby = sum(vnom.mdot_e(e.by, 1:57),'all');

%% Figure
t = 0:10:56*10;

figure
hold on
plot(t, vnom.mI(1:57),'Color',[0.6350 0.0780 0.1840],'Linewidth',1.5)
plot(t, vsim_all.mI,'Color',[0 0.4470 0.7410],'Linewidth',1.5)
plot(t, sum(vnom.mdot_e(e.by,1:57)),'--','Color',[0.6350 0.0780 0.1840],'Linewidth',1.5)
plot(t, sum(vsim_all.mdot_e(e.by,:)),'--','Color',[0 0.4470 0.7410],'Linewidth',1.5)

ylabel('Mass flow rate [kg/s]')
legend('Nominal $\dot{m}_I$','Flexible $\dot{m}_I$','Nominal $\dot{m}_{by}$','Flexible $\dot{m}_{by}$')
ylim([0 60]);
xlabel('Time [min]')
xlim([0 560])
box on; grid on; hold off