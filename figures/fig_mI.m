function fig_mI(sol_opt, sol_nom, params,e)
%FIG_MI Plot plant flow rate.
%
%   FIG_MI(sol_opt, sol_nom, params,e)
%
%   DESCRIPTION: Plots the plant mass flow rate in both the nominal and
%   optimized cases.
%
%   INPUTS:
%       sol_opt     - Cell of simulation solutions with flexibility. 
%       sol_nom     - Cell of simulation solutions without flexibility. 
%       params  - Structure of parameters.
%       e   - Structure of edge information.

%% Legend Vector
t = datetime(2018,2,1,0,0,0)+seconds(0:params.dt:24*60*60+10*60);
t = t(1:end-1);

%%
figure
hold on
plot(t(1:end-1), sol_nom.mI,'Color',[0.6350 0.0780 0.1840],'Linewidth',1.5)
plot(t, sol_opt.mI,'Color',[0 0.4470 0.7410],'Linewidth',1.5)
plot(t(1:end-1), sum(sol_nom.mdot_e(e.by,:)),'--','Color',[0.6350 0.0780 0.1840],'Linewidth',1.5)
plot(t, sum(sol_opt.mdot_e(e.by,:)),'--','Color',[0 0.4470 0.7410],'Linewidth',1.5)

ylabel('Mass Flow Rate [kg/s]','FontSize',14)
legend('Nominal $\dot{m}_I$','Flexible $\dot{m}_I$','Nominal $\dot{m}_{by}$','Flexible $\dot{m}_{by}$','FontSize',14)
xtickformat('HH:mm')
xlim([t(1) t(end)])
ylim([0 55])
ax = gca;
ax.FontSize = 14;
ax.XAxis.SecondaryLabel.Visible='off';
xlabel('Time','FontSize',14)
box on; grid on; hold off
end

