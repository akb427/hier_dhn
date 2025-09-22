function fig_mdot_u(sol,params,e,n)
%FIG_MDOT_U Plot user flow rates.
%
%   FIG_MDOT_U(v,params,e,n)
%
%   DESCRIPTION: Plots the mass flow rates delivered to the users, split
%   into residential and commercial buildings
%
%   INPUTS:
%       sol     - Cell of simulation solutions with flexibility. 
%       params  - Structure of parameters.
%       e   - Structure of edge information.
%       n   - Structure of sizes.

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

t = datetime(2018,2,1,0,0,0)+seconds(0:params.dt:24*60*60+10*60);
t = t(1:end-1);
m_res = zeros(n.res, size(sol.mdot_e,2));
m_com = zeros(n.res, size(sol.mdot_e,2));

for i = 1:n.res
    m_res(i,:) = sol.mdot_e(e.res(i)==e.u,:);
end
m_res = m_res(idx_res,:);

for i = 1:n.com
    m_com(i,:) = sol.mdot_e(e.com(i)==e.u,:);
end
m_com = m_com(idx_com,:);
%% Plot

figure('Name','mdot_u')
tiledlayout(1,2,'TileSpacing','tight');
%set(gcf,'Position',[2685,440.3333333333333,894.6666666666665,419.9999999999999])

nexttile
hold on
plot(t,m_res(1:7,:),'Linewidth',1);
plot(t,m_res(8:end,:),'--','Linewidth',1);
legend(leg.res,'Location','Northeast','Autoupdate','off','FontSize',14)
title('Residential','FontSize',16)
ax = gca;
ax.FontSize = 14;
xlabel('Time','FontSize',14)
xtickformat('HH:mm')
xlim([t(1) t(end)])
ax.XAxis.SecondaryLabel.Visible='off';
ylabel('Mass Flow to User [$kg/s$]', 'FontSize',14)
box on; grid on; hold off

nexttile
hold on
plot(t,m_com,'Linewidth',1);
legend(leg.com,'Location','Northeast','Autoupdate','off','FontSize',14)
title('Commercial','FontSize',16)
ax = gca;
ax.FontSize = 14;
xlabel('Time','FontSize',16)
xtickformat('HH:mm')
xlim([t(1) t(end)])
ax.XAxis.SecondaryLabel.Visible='off';
ylabel('Mass flow rate [$kg/s$]','FontSize',14)
box on; grid on; hold off
end

