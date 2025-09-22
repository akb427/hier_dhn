function fig_intQ(sol,params,params_all,e,n)
%FIG_INTQ  Plot used flexibility.
%
%   fig_intQ(sol,params,params_all,e,n)
%
%   DESCRIPTION: Plots the flexibility used by the buildings in the network
%   split into residential and commercial with upper and lower limits
%   labeled. 
%
%   INPUTS:
%       sol  - Cell of simulation solutions with flexibility. 
%       params      - Structure of parameters.
%       params_all  - Structure of network-wide parameters.
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

n_comp = size(sol.Qp,2);
t = datetime(2018,2,1,0,0,0)+seconds(0:params.dt_T:n_comp*params.dt_T);
t = t(1:end-1);
Qres = zeros(n.res, n_comp);
Qcom = zeros(n.res, n_comp);
Cres_u = zeros(n.res, n_comp);
Ccom_u = zeros(n.com, n_comp);
Cres_l = zeros(n.res, n_comp);
Ccom_l = zeros(n.com, n_comp);
for i = 1:n.res
    Qres(i,:) = sol.intQ(e.res(i)==e.u,:);
    Cres_u(i,:) = params_all.Cap_u(e.res(i)==e.u,1:n_comp);
    Cres_l(i,:) = params_all.Cap_l(e.res(i)==e.u,1:n_comp);
end
Qres = Qres(idx_res,:);
Cres_u = Cres_u(idx_res,:);
Cres_l = Cres_l(idx_res,:);

for i = 1:n.com
    Qcom(i,:) = sol.intQ(e.com(i)==e.u,:);
    Ccom_u(i,:) = params_all.Cap_u(e.com(i)==e.u,1:n_comp);
    Ccom_l(i,:) = params_all.Cap_l(e.com(i)==e.u,1:n_comp);
end
Qcom = Qcom(idx_com,:);
Ccom_u = Ccom_u(idx_com,:);
Ccom_l = Ccom_l(idx_com,:);

%% Plot

clr = lines(7);
figure('Name','intQ')
tiledlayout(1,2,'TileSpacing','tight');
%set(gcf,'Position',[2685,440.3333333333333,894.6666666666665,419.9999999999999])

nexttile
hold on
plot(t,Qres(1:7,:),'Linewidth',1.5);
plot(t,Qres(8:end,:),'--','Linewidth',1.5);
legend(leg.res,'Location','Northeast','Autoupdate','off','FontSize',14)
for i = 1:7
    stairs(t,Cres_u(i,:),'Color',clr(i,:));
    stairs(t,Cres_l(i,:),'Color',clr(i,:));
end
for i = 1:4
    stairs(t,Cres_u(i+7,:),'--','Color',clr(i,:));
    stairs(t,Cres_l(i+7,:),'--','Color',clr(i,:));
end
title('Residential','FontSize',16)
ax = gca;
ax.FontSize = 14;
xlabel('Time','FontSize',14)
xtickformat('HH:mm')
xlim([t(1) t(end)])
ax.XAxis.SecondaryLabel.Visible='off';
ylabel('Used Flexibility [$MJ$]','FontSize',14)
box on; grid on; hold off

nexttile
hold on
plot(t,Qcom,'Linewidth',1.5);
legend(leg.com,'Location','Northeast','Autoupdate','off','FontSize',14)
for i = 1:7
    stairs(t,Ccom_u(i,:),'Color',clr(i,:));
    stairs(t,Ccom_l(i,:),'Color',clr(i,:));
end
ax = gca;
ax.FontSize = 14;
title('Commercial','FontSize',16)
xlabel('Time','FontSize',14)
ylabel('Used Flexibility [$MJ$]','FontSize',14)
xtickformat('HH:mm')
xlim([t(1) t(end)])
ax.XAxis.SecondaryLabel.Visible='off';
box on; grid on; hold off
end

