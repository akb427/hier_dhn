function fig_tempprofiles(tp,temp_prof)
%FIG_TEMPPROFILES  Plots allowed temperature deviations.
%
%   FIG_TEMPPROFILES(tp,temp_prof)
%
%   DESCRIPTION: Plots the allowed temperature deviations of the buildings
%   in the network by building type. 
%
%   INPUTS:
%       tp  - Structure of timings.
%       temp_prof  - Structure of allowable temperature deviations by building type 

%% Create Temperatures
n_prof = 4;
T = 2*ones(n_prof,size(temp_prof.res,2));
T(1,temp_prof.res) = 4;
T(2,temp_prof.com) = 4;
T(3,temp_prof.retail) = 4;
T(4,temp_prof.med) = 4;

leg = {"Residential", "Commercial", "Retail", "Medical"};
clr = lines(n_prof);
%% Plot

% figure('Name','T Profile')
% set(gcf,'Position',[275,153,827.3333333333333,420])
% hold on
% for i = 1:n_prof
%     stairs(tp.time_T_day,T(i,:),'Color',clr(i,:),'Linewidth',1.5);
% end
% L = legend(leg,'AutoUpdate','off');
% L.Location = 'Northwest';
% L.NumColumns = 2;
% L.FontSize = 14;
% for i = 1:n_prof
%     stairs(tp.time_T_day,-T(i,:),'Color',clr(i,:),'Linewidth',1.5);
% end
% ax = gca;
% xlabel('Time','FontSize',14)
% ax.XAxis.SecondaryLabel.Visible='off';
% ax.FontSize = 14;
% ylabel('Flexible Temperature Profiles [$^{\circ} C$]','FontSize',14)
% ylim([-6 6]);
% box on; grid on; hold off

%% Tiled plot
figure('Name','T Profile2')
tiledlayout(2,2);%,'TileSpacing','tight'
set(gcf,'Position',[275,153,827.3333333333333,420])
for i = 1:n_prof
    nexttile
    hold on
    stairs(tp.time_T_day,T(i,:),'Color',clr(i,:),'Linewidth',1.5);
    stairs(tp.time_T_day,-T(i,:),'Color',clr(i,:),'Linewidth',1.5);
    ax = gca;
    xlabel('Time','FontSize',14)
    title(leg{i},'FontSize',16)
    ax.XAxis.SecondaryLabel.Visible='off';
    ax.FontSize = 14;
    ylabel('Profile [$^{\circ} C$]','FontSize',14)
    ylim([-6 6]);
    box on; grid on; hold off
end



end