function fig_graph_ppt(G,n,e,nd)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Create labels

nd_label = strings(n.n,1);
nd_label(nd.root,1) = "root";
nd_label(nd.term,1) = "term";

edg_label = strings(n.e,1);
edg_label(e.u) = e.u;

%% Plot figure
figure('Name','Sample Graph')
set(gcf,'Position',[2495,200,658,639])
hold on
clr = [[0.6350 0.0780 0.1840];[0 0.4470 0.7410];[0 0 0];[0 0 0]];
lnsty = ['-','-','-',':'];
for i =1:4
    plot(nan, nan,'Color',clr(i,:),'LineWidth',4,'LineStyle',lnsty(i))
end
L = legend('Feeding\quad', 'Return\quad', 'Bypass\quad', 'User','AutoUpdate','off','Orientation','horizontal','Location','southoutside');
L.ItemTokenSize(1) = 15;
L.FontSize = 20;
%L.Units = 'pixel';
%L.Position(3) = 260;
h = plot(G, 'NodeLabel',nd_label,'EdgeLabel', edg_label,'NodeColor', 'k', 'EdgeColor','k','interpreter','latex','layout','layered','NodeFontSize',20,'EdgeFontSize',20,'Linewidth',4,'MarkerSize',7,'ArrowPosition',.6);
highlight(h,'Edges',e.hot, 'EdgeColor',clr(1,:),'LineWidth',4)
highlight(h,'Edges',e.cold, 'EdgeColor',clr(2,:),'LineWidth',4)
highlight(h,'Edges',e.by, 'EdgeColor',clr(3,:),'LineWidth',4)
highlight(h,'Edges',e.u, 'EdgeColor',clr(4,:),'LineWidth',4,'LineStyle',':')
h.ArrowSize = 15;
set(h,'EdgeAlpha',1)

set(gca,'xtick',[],'ytick',[])
box on
hold off
end

