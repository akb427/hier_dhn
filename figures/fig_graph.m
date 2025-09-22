function fig_graph(G,n,e,v)
%FIG_GRAPH Plots network graph.
%
%   FIG_GRAPH(G,n,e,v)
%
%   DESCRIPTION: Plots the network graph with edge sets labeled as hot,
%   cold, bypass and user. 
%
%   INPUTS:
%       G   - Graph of network.
%       n   - Structure of sizes.
%       e   - Structure of edge information.
%       v   - Structure of node information.

%% Create labels

nd_label = strings(n.v,1);
nd_label(v.root,1) = "root";
nd_label(v.term,1) = "term";

edg_label = strings(n.e,1);
edg_label(e.u) = e.u;

%% Plot figure
figure('Name','Sample Graph')
%set(gcf,'Position',[2389,222,475,460])
hold on
clr = [[0.6350 0.0780 0.1840];[0 0.4470 0.7410];[0.4660 0.6740 0.1880];[0 0 0]];
lnsty = ['-','-','-',':'];
for i =1:4
    plot(nan, nan,'Color',clr(i,:),'LineWidth',2,'LineStyle',lnsty(i))
end
L = legend('Feeding\quad', 'Return\quad', 'Bypass\quad', 'User','AutoUpdate','off','Orientation','horizontal','Location','southoutside');
L.ItemTokenSize(1) = 15;
L.FontSize = 14;
%L.Units = 'pixel';
%L.Position(3) = 260;
h = plot(G, 'NodeLabel',nd_label,'EdgeLabel', edg_label,'NodeColor', 'k', 'EdgeColor','k','interpreter','latex','layout','layered','NodeFontSize',14,'EdgeFontSize',14,'Linewidth',1,'MarkerSize',3,'ArrowPosition',.6);
highlight(h,'Edges',e.hot, 'EdgeColor',clr(1,:),'LineWidth',1.5)
highlight(h,'Edges',e.cold, 'EdgeColor',clr(2,:),'LineWidth',1.5)
highlight(h,'Edges',e.by, 'EdgeColor',clr(3,:),'LineWidth',1.5)
highlight(h,'Edges',e.u, 'EdgeColor',clr(4,:),'LineWidth',1.5,'LineStyle',':')
set(h,'EdgeAlpha',1)

set(gca,'xtick',[],'ytick',[])
box on
hold off
end

