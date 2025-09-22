function fig_graphrd(G,n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Create labels
n.ered = numedges(G);
edg_label = strings(n.ered,1);
for i = 1:n.ered
    if ~isnan(G.Edges.Names(i))
        edg_label(i) =strcat('$\mathcal{G}_{',num2str(G.Edges.Names(i)),'}$');
    end
end
%% Plot figure
figure('Name','Reduced Graph')
set(gcf,'Position',[410,161,328,356])
hold on
clr = [[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880];[0.3010 0.7450 0.9330]];
h = plot(G, 'NodeLabel',{},'EdgeLabel',edg_label,'NodeColor', 'k', 'EdgeColor','k','interpreter','latex','layout','layered','NodeFontSize',12,'EdgeFontSize',12,'Linewidth',1,'MarkerSize',3,'ArrowPosition',.6);
for i = 1:n.sg
    highlight(h,'Edges',find(G.Edges.Names==i), 'EdgeColor',clr(i,:),'LineWidth',1.5)
end
set(h,'EdgeAlpha',1)

set(gca,'xtick',[],'ytick',[])
box on
hold off
end

