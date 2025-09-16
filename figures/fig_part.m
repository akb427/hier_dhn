function fig_part(G,Gred,n,e,sG)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Create labels
edg_label1 = strings(n.e,1);
edg_label1(e.u) = e.u;

n.ered = numedges(Gred);
edg_label2 = strings(n.ered,1);
for i = 1:n.ered
    if ~isnan(Gred.Edges.Names(i))
        edg_label2(i) =strcat('$\mathcal{G}_{',num2str(Gred.Edges.Names(i)),'}$');
    end
end

%% Plot figure
figure('Name','Partition')
set(gcf,'Position',[360,169,727,348.6666666666665])
tiledlayout(1,3,'TileSpacing','tight');

nexttile([1 2])
hold on
clr = [[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880];[0.3010 0.7450 0.9330]];
leg = strings(1,n.sg);
for i =1:n.sg
    plot(nan, nan,'Color',clr(i,:),'LineWidth',2)
    leg(i) = strcat('Partition',32,num2str(i));
end
L = legend(leg,'AutoUpdate','off','FontSize',14,'Orientation','horizontal','Interpreter','latex','NumColumns',5);
L.Layout.Tile = 'south';
L.ItemTokenSize(1) = 15;
%L.Units = 'pixel';
%L.Position(3) = 260;
h = plot(G, 'NodeLabel',{},'EdgeLabel',edg_label1,'NodeColor', 'k', 'EdgeColor','k','interpreter','latex','layout','layered','NodeFontSize',14,'EdgeFontSize',14,'Linewidth',1,'MarkerSize',3,'ArrowPosition',.6);
for i = 1:n.sg
    highlight(h,'Edges',sG{i}.Edges.Names, 'EdgeColor',clr(i,:),'LineWidth',1.5)
end
set(h,'EdgeAlpha',1)

set(gca,'xtick',[],'ytick',[])
box on
hold off

nexttile([1 1])
hold on
clr = [[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880];[0.3010 0.7450 0.9330]];
h = plot(Gred, 'NodeLabel',{},'EdgeLabel',edg_label2,'NodeColor', 'k', 'EdgeColor','k','interpreter','latex','layout','layered','NodeFontSize',14,'EdgeFontSize',14,'Linewidth',1,'MarkerSize',3,'ArrowPosition',.6);
for i = 1:n.sg
    highlight(h,'Edges',find(Gred.Edges.Names==i), 'EdgeColor',clr(i,:),'LineWidth',1.5)
end
set(h,'EdgeAlpha',1)

set(gca,'xtick',[],'ytick',[])
box on
hold off
end

