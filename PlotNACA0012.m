close all;
load('Results\PlotData_NACA0012', 'y_pred', 'test_y')
figure1 = figure('Unit','Centimeters','Position',[10 10 10 8]);
axes1 = axes('Parent',figure1);
set(axes1,'FontName','Times New Roman');
hold on;box on
scstr={'ko','md','g<','b^','rs'};
name1 = {'Kriging','GEK','WGEK','PGEK','KGEK'};
miny = min([y_pred(:,1); y_pred(:,2); y_pred(:,3); y_pred(:,4); y_pred(:,5); test_y]);
maxy = max([y_pred(:,1);y_pred(:,2); y_pred(:,3); y_pred(:,4); y_pred(:,5); test_y]);
xlimdata = [miny maxy];
fplot(@(x) x, xlimdata,'Linewidth',1.5,'LineStyle','--','Color','k')
sct(1) = scatter(test_y,y_pred(:,1),scstr{1},'DisplayName',name1{1});
sct(2) = scatter(test_y,y_pred(:,2),scstr{2},'DisplayName',name1{1});
sct(3) = scatter(test_y,y_pred(:,5),scstr{3},'DisplayName',name1{2});
sct(4) = scatter(test_y,y_pred(:,3),scstr{4},'DisplayName',name1{3});
sct(5) = scatter(test_y,y_pred(:,4),scstr{5},'filled','DisplayName',name1{4});
xlabel({'Simulated lift coefficient'});
ylabel({'Predicted lift coefficient'});
xlim(xlimdata);ylim(xlimdata);
legend([sct(1) sct(2) sct(3) sct(4)  sct(5)],name1,'location','northwest','Box','on')
print(figure1,'Fig_NACA0012','-dtiff','-r400')
