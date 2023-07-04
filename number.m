clear all;
close all;

r1 = [1:0.1:20];
len = 0.2752;
Gr = 2.5;
rho = [0.1:0.1:10]';
N = floor((4*pi.*r1.*rho)/(len*Gr)+0.5);

figure;
imagesc(r1,rho, N);
colormap('summer');
a = colorbar;
hL = ylabel(a,'N');    
hL.Position = [0.8 0 0];
set(hL,'rotation',0);
hold on;
contour(r1,rho, N,'--r','ShowText','on','LineWidth',2);
xlabel('r_{2} (m)');
ylabel('\rho');
set(gca,'fontsize',16);
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle');
filename = 'elementnumber';
save = 1;

% [h, wd, ht] = tightfig();
if save == 1
    name1 = append(filename, '.fig');
    name2 = append(filename, '.pdf');
    saveas(gca, name1);
    exportgraphics(gca, name2);
end
% print -opengl -dpdf -r600 hst_1u_1ms