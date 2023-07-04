clear all;
close all;
%%
len = 0.2752;
s = len/10;
Gr = 2.5;
r1 = 1000;
theta1 = pi/3;
Pt = 1;
rho = 1;
r2 = 5;
theta2 = 2*pi/3;
d = sqrt(r1^2+r2^2-r1*r2*cos(theta2-theta1));

delta_r1 = [0:0.1:50];
delta_theta1 = (pi*[0:0.1:30]/360)';

ep_phi = zeros(length(delta_r1),length(delta_theta1));
for i = 1:length(delta_r1)
    for j = 1:length(delta_theta1)
        temp_r1 = delta_r1(i);
        temp_theta1 = delta_theta1(j);
        l = sqrt(r1^2+(r1+temp_r1)^2-2*r1*(r1+temp_r1)*cos(temp_theta1));
        beita = asin((r1+temp_r1)*sin(temp_theta1)/l) - acos((r1^2+d^2-r2^2)/(2*r1*d));
        if l == 0
            ep_d = 0;
        else 
            ep_d = sqrt(d^2 + l^2 - 2*d*l*cos(beita)) - d;
        end
        N = floor((4*pi*r2*rho)/(len*Gr)+0.5);
        ep_phi(i,j) = 2*pi/len*(temp_r1-N*s/2*sin(theta1)*temp_theta1-ep_d);
    end
end
%%
ep_phi = abs(pi-abs(mod(ep_phi,2*pi)-pi));
figure;
imagesc(delta_theta1,delta_r1, ep_phi);
colormap('summer');
a = colorbar;
a.Label.String = 'Phase error';
hold on;
% contour(delta_r1,delta_theta1, ep_phi','--r','ShowText','on','LineWidth',2);
xlabel('\Delta \theta_1');
ylabel('\Delta r_1 (m)');
 hYLabel = get(gca,'YLabel');
 set(hYLabel,'rotation',0,'VerticalAlignment','middle')
set(gca,'fontsize',16);
filename = 'phaseerror1';
save = 1;

% [h, wd, ht] = tightfig();
if save == 1
    name1 = append(filename, '.fig');
    name2 = append(filename, '.pdf');
    saveas(gca, name1);
    exportgraphics(gca, name2);
end
%%
h = 10000;
len = 0.2752;
s = len/10;
Gr = 2.5;
theta1 = pi*[30:0.1:90]/180;
Pt = 1;
rho = 1;
r2 = 5;
theta2 = 2*pi/3;
delta_theta1 = (pi*[-15:0.1:15]/180)';
ep_phi2 = zeros(length(theta1),length(delta_theta1));
for i = 1:length(theta1)
    for j = 1:length(delta_theta1)
        thi = theta1(i);
        temp_theta1 = delta_theta1(j);
        r1 = h*csc(thi);
        d = sqrt(r1^2+r2^2-r1*r2*cos(theta2-thi));
        if temp_theta1 == 0
            ep_phi2(i,j) = 0;
        elseif temp_theta1>0
            temp_r1 = h*csc(thi+temp_theta1) - r1;
            l = sqrt(r1^2+(r1+temp_r1)^2-2*r1*(r1+temp_r1)*cos(temp_theta1));
            beita = real(asin((r1+temp_r1)*sin(temp_theta1)/l) - acos((r1^2+d^2-r2^2)/(2*r1*d)));
            if l == 0
                ep_d = 0;
            else 
                ep_d = sqrt(d^2 + l^2 - 2*d*l*cos(beita)) - d;
            end
            N = floor((4*pi*r2*rho)/(len*Gr)+0.5);
            ep_phi2(i,j) = 2*pi/len*(temp_r1-N*s/2*sin(thi)*temp_theta1-ep_d);
        elseif temp_theta1<0
            temp_r1 = h*csc(thi+temp_theta1) - r1;
            l = sqrt(r1^2+(r1+temp_r1)^2-2*r1*(r1+temp_r1)*cos(abs(temp_theta1)));
            alpha = acos((r1^2+d^2-r2^2)/(2*r1*d));
            beita = pi-thi + alpha;
            ep_d = sqrt(d^2 + l^2 - 2*d*l*cos(beita)) - d;
            N = floor((4*pi*r2*rho)/(len*Gr)+0.5);
            ep_phi2(i,j) = 2*pi/len*(temp_r1-N*s/2*sin(thi)*temp_theta1-ep_d);
        end
    end
end
%%
ep_phi2 = abs(pi-abs(mod(ep_phi2,2*pi)-pi))/pi;
% ep_phi2 = (cos(mod(ep_phi2,2*pi)-pi)+1)/2;
figure;
imagesc(delta_theta1,theta1, ep_phi2);
colormap('summer');
a = colorbar;
hold on;
% contour(delta_r1,delta_theta1, ep_phi','--r','ShowText','on','LineWidth',2);
xlabel('\Delta \theta_1');
ylabel('\theta_1');
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
hL = ylabel(a,'D_\epsilon');    
hL.Position = [1 0 0];
set(hL,'rotation',0);
set(gca,'fontsize',16);
filename = 'phaseerror2';
save = 1;

% [h, wd, ht] = tightfig();
if save == 1
    name1 = append(filename, '.fig');
    name2 = append(filename, '.pdf');
    saveas(gca, name1);
    exportgraphics(gca, name2);
end
%%
h = 10000;
len = 0.2752;
s = len/10;
Gr = 2.5;
theta1 = pi*[30:0.1:90]/180;
Pt = 1;
rho = 1;
r2 = 5;
theta2 = 2*pi/3;
ep = [-10:0.1:10];
ep_phi2 = zeros(length(theta1),length(ep));
d_theta1 = [];
for i = 1:length(theta1)
    for j = 1:length(ep)
        thi = theta1(i);
        r1 = h*csc(thi);
        epi = ep(j);
        if epi == 0
            ep_phi2(i,j) = 0;
        elseif epi > 0
            temp_r1 =sqrt(epi^2 + r1^2 - 2*epi*r1*cos(thi)) - r1;
            temp_theta1 = asin(h/(temp_r1+r1)) - thi;
            d = sqrt(r1^2+r2^2-r1*r2*cos(theta2-thi));
            l = epi;
            beita = real(asin((r1+temp_r1)*sin(temp_theta1)/l) - acos((r1^2+d^2-r2^2)/(2*r1*d)));
            if l == 0
                ep_d = 0;
            else 
                ep_d = sqrt(d^2 + l^2 - 2*d*l*cos(beita)) - d;
            end
            N = floor((4*pi*r2*rho)/(len*Gr)+0.5);
            ep_phi2(i,j) = 2*pi/len*(temp_r1-N*s/2*sin(thi)*temp_theta1-ep_d);
        elseif epi < 0
            l = abs(epi);
            temp_r1 =sqrt(epi^2 + r1^2 - 2*l*r1*cos(pi-thi)) - r1;
            temp_theta1 = abs(asin(h/(temp_r1+r1)) - thi);
            d = sqrt(r1^2+r2^2-r1*r2*cos(theta2-thi));
            alpha = acos((r1^2+d^2-r2^2)/(2*r1*d));
            beita = pi-thi + alpha;
            ep_d = sqrt(d^2 + l^2 - 2*d*l*cos(beita)) - d;
            N = floor((4*pi*r2*rho)/(len*Gr)+0.5);
            ep_phi2(i,j) = 2*pi/len*(temp_r1-N*s/2*sin(thi)*temp_theta1-ep_d);
        end
        if thi ==pi/2 && epi ~= 0 && mod(epi,5) == 0
            d_theta1 = [d_theta1, temp_theta1];
        end
    end
end
%%
ep_phi2 = abs(pi-abs(mod(ep_phi2,2*pi)-pi))/pi;
% ep_phi2 = (cos(mod(ep_phi2,2*pi)-pi)+1)/2;
figure;
imagesc(ep,theta1, ep_phi2);
colormap('summer');
a = colorbar;
hold on;
contour(ep,theta1, ep_phi2,'--r','ShowText','on','LineWidth',2);
xlabel('\epsilon (m)');
ylabel('\theta_1');
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle');
hL = ylabel(a,'D_\epsilon');    
hL.Position = [1 0 0];
set(hL,'rotation',0);
set(gca,'fontsize',16);
filename = 'phaseerror3';
save = 1;

% [h, wd, ht] = tightfig();
if save == 1
    name1 = append(filename, '.fig');
    name2 = append(filename, '.pdf');
    saveas(gca, name1);
    exportgraphics(gca, name2);
end
%%
a = [0:0.1:2*pi]; 
b = sin(a);
b1 = sin(a) + sin(a);
b2 = sin(a) + sin(a + pi);
b3 = sin(a) + sin(a + pi + pi/2);
b4 = sin(a) + sin(a + pi - pi/2);
plot(a,b,a,b1,a,b2,a,b3,a,b4);
% plot(a, b, 'r--',a,sin(a + pi),'g*',a,b2,'b')
% legend("RFI","Reflected waveform with D_\epsilon = 0", "Received waveform")