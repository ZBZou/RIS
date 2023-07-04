clear all;
close all;
%%
%%Time specifications:
   Fs = 3e11;                   % samples per second
   dt = 1/Fs;                    % seconds per sample
   
   StopTime = 5e-8;             % seconds
   t = (0:dt:StopTime-dt)';     
   
   %%Sine wave:
   fc = 1090e6;
%    2.4e9; 
%    1090M
   x = sin(2*pi*fc*t);
%    plot(t,x);
   
   
%    sampleRate = 20e6;   
%    maxDopplerShift  = 0; % Maximum Doppler shift of diffuse components (Hz)
%    delayVector = (0:6:12)*1e-8; % Discrete delays of four-path channel (s)
%    gainVector  = [0 -6 -12]; % Average path gains (dB)
% 
%    % Configure a Rayleigh channel object
%    rayChan = comm.RayleighChannel( ...
%    'SampleRate',sampleRate, ...
%    'PathDelays',delayVector, ...
%    'AveragePathGains',gainVector, ...
%    'MaximumDopplerShift',maxDopplerShift, ...
%    'RandomStream','mt19937ar with seed', ...
%    'Seed',10, ...
%    'PathGainsOutputPort',true);  

   
 
   c = 3e8;                 % light speed
   lambda = c / fc;             % wave length
   e_size =  lambda/4;            % size of each element
   spacing = lambda / 10;       % interval between two elements
   
   h = 10000;
   theta1 = deg2rad(60);
   r1 = h/sin(theta1);

   r2 =  20;
   theta2 = acos(-12/20);
   
   G_Tx = 10^(11.2/10);
   G_r = 10^(4/10);    % gain of each element
   
   d = sqrt(r1*r1+r2*r2-r1*r2*cos(theta1-theta2))          ;      % direct path length
   
   G_Rx  =   10^(50/10);            % total gain of beamforming pattern
   G_Rx_ =   10^(50/10);
  

   Rho =    sqrt(G_Rx/G_Rx_);           % ratio gain 
   N =  round((4*pi*r1*r2*Rho)/(lambda*d*G_r));  % total number of elements
   
   P_Tx = 100;
   P_Rx = P_Tx*G_Tx*G_Rx*lambda*lambda/((4*pi*d)^2);      %  received power from direct path
   P_Rx_ = (P_Tx* N *N  * G_Tx*G_Rx_* G_r*G_r* lambda.^4) /(16.^2*pi.^4*r1.^2*r2.^2);       %received power from reflected path
   
   S = 0;
   S1= 0;
%    phase_delay_all = [];
   phase_delay_RIS = [];
%    for n = 1:N
%        d_n = r1+r2+(n-1)*spacing*(sin(pi/2-theta1)-sin(pi/2+theta2))-d;
% %        d_n = r1+distance+(n-1)*spacing*(sin(pi/2-theta1))-d;
%        phi_n = 2*pi*d_n/lambda; 
%        K = 0;
%        tau_n = ((2*K+1)*lambda-2*d_n)/(2*c);
%        while(tau_n<0)
%            K=K+1;
%            tau_n = ((2*K+1)*lambda-2*d_n)/(2*c);
%        end  
%        phase_n = c*tau_n/lambda*2*pi;
%        phase_RIS = mod(d+d_n,lambda)/lambda*2*pi;
%        phase_delay_all = mod(phase_RIS+phase_n, 2*pi);
%        phase_delay_RIS = [phase_delay_RIS,phase_n];
%        S = S + sqrt(P_Rx_/N/N)* sin(2*pi*fc*t + mod(phase_RIS+phase_n, 2*pi));
% %        S1 = S + (sqrt(P_Rx_/N/N)* sin(2*pi*fc*t + mod(phase_RIS, 2*pi))).^2;
%    end
%    
%    
% %    y = rayChan(x);
% 
%    Phase_rece = mod(d,lambda)/lambda*2*pi;
% %    x = sin(2*pi*fc*t + Phase_rece).^2;
%    x = sqrt(P_Rx) * sin(2*pi*fc*t + Phase_rece);
% %    x = awgn(x,5,'measured');
% %    S = awgn(S,5,'measured');
%    
% %    figure;
% % %    subplot(3,1,1);
% %    plot(t,S1);
% % %    title("received signal");   
%    
%    
% % %    figure;
% %    subplot(3,1,1);
% %    plot(t,x);
% %    title("received signal");
% %    
% % %    figure;
% %    subplot(3,1,2);
% %    plot(t,S);
% %    title("reflected signal");
% %    
% % %    figure;
% %    subplot(3,1,3);
% %    plot(t,S+x);
% %    title("added");
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = [-20:0.5:20];
Y = [0:0.5:20];
% Distance = zeros([2001 2001]);
Z = zeros([81 41]);
N1 = round(sqrt(N));
position_RIS = zeros([1 N1]);
ZdB = zeros([81 41]);
for k = 1:N1
    position_RIS(1,k) = (k-N1)*(spacing+e_size);
end

xposition_Tx = position_RIS(1,N1) + r1*cos(theta1);
yposition_Tx = r1 * sin(theta1);
xposition_Rx = position_RIS(1,N1) - r2*cos(pi-theta2);
yposition_Rx = r2 * sin(pi-theta2);
d = sqrt((xposition_Tx - xposition_Rx)^2 + (yposition_Tx - yposition_Rx)^2);
%%
% for i = 1:81
%     disp(i);
%     for j = 1:41
%         energy = 0;
%         energy2 = 0;
%         for n = 1:N1  
% % %             disp(j);
% %             distance = sqrt(Y(1,(41+1-j)).^2 + (position_RIS(1,n)-X(1,i)).^2);
% %             d_n = r1+distance-d;
% % %             +(N1+n-1)*spacing*(sin(pi/2-theta1));
% % %             sqrt(yposition_Tx.^2 + (position_RIS(1,n)-xposition_Tx).^2)+distance-d;
% % %             r1+distance-d+(N1+n-1)*spacing*(sin(pi/2-theta1));
% % %                (N1+n-1)*spacing*(sin(pi/2-theta1))-d;
% %                phi_n = 2*pi*d_n/lambda; 
% %                K = 0;
% %                tau_n = ((2*K+1)*lambda-2*d_n)/(2*c);
% %                while(tau_n<0)
% %                    K=K+1;
% %                    tau_n = ((2*K+1)*lambda-2*d_n)/(2*c);
% %                end  
% %                phase_n = c*tau_n/lambda*2*pi;
% % %                phase_RIS = mod(d+d_n,lambda)/lambda*2*pi;
% % %                phase_delay_all = mod(phase_RIS+phase_n, 2*pi);
% %                phase_delay_RIS = [phase_delay_RIS,phase_n];
% %                
% %        
% %             phase = phase_delay_RIS(1,n) + (r1+spacing*sin(pi/2-theta1)*(N1-n)*2*pi)/lambda+ distance/lambda*2*pi; 
% % %             phase = phase_delay_RIS(1,n) + sqrt(yposition_Tx.^2 + (position_RIS(1,n)-xposition_Tx).^2)/lambda+ distance/lambda*2*pi;            
% %             magnitude = (G_r*lambda*lambda)/(4*pi*distance)/(4*pi*distance);
% %             a = sqrt(magnitude)*sin(phase);
% %             energy = energy+ a;
% %             
% %  
% % %             dis = sqrt(yposition_Tx.^2 + (position_RIS(1,n)-xposition_Tx).^2);
% % %             b = sqrt(G_r*lambda*lambda/(4*pi*dis)/(4*pi*dis));
% % % %             b = sqrt(G_r/dis/dis);
% % %             energy = energy+ a*b;
% %         
% %         end
% %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %         for n = 1:N1
% %             dis = sqrt(yposition_Tx.^2 + (position_RIS(1,n)-xposition_Tx).^2);
% %             b = sqrt(P_Tx*G_Tx*G_r*lambda*lambda/(4*pi*dis)/(4*pi*dis));
% % %             b = sqrt(G_r/dis/dis);
% %             energy2 = energy2+ b;
% %         end
% % %         Z(i,j) = (abs(energy*energy2))^2;
% % %         ZdB(i,j) = 10*log10((abs(energy*energy2))^2);
% %         ZdB(i,j) = 10*log10(abs(energy)^2);
%         %             disp(j);
%             distance = sqrt(Y(1,(41+1-j)).^2 + (position_RIS(1,n)-X(1,i)).^2);
% %             Distance(i,j) = distance;
% 
%                d_n = r1+distance-d;
% %                (N1+n-1)*spacing*(sin(pi/2-theta1))-d;
%                phi_n = 2*pi*d_n/lambda; 
%                K = 0;
%                tau_n = ((2*K+1)*lambda-2*d_n)/(2*c);
%                while(tau_n<0)
%                    K=K+1;
%                    tau_n = ((2*K+1)*lambda-2*d_n)/(2*c);
%                end  
%                phase_n = c*tau_n/lambda*2*pi;
% %                phase_RIS = mod(d+d_n,lambda)/lambda*2*pi;
% %                phase_delay_all = mod(phase_RIS+phase_n, 2*pi);
%                phase_delay_RIS = [phase_delay_RIS,phase_n];
%                
%        
%             phase = phase_delay_RIS(1,n) + (r1+spacing*sin(pi/2-theta1)*(N1-n)*2*pi)/lambda+ distance/lambda*2*pi;            
%             magnitude = (G_r*lambda*lambda)/(4*pi*distance)/(4*pi*distance);
%             a = sqrt(magnitude)*sin(phase);
%             energy = energy+ a;
%         end
%         ZdB(i,j) = 10*log10(abs(energy)^2);
%     end
% end
% phase_delay_RIS = [];
% for n = 1:N1  
% %   disp(j);
%     dis = sqrt(yposition_Tx.^2 + (position_RIS(1,n)-xposition_Tx).^2);
%     distance_rx = sqrt(yposition_Rx.^2 + (position_RIS(1,n)-xposition_Rx).^2);
%     d_n = dis+distance_rx-d; 
% %                +(N1-n)*spacing*(sin(pi/2-theta1))-d;
% %                d_n = dis+distance-d;
% %                +(N1+n-1)*spacing*(sin(pi/2-theta1))-d;
%    phi_n = mod(2*pi*(dis+distance_rx)/lambda,2*pi) - mod(2*pi*(d)/lambda,2*pi); 
%    if phi_n <= pi
%        tau_n = (pi - phi_n)*lambda/(2*pi*c);
%    else
%        tau_n = (3*pi - phi_n)*lambda/(2*pi*c);
%    end
%        
% %                K = 0;
% %                tau_n = ((2*K+1)*lambda-2*d_n)/(2*c);
% %                while(tau_n<0)
% %                    K=K+1;
% %                    tau_n = ((2*K+1)*lambda-2*d_n)/(2*c);
% %                end  
%        phase_n = c*tau_n/lambda*2*pi;
% %                phase_RIS = mod(d+d_n,lambda)/lambda*2*pi;
% %                phase_delay_all = mod(phase_RIS+phase_n, 2*pi);
%       
%        
% %             dis = sqrt(yposition_Tx.^2 + (position_RIS(1,n)-xposition_Tx).^2);
% %             phase = phase_delay_RIS(1,n) + dis/lambda*2*pi+ distance/lambda*2*pi;
%     phase = phase_n + mod(((dis+distance_rx)*2*pi)/lambda,2*pi);            
% %     magnitude = (G_r*G_Rx_*lambda*lambda)*(P_Tx*G_Tx*G_r*lambda*lambda)/(4*pi*dis)^2/(4*pi*distance_rx)^2;
%     
%     
% %             b = sqrt(P_Tx*G_Tx*G_r*lambda*lambda/(4*pi*dis)/(4*pi*dis));
%     
% %     a = sqrt(magnitude*N1*N1)*sin(phase);
% %     energy = energy+ a;
%     phase_delay_RIS = [phase_delay_RIS,phase_n]; 
% end
%%
T = zeros([81 41]);
TdB = zeros([81 41]);
R = zeros([81 41]);
RdB = zeros([81 41]);
rx = [17,9];
for i = 1:81
%     disp(i);
    for j = 1:41
        a_ris = 0;
        energy = 0;
        energy2=0;
%         disp(j);
        distance_t = sqrt((Y(1,(41+1-j))-yposition_Tx).^2 + (xposition_Tx-X(1,i)).^2);
        phase_t = mod(distance_t*2*pi/lambda,2*pi);         
        magnitude_t = sqrt((G_Tx*G_Rx*P_Tx*lambda*lambda)/(4*pi*distance_t)/(4*pi*distance_t));
%         magnitude = 1/distance/distance;
        a_t = magnitude_t*sin(phase_t);
        T(i,j) = abs(a_t)^2;
        TdB(i,j) = 10*log10(abs(a_t)^2);
        for n = 1:N1  
%             disp(j);
            distance = sqrt(Y(1,(41+1-j)).^2 + (position_RIS(1,n)-X(1,i)).^2);
            dis = sqrt(yposition_Tx.^2 + (position_RIS(1,n)-xposition_Tx).^2);
            distance_rx = sqrt(yposition_Rx.^2 + (position_RIS(1,n)-xposition_Rx).^2);
            d_n = dis+distance_rx-d; 
            phi_n = mod(2*pi*(dis+distance_rx)/lambda,2*pi) - mod(2*pi*(d)/lambda,2*pi); 
            if phi_n <= pi
                tau_n = (pi - phi_n)*lambda/(2*pi*c);
            else
                tau_n = (3*pi - phi_n)*lambda/(2*pi*c);
            end
            phase_n = c*tau_n/lambda*2*pi;
%             d_n = dis+distance-d; 
% %                +(N1-n)*spacing*(sin(pi/2-theta1))-d;
% %                d_n = dis+distance-d;
% %                +(N1+n-1)*spacing*(sin(pi/2-theta1))-d;
%                phi_n = mod(2*pi*(dis+distance)/lambda,2*pi) - mod(2*pi*(d)/lambda,2*pi); 
%                if phi_n <= pi
%                    tau_n = (pi - phi_n)*lambda/(2*pi*c);
%                else
%                    tau_n = (3*pi - phi_n)*lambda/(2*pi*c);
%                end
               
%                K = 0;
%                tau_n = ((2*K+1)*lambda-2*d_n)/(2*c);
%                while(tau_n<0)
%                    K=K+1;
%                    tau_n = ((2*K+1)*lambda-2*d_n)/(2*c);
%                end  
%                phase_n = c*tau_n/lambda*2*pi;
%                phase_RIS = mod(d+d_n,lambda)/lambda*2*pi;
%                phase_delay_all = mod(phase_RIS+phase_n, 2*pi);
              
               
%             dis = sqrt(yposition_Tx.^2 + (position_RIS(1,n)-xposition_Tx).^2);
%             phase = phase_delay_RIS(1,n) + dis/lambda*2*pi+ distance/lambda*2*pi;
            phase = phase_n + mod(((dis+distance)*2*pi)/lambda,2*pi);    
            phase = mod(phase,2*pi);
            magnitude = sqrt((G_r*G_Rx_*lambda*lambda)*(P_Tx*G_Tx*G_r*lambda*lambda)/(4*pi*r1)^2/(4*pi*r2)^2);
            
            
%             b = sqrt(P_Tx*G_Tx*G_r*lambda*lambda/(4*pi*dis)/(4*pi*dis));
            
            a = N1*magnitude*sin(phase);
            a_ris = a_ris+ a;
        end
%         energy=energy  - (N-N1^2)*magnitude*sin(phase);
%         phase_delay_RIS = [phase_delay_RIS,phase_n]; 

%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         for n = 1:N1
%             dis = sqrt(yposition_Tx.^2 + (position_RIS(1,n)-xposition_Tx).^2);
%             b = sqrt(P_Tx*G_Tx*G_r*lambda*lambda/(4*pi*dis)/(4*pi*dis));
% %             b = sqrt(G_r/dis/dis);
%             energy2 = energy2+ b;
%         end
        Z(i,j) = abs(a_ris)^2;
        ZdB(i,j) = 10*log10((abs(a_ris)^2));
        a_r = a_ris + a_t;
        R(i,j) = abs(a_r)^2;
        RdB(i,j) = 10*log10((abs(a_r)^2));
    end
end
%%

% x1 = [-12.1:0.001:-11.9];
% y1 = [15.9:0.001:16.1];
x1 = [-20:0.1:20];
y1 = [0:0.1:20];
T = zeros([length(x1) length(y1)]);
TdB = zeros([length(x1) length(y1)]);
R = zeros([length(x1) length(y1)]);
RdB = zeros([length(x1) length(y1)]);
rx = [17,9];
for i = 1:length(x1)
%     disp(i);
    for j = 1:length(y1)
        a_ris = 0;
        energy = 0;
        energy2 = 0;
%         disp(j);
        distance_t = sqrt((y1(j)-yposition_Tx).^2 + (xposition_Tx-x1(i)).^2);
        phase_t = mod(distance_t*2*pi/lambda,2*pi);         
        magnitude_t = sqrt((G_Tx*G_Rx*P_Tx*lambda*lambda)/(4*pi*distance_t)/(4*pi*distance_t));
%         magnitude = 1/distance/distance;
        a_t = magnitude_t*sin(phase_t);
        T(i,j) = abs(a_t)^2;
        TdB(i,j) = 10*log10(abs(a_t)^2);
        for n = 1:N1  
%             disp(j);
            distance = sqrt(y1(j)^2 + (position_RIS(1,n)-x1(i))^2);
            dis = sqrt(yposition_Tx.^2 + (position_RIS(1,n)-xposition_Tx).^2);
            distance_rx = sqrt(yposition_Rx.^2 + (position_RIS(1,n)-xposition_Rx).^2);
            d_n = dis+distance_rx-d; 
            phi_n = mod(2*pi*(dis+distance_rx)/lambda,2*pi) - mod(2*pi*(d)/lambda,2*pi); 
            if phi_n <= pi
                tau_n = (pi - phi_n)*lambda/(2*pi*c);
            else
                tau_n = (3*pi - phi_n)*lambda/(2*pi*c);
            end
            phase_n = c*tau_n/lambda*2*pi;
%             d_n = dis+distance-d; 
% %                +(N1-n)*spacing*(sin(pi/2-theta1))-d;
% %                d_n = dis+distance-d;
% %                +(N1+n-1)*spacing*(sin(pi/2-theta1))-d;
%                phi_n = mod(2*pi*(dis+distance)/lambda,2*pi) - mod(2*pi*(d)/lambda,2*pi); 
%                if phi_n <= pi
%                    tau_n = (pi - phi_n)*lambda/(2*pi*c);
%                else
%                    tau_n = (3*pi - phi_n)*lambda/(2*pi*c);
%                end
               
%                K = 0;
%                tau_n = ((2*K+1)*lambda-2*d_n)/(2*c);
%                while(tau_n<0)
%                    K=K+1;
%                    tau_n = ((2*K+1)*lambda-2*d_n)/(2*c);
%                end  
%                phase_n = c*tau_n/lambda*2*pi;
%                phase_RIS = mod(d+d_n,lambda)/lambda*2*pi;
%                phase_delay_all = mod(phase_RIS+phase_n, 2*pi);
              
               
%             dis = sqrt(yposition_Tx.^2 + (position_RIS(1,n)-xposition_Tx).^2);
%             phase = phase_delay_RIS(1,n) + dis/lambda*2*pi+ distance/lambda*2*pi;
            phase = phase_n + mod(((dis+distance)*2*pi)/lambda,2*pi);    
            phase = mod(phase,2*pi);
            magnitude = sqrt((G_r*G_Rx_*lambda*lambda)*(P_Tx*G_Tx*G_r*lambda*lambda)/(4*pi*r1)^2/(4*pi*r2)^2);
            
            
%             b = sqrt(P_Tx*G_Tx*G_r*lambda*lambda/(4*pi*dis)/(4*pi*dis));
            
            a = N1*magnitude*sin(phase);
            a_ris = a_ris+ a;
        end
       
%         energy=energy  - (N-N1^2)*magnitude*sin(phase);
%         phase_delay_RIS = [phase_delay_RIS,phase_n]; 

%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         for n = 1:N1
%             dis = sqrt(yposition_Tx.^2 + (position_RIS(1,n)-xposition_Tx).^2);
%             b = sqrt(P_Tx*G_Tx*G_r*lambda*lambda/(4*pi*dis)/(4*pi*dis));
% %             b = sqrt(G_r/dis/dis);
%             energy2 = energy2+ b;
%         end
        Z(i,length(y1)-j+1) = abs(a_ris)^2;
        ZdB(i,length(y1)-j+1) = 10*log10((abs(a_ris)^2));
        a_r = a_ris + a_t;
        R(i,length(y1)-j+1) = abs(a_r)^2;
        RdB(i,length(y1)-j+1) = 10*log10((abs(a_r)^2));
    end
end

%%
figure;
imagesc(x1,y1,RdB'); 
colormap('hot');
colorbar('vert');
set(gca,'YTick',[15.9 15.95 16 16.05 16.1],'FontSize', 16);
set(gca,'YTickLabel',[16.1 16.05 16 15.95 15.9],'FontSize', 16);
set(gca,'XTick',[-12.1 -12.05 -12 -11.95 -11.9],'FontSize', 16);
set(gca,'XTickLabel',[-12.1 -12.05 -12 -11.95 -11.9],'FontSize', 16);
xlabel('X(m)','FontSize', 16);
ylabel('Y(m)','FontSize', 16,'rotation',0);
% set(gca,'YTick',[0:1:10]);
% set(gca,'YTickLabel',[10 9 8 7 6 5 4 3 2 1 0] );
hold on;
ypoint = 16.1 - sin(pi-theta2)*r2;
xpoint = position_RIS(1,N1) - cos(pi-theta2)*r2;
label_h = ylabel('Y(m)');
label_h.Position = [-12.12   16.04    1.0000]
caxis([-80 -25])
filename = 'zoomindarkzone';
save = 1;

% [h, wd, ht] = tightfig();
if save == 1
    name1 = append(filename, '.fig');
    name2 = append(filename, '.pdf');
    saveas(gca, name1);
    exportgraphics(gca, name2);
end
% plot(xpoint,ypoint,'r*')
% x = [0 0 0];
% y = [.8 .7 .6];
% labels = {'Telescope Receiver'};
plot(xpoint,ypoint,'o','LineWidth',1)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
imagesc(X,Y,TdB'); 
colormap('hot');
colorbar('vert');
set(gca,'YTick',[0:5:20],'FontSize', 16);
set(gca,'YTickLabel',[20 15 10 5 0],'FontSize', 16);
xlabel('X(m)','FontSize', 16);
ylabel('Y(m)','FontSize', 16,'rotation',0);
% set(gca,'YTick',[0:1:10]);
% set(gca,'YTickLabel',[10 9 8 7 6 5 4 3 2 1 0] );
hold on;
ypoint = 20 - sin(pi-theta2)*r2;
xpoint = position_RIS(1,N1) - cos(pi-theta2)*r2;
label_h = ylabel('Y(m)');
label_h.Position = [-24.1156   13.0000    1.0000];
caxis([-80 -25])
filename = 'tx_energy';
save = 1;

% [h, wd, ht] = tightfig();
if save == 1
    name1 = append(filename, '.fig');
    name2 = append(filename, '.pdf');
    saveas(gca, name1);
    exportgraphics(gca, name2);
end
% plot(xpoint,ypoint,'r*')
% x = [0 0 0];
% y = [.8 .7 .6];
% labels = {'Telescope Receiver'};
plot(xpoint,ypoint,'o','LineWidth',1)
% text(xpoint,ypoint,labels,'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize', 16,'Color', 'white');
% caxis([-30, -100]);
% text(xpoint,ypoint,labels,'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize', 16,'Color', 'white')
%%
% 
% 
figure;
imagesc(x1,y1,ZdB'); 
colormap('hot');
colorbar('vert');
set(gca,'YTick',[0:5:20],'FontSize', 16);
set(gca,'YTickLabel',[20 15 10 5 0],'FontSize', 16);
xlabel('X(m)','FontSize', 16);
ylabel('Y(m)','FontSize', 16,'rotation',0);
% set(gca,'YTick',[0:1:10]);
% set(gca,'YTickLabel',[10 9 8 7 6 5 4 3 2 1 0] );
hold on;
ypoint = 20 - sin(pi-theta2)*r2;
xpoint = position_RIS(1,N1) - cos(pi-theta2)*r2;
label_h = ylabel('Y(m)');
label_h.Position = [-24.1156   13.0000    1.0000];
% caxis([-30, -100]);
% plot(xpoint,ypoint,'r*')
% x = [0 0 0];
% y = [.8 .7 .6];
% labels = {'Telescope Receiver'};
caxis([-80 -25])
plot(xpoint,ypoint,'o','LineWidth',1);
filename = 'ris_engergy';
save = 1;

% [h, wd, ht] = tightfig();
if save == 1
    name1 = append(filename, '.fig');
    name2 = append(filename, '.pdf');
    saveas(gca, name1);
    exportgraphics(gca, name2);
end
% text(xpoint,ypoint,labels,'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize', 16)
%%
% 
% 
% figure;
% % [X1,Y1] = meshgrid(X,Y);
% contour(ZdB'); 
% colormap('hot');
% c = colorbar('vert');
% set(gca,'YTick',[0:5:20],'FontSize', 16);
% set(gca,'YTickLabel',[20 15 10 5 0],'FontSize', 16);
% xlabel('X(m)','FontSize', 16);
% ylabel('Y(m)','FontSize', 16,'rotation',0);
% % set(gca,'YTick',[0:1:10]);
% % set(gca,'YTickLabel',[10 9 8 7 6 5 4 3 2 1 0] );
% hold on;
% ypoint = 20 - sin(pi-theta2)*r2;
% xpoint = position_RIS(1,N1) - cos(pi-theta2)*r2;
% label_h = ylabel('Y(m)');
% label_h.Position = [-24.1156   13.0000    1.0000];
% % caxis([-30, -100]);
% % plot(xpoint,ypoint,'r*')
% % x = [0 0 0];
% % y = [.8 .7 .6];
% % labels = {'Telescope Receiver'};
% caxis([-100 -20])
% plot(xpoint,ypoint,'o','LineWidth',1);
% text(xpoint,ypoint,labels,'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize', 16)
%%
figure;
imagesc(x1,y1,RdB'); 
colormap('hot');
colorbar('vert');
set(gca,'YTick',[0:5:20],'FontSize', 16);
set(gca,'YTickLabel',[20 15 10 5 0],'FontSize', 16);
xlabel('X(m)','FontSize', 16);
ylabel('Y(m)','FontSize', 16,'rotation',0);
% set(gca,'YTick',[0:1:10]);
% set(gca,'YTickLabel',[10 9 8 7 6 5 4 3 2 1 0] );
hold on;
ypoint = 20 - sin(pi-theta2)*r2;
xpoint = position_RIS(1,N1) - cos(pi-theta2)*r2;
label_h = ylabel('Y(m)');
label_h.Position = [-24.1156   13.0000    1.0000];
% caxis([-30, -100]);
% plot(xpoint,ypoint,'r*')
% x = [0 0 0];
% y = [.8 .7 .6];
% labels = {'Telescope Receiver'};
plot(xpoint,ypoint,'go','LineWidth',1.2)
caxis([-80 -25])
filename = 'dark_zone';
save = 0;

% [h, wd, ht] = tightfig();
if save == 1
    name1 = append(filename, '.fig');
    name2 = append(filename, '.pdf');
    saveas(gca, name1);
    exportgraphics(gca, name2);
end
% 'o','LineWidth',1)
% text(xpoint,ypoint,labels,'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize', 16)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    X=-5:0.01:110;
%    Y=-5:0.01:110;
%    L_X=length(X);
%    Z=zeros(1,L_X);
%    
%    [X1,Y1,Z1]=griddata(X,Y,Z,linspace(min(X),max(X),150)',...
%        linspace(min(Y),max(Y),150),'cubic');
%    [X2,Y2,Z2]=griddata(X,Y,Z,linspace(min(X),max(X),150)',...
%        linspace(min(Y),max(Y),150),'cubic');
%    [X3,Y3,Z3]=griddata(X,Y,Z,linspace(min(X),max(X),150)',...
%        linspace(min(Y),max(Y),150),'cubic');
%    Num=length(X3);
% 
%     for i=1:Num
%         for j=1:Num
%             r1=sqrt((X1(i,j)-58.601)^2+(Y1(i,j)-100)^2);
%             r2=sqrt((X2(i,j)-0.868)^2+(Y2(i,j)-0)^2);
%             Z1(i,j)= sqrt(P_Rx) * sin(r1*fc + Phase_rece);
%             Z2(i,j)= N*sqrt(P_Rx_/N/N)* sin(fc*r2 + phase_delay_all);
%             Z3(i,j)=Z1(i,j)+Z2(i,j);
%         end
%     end

   
%    figure
%    surf(X,Y,S+x);
%    colorbar('vert')
%    view(-8,87)

%     figure();
%     surf(X1,Y1,Z1);
%     colorbar('vert')
%     view(-8,87)
%     figure();
%     surf(X2,Y2,Z2);
%     colorbar('vert')
%     view(-8,87)
%     
%     figure
%     surf(X3,Y3,Z3);%
%     colorbar('vert')
%     view(-8,87)
% 
%    figure
%    contourf(X3,Y3,Z3),shading interp;
%    colorbar('vert')
   
   
%    xlabel('time (in seconds)');
%    title('Signal versus Time');
%    zoom xon;

