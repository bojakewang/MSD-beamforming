% BeamP_ULA
function [weight, pos, BeamP] = BeamP_ULA(narr, phi_s, theta_s, phi, fs, fmin, fmax, Nfft, cp)

% parameter

f = linspace(fmin,fmax,Nfft); %确定工作的范围
lamda = cp/fmax;
%position of the sensors
px = zeros(1,narr);
py = zeros(1,narr);

for ii = 1:narr
    py(ii) = lamda*(ii-1)/2;
end

% pos保存了传感器的坐标
pos(1,:) = px;  
pos(2,:) = py;

% 信号方向的单位向量，二维
u0 = [cos((pi/180)*phi_s);sin((pi/180)*phi_s)];
tao0 = (u0.'*pos)/cp; %1*narr
% 常规波束形成权矢量
weight = exp(1j*2*pi*(tao0.')*f);  %narr*(length(f));

%% Beampattern
Beampattern = zeros(length(f),length(phi));

 for ii =  1:length(phi)
     u = [cos((pi/180)*phi(ii));sin((pi/180)*phi(ii))];
     tao = (u.'*pos)/cp;
     v = exp(1j*2*pi*(tao.')*f);  %narr*(length(f));
     Beampattern_tmp = weight'*v;
     Beampattern(:,ii) = diag(Beampattern_tmp);
 end
 Beampattern = abs(Beampattern);  %length(f)*length(phi)
 BeamP = Beampattern/max(max(Beampattern));


% 
% figure;
% surf(phi,f,20*log10( Beampattern ));
% shading interp;
% caxis([-80 0]);
% colorbar; colormap('jet');
% set(gca,'XLim',[-90,90],'YLim',[fmin,fmax],'ZLim',[-80,0]);
% xlabel('DOA angle \theta','FontName','Times New Roman','FontSize',10);
% ylabel('Frequency (Hz)','FontName','Times New Roman','FontSize',10);
% zlabel('Beampattern(dB)','FontName','Times New Roman','FontSize',10);