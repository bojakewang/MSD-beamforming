
function [weight, pos, BeamP] = BeamP_UCA(narr, phi_s, theta_s, phi, radius, fs, Nfft, cp)
% clear all
% 
% % beampattern for UCA
% 
% cp = 350;
% narr = 5;
% freq_h = 4e3;
% freq_l = 4e2;
% fs = 8e3;
% Nfft = 500;
% f = linspace(1, Nfft, Nfft)*fs/Nfft;
% 
% radius = 0.0372159728654;  % wavelength/(4*sin(pi/narr))
% phi_s = 0;
% theta_s = 90;
% 
% phi = linspace(-180, 180, 361);

phi_narr = 360/narr;
pos = zeros(2,narr);
for ii = 1:narr
    pos(1, ii) = radius*cos((phi_narr/2+phi_narr*ii)*pi/180); % for x position
    pos(2, ii) = radius*sin((phi_narr/2+phi_narr*ii)*pi/180); % for y position
end


w_CBF = zeros(narr,Nfft);  
for ii = 1:Nfft
   f = floor(fs/Nfft)*(ii);
   a0 = [-cos((pi/180)*phi_s);-sin((pi/180)*phi_s)];
   tao0 = (a0.'*pos)/cp; 
   w_CBF(:,ii) = exp(1j*2*pi*(tao0.')*f)/narr; 
end

%% 波束模式图
f = (1:Nfft)*(fs/Nfft);
Beampattern = zeros(Nfft,length(phi));
for ii = 1:length(phi)
    a = [-cos((pi/180)*phi(ii));-sin((pi/180)*phi(ii))];
    tao = (a.'*pos)/cp;
    v = exp(-1j*2*pi*(tao.')*f);  % N*(length(f))
    Beampattern_tmp = w_CBF.'*v;
    Beampattern(:,ii) = diag(Beampattern_tmp);
end
 Beampattern = abs(Beampattern);  %length(f)*length(phi)
 Beampattern = Beampattern/max(max(Beampattern));
 
 weight = w_CBF; 
 BeamP = Beampattern;

 
%  Beampattern = 20*log10( Beampattern );
% % 波束模式图
% figure;
% mesh(phi,f,Beampattern);
% shading interp;
% % cpaxis([-40 0]);
% colorbar; colormap('jet');
% set(gca,'XLim',[-180,180],'YLim',[freq_l,freq_h],'ZLim',[-60,0]);
% xlabel('DOA angle \theta','FontName','Times New Roman','FontSize',10);
% ylabel('Frequencpy (Hz)','FontName','Times New Roman','FontSize',10);
% zlabel('Beampattern(dB)','FontName','Times New Roman','FontSize',10);
% title('常规宽带波束形成');