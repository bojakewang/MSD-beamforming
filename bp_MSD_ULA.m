% Modal Subspace Decomposition code
% "Generalized broadband beamforming using a Modal Subspace Decomposition"
% he, wang  2017-07-15

clear all; 

narr = 10;  
phi_arr = 360/narr;
stm = 32;
M = stm*narr;
cp = 350; % in the air
freq_l = 4e2; freq_h = 4e3; fs = 2*freq_h;

wavel = cp/freq_h;
radius = 0.0372159728654;  % wavel/(4*sin(pi/narr))
% source
phi_s = 0;
theta_s = 90;
Nfft = 300;
f = linspace(1, Nfft, Nfft)/Nfft*fs;
f_ind = find(f>=freq_l & f<=freq_h);
f_s = f(f_ind);

% f_s = linspace(freq_l, freq_h, Nfft);
k_s = 2*pi*f_s/cp;
dk_s = k_s(2) - k_s(1);
phi = linspace(-180, 180, 361);
dphi = 1;
% modal = M; % number of the modal using
[weight, pos, BeamP] = BeamP_ULA(narr, phi_s, theta_s, phi, fs, freq_l, freq_h, Nfft, cp);

%% compute the operator spectrum ...
ZZ = zeros(M, M);
sample_t = 1/fs*repmat(1:stm, 1, narr);
sample_s = kron(1:narr, ones(1, stm));
for ith = 1:M
    for jth = 1:M 
        tmp_t = k_s.*exp(1j*k_s*cp*(sample_t(jth)-sample_t(ith)));
        dist = norm( pos(:, sample_s(ith)) - pos(:, sample_s(jth)) );  %%%%
        tmp_phi = besselj(0, k_s*dist);
        ZZ(ith, jth) = 2*pi*sum( tmp_t.*tmp_phi )*dk_s;
    end
end
[LVect, Eigen, RVect] = eig(ZZ);
Eigen = diag(Eigen);

[Eigen, loc] = sort(Eigen, 'descend');
VectModes = RVect(:, loc);

% figure;
% plot(Eigen, 'b-*', 'Linewidth', 1.2)
% xlabel('number n'); ylabel('eigenvalue \lambda');
% set(gcf, 'position', [0,0,450, 400])
% set(gca, 'Position', [0.12, .11, 0.85, 0.8150])
% return
%% compute continues Modes ...

threshold = 3e3;
loc_hat = find(Eigen>threshold);
modal = numel(loc_hat)

% modal = M

% modal could be determinted by the eigenvalue
opA = zeros(M, length(phi));
ContiModes = zeros(length(k_s), length(phi), modal);
for ith = 1:modal % ----.> n
    for jth = 1:length(k_s)  % ----.> k
       for kth = 1:M  % ----.> m
          opA(kth, :) = exp(1j*k_s(jth)*(cp*sample_t(kth)+radius*cos(sample_s(kth)*phi_arr*pi/180 - phi*pi/180))); 
       end
       tmp = VectModes(:, ith)*ones(1, length(phi));
       ContiModes(jth, :, ith) = 1/sqrt(abs(Eigen(ith)))*sum(opA.*tmp, 1);
    end
    ContiModes(:, :, ith) = ContiModes(:, :, ith)/abs(max(max(ContiModes(:, :, ith))));    
end

%% compute modal coefficients
CModes = ContiModes;
% first kind of desire BeamPattern
Wdes = BeamP(f_ind, :);

% % second kind of desire beampattern
% Wdes = zeros(Nfft, length(phi));
% Wdes(:, 171:190) = 1;
% Wdes(ceil(Nfft/3):end, 171:190) = 1;
% Wdes = Wdes(f_ind, :);

% % third kind of desire beampattern
% Wdes = zeros(Nfft, length(phi));
% Wdes(ceil(Nfft/3.5):ceil(Nfft/3), 101:110) = 1;
% Wdes(ceil(Nfft/4.5):ceil(Nfft/4), 201:210) = 1;
% Wdes = Wdes(f_ind, :);

% figure;
% mesh(phi,f(f_ind),(Wdes));
% shading interp;
% colorbar; colormap('jet');
% set(gca,'XLim',[-180,180],'YLim',[freq_l,freq_h]);
% xlabel('\phi azimuth angle ','FontName','Microsoft YaHei','FontSize',12);
% ylabel('Frequencpy (Hz)','FontName','Microsoft YaHei','FontSize',12);
% zlabel('Beampattern(dB)','FontName','Microsoft YaHei','FontSize',12);
% title('Wdes','FontName','Microsoft YaHei','FontSize',12);

k_ach = k_s;
bn = zeros(modal, 1);
for ith = 1:modal
   Un = CModes(:, :, ith);
   tmp = sum( Wdes.*conj(Un), 2)*dphi*pi/180;
   bn(ith) = sum( tmp.*k_ach.' )*dk_s;
end

%% recompute the BeamPattern, i.e. Wach
tmp = reshape(bn, 1, 1, modal);
bn_rep = repmat(tmp, [length(f_ind), length(phi), 1]); % 
Wach = sum( bn_rep.*CModes, 3 );
Wach = abs(Wach);
Wach = Wach/max(max(Wach));

%% results ...

figure;
plot(sort(Eigen, 'descend'), 'b-*', 'Linewidth', 1.2)
xlabel('number n','FontName','Microsoft YaHei','FontSize',12); 
ylabel('Eigenvalue \lambda','FontName','Microsoft YaHei','FontSize',12);
title('Eigenvalue distribution','FontName','Microsoft YaHei','FontSize',12);
set(gcf, 'position', [500,500,450, 400])
set(gca, 'Position', [0.12, .11, 0.85, 0.8150])

figure;
mesh(phi,f(f_ind),20*log10(Wach));
shading interp;
caxis([-40 0]);
colorbar; colormap('jet');
set(gca,'XLim',[-180,180],'YLim',[freq_l,freq_h]);
xlabel('\phi azimuth angle ','FontName','Microsoft YaHei','FontSize',12);
ylabel('Frequencpy (Hz)','FontName','Microsoft YaHei','FontSize',12);
zlabel('Beampattern (dB)','FontName','Microsoft YaHei','FontSize',12);
title('Wach','FontName','Microsoft YaHei','FontSize',12);
% set(gcf, 'position', [500,500,450, 400])
% set(gca, 'Position', [.13, .11, 0.85, 0.8150])

figure;
mesh(phi,f(f_ind),20*log10(Wdes));
shading interp;
caxis([-40 0]);
colorbar; colormap('jet');
set(gca,'XLim',[-180,180],'YLim',[freq_l,freq_h]);
xlabel('\phi azimuth angle ','FontName','Microsoft YaHei','FontSize',12);
ylabel('Frequencpy (Hz)','FontName','Microsoft YaHei','FontSize',12);
zlabel('Beampattern (dB)','FontName','Microsoft YaHei','FontSize',12);
title('Wdes','FontName','Microsoft YaHei','FontSize',12);
% set(gcf, 'position', [500,500,450, 400])
% set(gca, 'Position', [0.13, .11, 0.85, 0.8150])

%% recompute the weight, i.e. w

