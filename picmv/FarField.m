clc;clear;

% setting
fc = 1e9;
lamb = physconst('LightSpeed')/fc;
N = 30;
PosX = [0:N-1]*lamb/2;

nTheta = 361;
AngleSet = linspace(0,180,nTheta);
SVSet = SteerVec( PosX', 0, lamb, AngleSet, 0 );

Theta0 = 90;

% interference

PhiSet = [87,102];
K = length(PhiSet);
NearSet = -1:0.5:1;

%% static beam
idx = AngleSet==Theta0;
w0 = SVSet(:,idx)/norm(SVSet(:,idx))^2;
Pattern0 = 20*log10(abs(w0'*SVSet));
%% robust

% cal. H_Theta
H_Theta = struct();
idx = find(AngleSet==Theta0);
H_Theta.H = SVSet(:,idx);
H_Theta.c = 1e-1;
Size_Theta = 1;


% cal. H_Theta
H_Phi = struct();
Size_Phi = length(NearSet);
for k = 1:K
    tmp = PhiSet(k) + NearSet;
    for j = 1:Size_Phi
        
        idx = find(AngleSet == tmp(j));
        H_Phi(k).H(:, j) = SVSet(:,idx);
        H_Phi(k).c = 1e0;
    end
end

%% admm for picmv
R = eye(N);
delta = 1e-9;
rho = 1e3;
mu = 1e4;
gamma = ones(K,1);


[ wpimcv, t ] = PICMV_New_ADMMBis( rho, mu, R, gamma, delta, H_Theta, H_Phi, N, K, Size_Theta, Size_Phi  );
% [ wpimcv, t ] = PICMV_ADMMGold( rho, mu, R, gamma, H_Theta, H_Phi, N, K, Size_Theta, Size_Phi);
Pattern = 20*log10(abs(wpimcv'*SVSet));

disp(['objective function', num2str(t)]);



%% figures
figure;
hold on;grid on;
plot(AngleSet,Pattern,'r-', 'linewidth',1);
% plot(AngleSet,Pattern0,'b-', 'linewidth',1);
for k =1:K
    plot([PhiSet(k),PhiSet(k)], [-100,20], 'm--', 'linewidth',2);
end
ylim([-100,5]);











