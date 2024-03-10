function [ w, t, delta_theta, delta_phi,  epsilon, ActiveSourceCons, ActiveInterfCons] = PICMV_ADMMGold( Rho, mu, R, gamma, H_Theta, H_Phi, M, K, Size_Theta, Size_Phi )
%ADMMGOLD Summary of this function goes here
%   Detailed explanation goes here
%%   input parameter:
%   Rho -- penalty parameter of ADMM 
%   mu -- trade-off parameter of PICMV
%   R -- correlation matrix
%   gamma -- penalty parameter of PICMV, make difference between intf.
%   H_Theta -- parameters of target source protection constraints
%   H_Phi -- parameters of intrf. rejection constraints
%   M -- total number of microphones
%   K -- number of interfering sources
%   Size_Theta -- number of constraints for target source, e.g., 1 or 3
%   Size_Phi -- number of constraints for each interfering source, e.g., 1
%   or 3

%% output parameters:
%   t -- maximum value of \gamma_k\epsilon_k
%   w -- beamformer
%   delta_theta,delta_phi -- auxillary variables
%   epsinlon -- opt. variable in intf. rejection constraints
%   ActiveSourceCons -- activity for target source protection constraints
%   ActiveInterfCons -- activity for interference rejection constraints

eps = 1e-5; % numerical precision
MaxIteration = 1000; % maximun number of iteration
% set some convenient parameters
R = (R+R')/2;

% D
HTheta = H_Theta.H;
HPhi = zeros(M,K*Size_Phi);
for k = 1:K
    HPhi(:,(k-1)*Size_Phi+1:k*Size_Phi) = H_Phi(k).H;
    cPhi ((k-1)*Size_Phi+1:k*Size_Phi,1) = H_Phi(k).c;
end
% HAll = [HTheta,HPhi]*[HTheta,HPhi]';
% [~,D] = eig(R);
% MaxR = max(diag(D));
% [~,D] = eig(HAll);
% MaxH = max(abs(diag(D)));
% Rho = 1e3*MaxR;

A = R+ Rho/2*[HTheta,HPhi]*[HTheta,HPhi]';
A = (A+A')/2;
% [~,D] = eig(A);
% MaxEv = abs(max(diag(D)));
% if min(diag(D))<=0
%     A = A + 2*abs(min(diag(D)))*eye(2*M);
% end
% if flag_reg
% epsl = 1e-3;
% A = A + epsl*MaxEv*eye(M,M);
% end

% initial variables and multipliers
w = zeros(M,1);
epsilon = zeros(K,1);
delta_theta = ones(Size_Theta,1);
delta_phi = ones(K*Size_Phi,1);
lambda_theta = ones(Size_Theta,1);
lambda_phi = ones(K*Size_Phi,1);
ActiveSourceCons = zeros(Size_Theta,1);
ActiveInterfCons = zeros(K*Size_Phi,1);
%% iteration starts
no_ite = 1; % count numer of iteration
er = 1; % feability gap, as iteration goes, er goes to zero
Ainv = inv(A);
% for ite = 1:no_ite
while er>eps& no_ite<MaxIteration
    no_ite = no_ite +1;
    Delta = [delta_theta;delta_phi];
    Epsilon = epsilon;
    wbar = w;
    %% the first block
    % update w
    b = HTheta*(lambda_theta-Rho*delta_theta)+...
        HPhi*(lambda_phi-Rho*delta_phi);
    b = b/2;
    w =  -A\b;
    
    %% the second block
    % update delta for Theta
    temp = lambda_theta + Rho*HTheta'*w;
    delta_theta = temp/Rho;
    idx = find(abs(temp-Rho)>Rho*H_Theta.c);
    if ~isempty(idx)
        delta_theta(idx) = 1 + H_Theta.c(idx).*(temp(idx)-Rho)./abs(temp(idx)-Rho);
    end
    
    % update for delta_Phi, t
    gamma_temp = kron(gamma,ones(Size_Phi,1));
    temp = lambda_phi + Rho*HPhi'*w;
    t_min = 0;
    t_max = max(gamma_temp.*abs(temp/Rho).^2./cPhi.^2);
    a = t_min; b = t_max+1;
    while abs(b-a)>eps/10
        mid = (a+b)/2;
        fmid = funcpoly2( mid, mu, gamma_temp,temp, cPhi, Rho );
        if fmid<0
            a = mid;
        elseif fmid > 0
            b = mid;
        else
            a = mid;
            b = mid;
        end
    end
    t = (a+b)/2;
    
    % updata delta_phi
    delta_phi = temp/Rho;
    idx = find(abs(delta_phi)>sqrt(t./gamma_temp).*cPhi);
    if ~isempty(idx)
        delta_phi(idx) = temp(idx)./abs(temp(idx)).*sqrt(t./gamma_temp(idx)).*cPhi(idx);
    end
    
    % compute epsilon
%     temp = lambda_phi + Rho*HPhi'*w;
    % updata epsilon
    for k = 1:K
        abs_temp = abs(temp((k-1)*Size_Phi+1:k*Size_Phi))./cPhi((k-1)*Size_Phi+1:k*Size_Phi);
        epsilon(k) = max(abs_temp.^2/Rho^2);
    end
    idx = find(abs(temp)>Rho*sqrt(t./gamma_temp).*cPhi);
    if ~isempty(idx)
        for i = 1:length(idx)
            idx_k = ceil(idx(i)/Size_Phi);
            epsilon(idx_k) = t/gamma(idx_k);
        end
    end
   
    
    %% multipliers
    lambda_theta = lambda_theta + Rho*(HTheta'*w - delta_theta);
    lambda_phi = lambda_phi + Rho*(HPhi'*w - delta_phi);
    
    %% primal and dual feasibility
    % primal feasibility
    er = max(abs(HPhi'*w - delta_phi));
    er = max([er;abs(HTheta'*w - delta_theta)]);
     % primal iteration
    er = max([er;abs([delta_theta;delta_phi]-Delta)]);
    er = max([er;abs(Epsilon-epsilon)]);
    er = max([er;abs(w-wbar)]);
 % dual iteration
    er = max([er;Rho*abs(HPhi'*w - delta_phi)]);
    er = max([er;Rho*abs(HTheta'*w - delta_theta)]);
end
disp('P-ICMV: ');
no_ite
er

% save activity of constraints

% idx = find(abs(delta_theta-1)<H_Theta.c*0.99);
% if ~isempty(idx)
% ActiveSourceCons(idx,1) = 0;
% end
% 
% temp_ep = repmat(epsilon',Size_Phi,1);
% idx = find(abs(delta_phi)<sqrt(vec(temp_ep)).*cPhi*0.99);
% if ~isempty(idx)
% ActiveInterfCons(idx,1) = 0;
% end


end


