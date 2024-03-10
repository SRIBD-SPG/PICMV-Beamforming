function [ w, t ] = PICMV_New_ADMMBis( rho, mu, R, gamma, delta, H_Theta, H_Phi, M, K, Size_Theta, Size_Phi  )
%PICMV_NEW_ADMM Summary of this function goes here
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

eps = 1e-4; % numerical precision
MaxIteration = 200; % maximun number of iteration
% set some convenient parameters
R = (R+R')/2;

% Convenient Para. 
HTheta = H_Theta.H;
cTheta = H_Theta.c;
HPhi = zeros(M,K*Size_Phi);
cPhi = zeros(K*Size_Phi,1);
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
alp_y = (Size_Theta + K*Size_Phi)*rho/2;
A = R+ rho/2*[HTheta,HPhi]*[HTheta,HPhi]';
A = (A+A')/2;
[V, D] = eig(A);
if min(diag(D))<=0
%     A = A + abs(min(diag(D)))*eye(M);
    D = D + 1.001*abs(min(diag(D)))*eye(M);
end
A = V*D*V';

% [~,Dla] = eig(blkdiag(A,alp_y));
% la = abs(max(diag(Dla)));
% if min(diag(D))<=0
%     A = A + abs(min(diag(D)))*eye(M);
%     la = la + abs(min(diag(D)));
% end
% A = la*eye(M);
% initial variables and multipliers
w = zeros(M,1);
% epsilon = zeros(K,1);
z_theta = zeros(Size_Theta,1);
y_theta = zeros(Size_Theta,1);
z_phi = zeros(K*Size_Phi,1);
y_phi = zeros(K*Size_Phi,1);
y = 0;
t = 0;

lamb_theta = ones(Size_Theta,1);
lamb_phi = ones(K*Size_Phi,1);
eta_theta = ones(Size_Theta,1);
eta_phi = ones(K*Size_Phi,1);

%% iteration starts
no_ite = 1; % count numer of iteration
er = 1; % feability gap, as iteration goes, er goes to zero
% Ainv = inv(A);
% for ite = 1:no_ite
while er>eps & no_ite<MaxIteration
    no_ite = no_ite +1
    z_all = [z_theta;z_phi];
%     Epsilon = epsilon;
    wbar = w;
    ybar = y;
    tbar = t;
    %% the first block
    % update w
    if no_ite > 1
        b = HTheta*(conj(lamb_theta)-rho*conj(z_theta))+...
            HPhi*(conj(lamb_phi)-rho*conj(z_phi));
        %     b = b + 2*A*w - 2*la*w;
        
        %     eta_theta
        be = -(sum(y_phi) + sum(y_theta))*rho  + sum(eta_phi) + sum(eta_theta);
        %     be = be + 2*alp_y*y - 2*la*y;
        %     [w, y] = CheckSubProblemyw(la, b, A, la, be);
%         tic
        [w, y] = CheckSubProblemywBis(b, V, D, alp_y, be);
%         toc
    else
        b = HTheta*(conj(lamb_theta)-rho*conj(z_theta))+...
            HPhi*(conj(lamb_phi)-rho*conj(z_phi));
        w = -A\b/2;
    end
%     b = b/2;
%     w =  -Ainv*b;
    
    %% the second block
    % update delta for Theta
    b = -lamb_theta - rho*(w'*HTheta).';
%     a = rho/2;
    be = - eta_theta - rho*y;
%     tic
    for i = 1:Size_Theta
        [z_theta(i), y_theta(i)] = CheckSubProblemTheta(rho/2, b(i), rho/2, be(i), cTheta(i), delta, 1);
    end
%     toc
    %% the third block
    % update for delta_Phi, t
    
    gamma_temp = kron(gamma,ones(Size_Phi,1));
    b = - lamb_phi - rho*(w'*HPhi).';
    be = - eta_phi - rho*y;
    a = rho/2*ones(size(b));
    alp = rho/2*ones(size(b));
%     tic
    [z_phi, y_phi, t] = CheckSubProblemt(a, b, alp, be, delta, cPhi, gamma_temp, mu, K*Size_Phi);
%     toc
    
    % compute epsilon
%     temp = lambda_phi + Rho*HPhi'*w;
    % updata epsilon
%     for k = 1:K
%         abs_temp = abs(temp((k-1)*Size_Phi+1:k*Size_Phi))./cPhi((k-1)*Size_Phi+1:k*Size_Phi);
%         epsilon(k) = max(abs_temp.^2/rho^2);
%     end
%     idx = find(abs(temp)>rho*sqrt(t./gamma_temp).*cPhi);
%     if ~isempty(idx)
%         for i = 1:length(idx)
%             idx_k = ceil(idx(i)/Size_Phi);
%             epsilon(idx_k) = t/gamma(idx_k);
%         end
%     end
   
    
    %% multipliers
    lamb_theta = lamb_theta + rho*((w'*HTheta).' - z_theta);
    lamb_phi = lamb_phi + rho*((w'*HPhi).' - z_phi);
    
    eta_theta = eta_theta + rho*(y - y_theta);
    eta_phi = eta_phi + rho*(y - y_phi);
    
    %% primal and dual feasibility
    % primal feasibility
    er = max(abs((w'*HPhi).' - z_phi));
    er = max([er;abs((w'*HTheta).' - z_theta)]);
    % primal iteration
    er = max([er;abs([z_theta;z_phi]-z_all)]);
%     er = max([er;abs(Epsilon-epsilon)]);
    er = max([er;abs(w-wbar)]);
    er = max([er;abs(t - tbar)]);
    er = max([er;abs(y - ybar)]);
    % dual iteration
%     er = max([er;rho*abs((w'*HPhi).' - z_phi)]);
%     er = max([er;rho*abs((w'*HTheta).' - z_theta)]);
    er
end
disp('P-ICMV bisection: ');
no_ite
er

% save activity of constraints
% 
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

