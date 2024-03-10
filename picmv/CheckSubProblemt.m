function [z,y,t_opt] = CheckSubProblemt(a, b, alpha, beta, delta, c, gamma, mu, N)

% K = 2;
% S = 3;
% mu = 1e-3;
% 
% gamma = 1e0*abs(randn(K,1));
% 
% a = 1e0*abs(randn(K*S,1));
% b = 1e0*(randn(K*S,1) + sqrt(-1)*randn(K*S,1));
% 
% alpha = 1e0*abs(randn(K*S,1));
% beta = 1e0*(randn(K*S,1));
% 
% delta = 1e-1*abs(randn());
% 
% c = abs(randn(K*S,1));
% y = min(c/delta,randn());


%% closed form for t
% convenient parameter
c = c./gamma;
delta2 = delta.^2;
vb = abs(b);
alpha2c = (2.*alpha.*c);
betad = beta.*delta;
betac = beta.*c;
ac = a.*c;
dinom = alpha + a.*delta2;

t0 = -(betad + delta2.*vb)./alpha2c;
t1 = - (betad - alpha./a.*vb)./alpha2c;

% sort_t = sort([t0;t1])';

a1 = alpha2c.*c./delta2;
b1 = betac./delta;
% g1 = a1(1)*sort_t;
% g1 = a1*sort_t + b1*ones(size(sort_t));

a2 = alpha2c.*ac./dinom;
b2 = (ac.*betad - alpha2c.*vb/2)./dinom;
% g2 = a2*sort_t + b2*ones(size(sort_t));
% g2 = a2(1)*sort_t;

% tb1 = (b2-b1)/(a1(1)-a2(1));
% tb2 = -b2/a2(1);
% g = zeros(1,2*N);

a1c = a1(1);
a2c = a2(1);

%% bisection search for t
epst = 1e-5;
if ft(min(t0),a1c,a2c,b1,b2)>-mu
    t_opt = -(sum(b1)+mu)/(sum(a1));
else
    tmax = max(t1);
    tmin = min(t0);

    while abs(tmin-tmax)>epst
        tmid = (tmin+tmax)/2;
        tmp = ft(tmid,a1c,a2c,b1,b2);
        if tmp<-mu
            tmin = tmid;
        elseif tmp>-mu
            tmax = tmid;
        else
            break
        end
    end
    t_opt = (tmin+tmax)/2;
end


%% closed-form for t
% 
% for i = 1:N
%     gtmp= (a1c*sort_t + b1(i)) .* (sort_t<=tb1(i)) + (a2c*sort_t + b2(i)) .* (sort_t>tb1(i) & sort_t<=tb2(i)) + zeros(1,2*N);
%     g = g + gtmp;
% %     tmp = [g1+b1(i);g2 + b2(i);zeros(1,2*N)];
% %     g = g + min(tmp,[],1);
% end
% 
% 
% idx = find(g <= -mu);
% 
% if isempty(idx)
%     t_opt = -(sum(b1)+mu)/(sum(a1));
% else
%     idx = idx(end);
%     idx1 = find(t0>=sort_t(idx+1));
%     idx2 = find(t1>=sort_t(idx+1) & t0 <=sort_t(idx));
%     
%     sa = sum(a1(idx1)) + sum(a2(idx2));
%     sb = sum(b1(idx1)) + sum(b2(idx2));
%     
%     t_opt = -(sb + mu)/(sa);
%     
% end

% figure;
% hold on;
% plot(sort_t, g, 'b-o');
% plot([t_opt,t_opt],[min(g),max(g)],'r--');
% plot([min(t),max(t)],[-mu,-mu],'b--');
%% closed-form for fixed t
tau = b./(2*a);
ep = exp(sqrt(-1)*angle(b));
abar = 2.*a.*delta2;
bbar = 2.*ac.*delta.*t_opt - delta.*vb;

y = -beta./(2*alpha);
z = -tau;

tmp = 2*alpha.*(bbar./abar) + beta + abar.*(bbar./abar) - bbar;
idx = find(tmp<0);
y(idx) = (bbar(idx) - beta(idx)) ./ ( 2*alpha(idx) + abar(idx));
z(idx) = -ep(idx).*(c(idx)*t_opt-delta*y(idx));

tmp = 2*alpha.*(c./delta*t_opt) + beta + abar.*(c./delta*t_opt) - bbar;
idx = find(tmp<0);
y(idx) = c(idx)/delta*t_opt;
z(idx) = 0;


    function v = ft(t,a1,a2,b1,b2)
        v = 0;
        for cnt = 1:length(b1)
            v = v+min([a1*t+b1(cnt),a2*t+b2(cnt),0]);
        end
        
    end


% Az = diag(a);
% Ay = diag(alpha);
% cvx_begin
% cvx_precision high
% variable z(N) complex
% variable y(N)
% variable t
% 
% minimize( mu*t + z'*Az*z + real(b'*z) + y'*Ay*y + beta'*y)
% subject to 
% for i = 1:N
%     abs(z(i))<= c(i)*t - delta*y(i)
% end
% cvx_end
% 
% t_opt, t
% z_opt, z
% y_opt,y

end
