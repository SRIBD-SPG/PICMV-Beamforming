function [z, y] = CheckSubProblemTheta(a, b, alpha, beta, c, delta, d)


% a = abs(randn());
% b = randn() + sqrt(-1)*randn();
% 
% alpha = abs(randn());
% beta = 1e0*randn();
% 
% d = randn() + sqrt(-1)*randn();
% c = abs(randn());
% delta = 1e-1;

% y = min(c/delta,randn());


%% closed-form
tau = b/(2*a);
r = abs(tau + d);
ep = exp(sqrt(-1)*angle(tau + d));
% z_opt = d + ep*delta*max((-r+c)/delta, y) - c*ep;

if -beta/(2*alpha)<=(c-r)/delta
%     disp('left side')
    y = -beta/(2*alpha);
    z = d + ep*delta*max((-r+c)/delta, y) - c*ep;
%     zp = d - ep*min(r,c-delta*y);
else
    tmp = (2*a*delta*(c-r) - beta)/(2*a*delta^2+2*alpha);
    if tmp<=c/delta
%         disp('inner')
        y = tmp;
        z= d + ep*delta*max((-r+c)/delta, y) - c*ep;
%         zp = d - ep*min(r,c-delta*y);
    else
%         disp('right side')
        y = c/delta;
        z = d + ep*delta*max((-r+c)/delta, y) - c*ep;
%         zp = d - ep*min(r,c-delta*y);
    end
end

% z,zp

% 
% cvx_begin
% cvx_precision high
% variable z complex;
% variable y;
% 
% minimize( a*(z'*z) + real(b'*z) + alpha*y*y + beta*y)
% subject to 
% abs(z-d)<= c -delta*y
% 
% cvx_end
% 
% z_opt, z
% y_opt, y


end