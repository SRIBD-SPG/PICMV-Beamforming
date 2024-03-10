function [w_opt, y_opt] = CheckSubProblemywBis(b, V, D, alpha, beta)

% M = 10;
% 
% b = 1e0*(randn(M,1) + sqrt(-1)*randn(M,1));
% alpha = 1e0*abs(randn());
% beta = 1e0*(randn());
% A = 1e0*(randn(M,M) +  sqrt(-1)*randn(M,M)); % upper bound
% A = A + A';
% [V, D] = eig(A);
% D = abs(D);



%%
M = length(b);

epsy = 1e-5;
A = V*D*V';
A = (A + A')/2;
%% 
bb = V'*b;
dd = 2*diag(D);
wy = @(y) -V*diag(1./(dd+2*alpha+beta/y))*bb;
fy = @(y)norm( ( bb./( (dd+2*alpha)*y+beta ) ) );

wbar = -A\b/2;
ybar = -beta/2/alpha;
if norm(b) <= beta
    y_opt = 0;
    w_opt = zeros(M,1);
elseif norm(wbar)<ybar
%     disp('this')
    y_opt = ybar;
    w_opt = wy(y_opt);
else
    ymax = norm(wbar);
    ymin = max(ybar,0);

    while abs(ymin-ymax)>epsy
        ymid = (ymin+ymax)/2;
        tmp = fy(ymid);
        if tmp<1
            ymax = ymid;
        elseif tmp>1
            ymin = ymid;
        else
            break
        end
    end
    y_opt = (ymin+ymax)/2;
    w_opt = wy(y_opt);
end

%%



%% bisection for y

% f = @(y)norm( (V*diag(1./(diag(D)+(alpha + beta/(2*y))))*V')*b )/2/y;
% wbar = -A\b/2;
% ybar = - beta/2/alpha;
% 
% 
% if norm(wbar)<= ybar
%     w_opt = wbar;
%     y_opt = ybar;
% else
%     
%     ymax = norm(wbar);
%     ymin = max(ybar,0);
%     
% %     ystep = linspace(ymin,ymax,100);
% %     for i = 1:length(ystep)
% %         fstep(i) = f(ystep(i));
% %     end
% %     figure;
% %     plot(ystep,fstep,'b-');
% 
%     while abs(ymin-ymax)>epsy
%         ymid = (ymin+ymax)/2;
%         tmp = f(ymid);
% %         fy(ymid)
%         if tmp<1
%             ymax = ymid;
%         elseif tmp>1
%             ymin = ymid;
%         else
%             break
%         end
%     end
% 
%     y_opt = (ymin+ymax)/2;
%     w_opt = -(V*diag(1./(diag(D)+(alpha + beta/(2*y_opt))))*V')*b /2;
% 
% end

% 
% y_opt, y_opt1
% w_opt,w_opt1


% cvx_begin
% cvx_precision high
% variable w(M) complex
% variable y
% minimize( w'*A*w + real(b'*w) + alpha*y*y + beta*y)
% subject to 
% norm(w) <= y
% cvx_end
% 
% w,w_opt
% y,y_opt
end