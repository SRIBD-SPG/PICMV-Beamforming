function [ fval ] = funcpoly2( t, mu, gamma_temp, temp, cPhi, Rho )
%FUNCPOLY Summary of this function goes here
%   Detailed explanation goes here

alpha = abs(temp)./(2*sqrt(t./gamma_temp).*cPhi)-Rho/2;
idx = find(alpha>0);
alpha = alpha.*cPhi.^2./gamma_temp;
if ~isempty(idx)
    fval = mu - sum(alpha(idx));
else
    fval = mu;
end

end

