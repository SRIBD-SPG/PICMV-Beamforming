function [ vec ] = SteerVec( x, z, lambda, theta, phi )
%STEERVEC Summary of this function goes here
%   Detailed explanation goes here
%   compute steering vector
%   x,y --  position of all elements on x-y plane
%   lambda -- wave length
%   theta,phi --  theta is azimuth angle, phi is elevation angle
% Uncertainty: U,Phi

vec = exp(-sqrt(-1)*2*pi/lambda*(x*(cosd(theta)*cosd(phi))+z*(ones(size(theta))*sind(phi)))); 
end

