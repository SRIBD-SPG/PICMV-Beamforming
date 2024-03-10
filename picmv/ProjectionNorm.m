function [ w ] = ProjectionNorm( w, value )
%PROJECTION Summary of this function goes here
%   Detailed explanation goes here
w = w/max(abs(w));
for i = 1:length(w)
    if abs(w(i)) ~= 0
    tmp = abs(abs(w(i)) - value);
    [~,idx] = min(tmp);
    w(i) = w(i)/abs(w(i))*value(idx);
    end
end

end

