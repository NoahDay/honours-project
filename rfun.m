function [y] = rfun(r,Q_c,f,ic_c,ind_c)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
z = (Q_c-r*diag(ones(1,length(ind_c))))\f;
y = z(ic_c) - 1;
end

