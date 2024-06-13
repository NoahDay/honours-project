function [y] = rfunhet(r,Q_c,f,ind_c,pi_c)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
z = (Q_c-r*diag(ones(1,length(ind_c))))\f;


h = [0.244, 0.334, 0.162, 0.159, 0.067, 0.034];
k = 1:6;
pi = h.*k./sum(k.*h); % household size distribution, probability that a 
% randomly selected person from the population belongs to a household

y = sum(pi_c'.*z) - 1;
end

