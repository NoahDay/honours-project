function [q,x] = QSIS(N, beta, gamma, eps)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

% X(1): Susceptible
% X(2): Infected


% parameters
% beta - transmission within household 
% gamma - recovery of each infected individual

% mu = beta*S*I/(N-1); % infection- logistic growth beta*I*(N-I)/(N-1)
% lambda = gamma*I;       % recovery

% number of states ordered (0,3),(1,2),(2,1),(3,0)



% Assigning coordinates to states
 
ss = 0:N;

ii = N:-1:0;


x = [ss;ii]; % [S, I]

q = sparse(length(ss),length(ss));
% Entering recovery and infection rates into Q

for i=1:N
    if ss(i) > 0 && ii(i) % checking that there is someone to be infected and someone to spread the disease
        search = x == [ss(i)-1; ii(i)+1];
        pos = sum(search)==2;
        q(i, find(pos)) = beta*ss(i)*ii(i)/(N-1); %+eps*ii(i);
    end
end

for i = 1:N
    if ii(i) > 0
        search = x == [ss(i)+1; ii(i)-1];
        pos = sum(search) == 2;
        q(i,find(pos)) = gamma*ii(i);
    end
end




for i = 1:N
    q(i,i) = -sum(q(i,:));
end


end




