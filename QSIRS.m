function [q,x] = QSIRS(N, beta, gamma, mu, epsilon)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

% X(1): Susceptible
% X(2): Infected


% parameters
% beta - transmission within household 
% gamma - recovery of each infected individual

% mu = beta*S*I/(N-1); % infection- logistic growth beta*I*(N-I)/(N-1)
% lambda = gamma*I;       % recovery

n = sum(1:N+1); % number of states ordered (0,3),(0,2),(0,1), (0,0),...,(1,2),(1,1),..

q = sparse(n,n);


% Assigning coordinates to states
 
ss = [];

for i = 0:N
    ss = [ss repelem(i,N+1-i)];
end

ii = [];

for i = N:-1:0
    ii = [ii i:-1:0];
end

x = [ss;ii]; % [S, I]

% Entering recovery and infection rates into Q

for i = 1:n
    if ii(i) > 0
        search = x == [ss(i); ii(i)-1];
        pos = sum(search) == 2;
        q(i,find(pos)) = gamma*ii(i);
    end
end

for i=1:n
    if ss(i) > 0 && ii(i) > 0 % checking that there is someone to be infected and someone to spread the disease
        search = x == [ss(i)-1; ii(i)+1];
        pos = sum(search)==2;
        q(i, find(pos)) = beta*ss(i)*ii(i)/(N-1) + ss(i)*epsilon;
    end
end


% now entering waning immunity

for i=1:n
    if N - ss(i) - ii(i) > 0  % checking that someone has recovered
        search = x == [ss(i)+1; ii(i)];
        pos = find(search(1,:).*search(2,:));
        q(i, pos) = mu*(N-sum(x(:,i)));
    end
end

% external infection 
for i=1:n
    if ss(i) > 0 && ii(i) >= 0 % checking that there is someone to be infected and someone to spread the disease
        search = x == [ss(i)-1; ii(i)+1];
        pos = sum(search)==2;
        q(i, find(pos)) = q(i, find(pos)) + ss(i)*epsilon;
    end
end


for i = 1:n
    q(i,i) = -sum(q(i,:));
end

end




