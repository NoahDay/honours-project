clear
% n household of size 1 --------------------------------------------------

clear

beta = 6; % internal infection
N = 3; % size of household
houses = 10;
maxhouses=30;
gamma = 1; % recovery
mu = 0.5; % waning immunity
alpha = 1; % external infection
tmax = 20;

% 1 person in each household
[Q, x] = Qgen(houses, alpha, gamma); % Q matrix for 10 houses of size 1, 
% note that beta is now gamma as gamma is extenal infection 

n = sum(1:N+1);
i_0 = find(x(2,:)==0); % finding the absorbing states
ss = x(1,:); % number susceptible for each state
ii = x(2,:); % number infected for each state

ind_c = find(ii>0); %index for i in c

n = sum(1:houses+1);
i_0 = find(x(2,:)==0); % finding the absorbing states

Q_c = Q(ind_c,ind_c); % transient states only

f = ii(ind_c)';


% % for calculating R* we set f = alpha
% f = alpha*ones(3,1);
% 
%  e_temp = -Q\f; 
%  
%  e = zeros(3,1); 
%  
%  pos = find(x(2,:)~=0);
%  
%  for i = 1:length(e_temp)
%      e(pos(i)) = e_temp(i);
%  end
%  
%  
%  e = -Q\f; 
 

 % solving yi(s)-----------------------------------------------------------
 
 s = 1;
 y = ones(n,1);
 

Qf = (Q_c-alpha*diag(f))
v=-gamma*(f==1);
y_c = Qf\v; % y_c as were only considering y_i for i in c
g = y_c; % as m=0 g(0) = y(1)

for i = 1:length(y_c)
    y(ind_c(i)) = y_c(i);
end

for m=2:(maxhouses+1)
    y_c(:,m) = (Qf)\((m-1)*(alpha*f).*y_c(:,m-1));%m-1 as our first col is for 0, m represents the mth derivative of y
    g(:,m) = y_c(:,m)*(((-1)^(m-1))/factorial(m-1));
end



initial_i = find(ii==1 & ss == N-1);

initial_c = find(ind_c==initial_i); % finding the position in C which corresponds to intial conditions in S
figure(2)
plot(0:maxhouses,g(initial_c,:))
title('PDF of the number of houses infected from 1 infected individual, $\textit{g(m)}$','Interpreter','latex')
ylabel('Probability','Interpreter','latex')
xlabel('Number of households, $m$','Interpreter','latex')

