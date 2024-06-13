clear
% n household of size 1 --------------------------------------------------

clear

beta = 0; % no internal internal infection
N = 1; % household of size 1
houses = 10; % number of initial houses
maxhouses = 30;
gamma = 1; % recovery
mu = 0.5; % waning immunity
alpha = 1; % external infection
tmax = 20;

% 1 person in each household
[Q, x] = Qgen(houses, alpha, gamma); % Q matrix for 10 houses of size 1, 
% note that beta is now gamma as gamma is extenal infection 

i_0 = find(x(2,:)==0); % finding the absorbing states
ss = x(1,:); % number susceptible for each state
ii = x(2,:); % number infected for each state

ind_c = find(ii>0); % index for C - transient states
n = sum(1:houses+1); % number of states

Q_c = Q(ind_c,ind_c); % transient states only
f = ii(ind_c)'; % f is the number of infected

% solving yi(s)-----------------------------------------------------------

s = 1;
y = ones(n,1);

Q_f = (Q_c-alpha*diag(f)); % solving from proposition 1
v = -gamma*(f==1); 
y_c = Q_f\v; % y_c as were only considering y_i for i in c
g = y_c; % as m=0 g(0) = y(1)

y(ind_c) = y_c;

% solving the derivative of proposition 1
for m=2:(maxhouses+1)
    y_c(:,m) = (Q_f)\((m-1)*(alpha*f).*y_c(:,m-1));
    % m-1 as our first col is for 0, m represents the mth derivative of y
    g(:,m) = y_c(:,m)*(((-1)^(m-1))/factorial(m-1));
end

ic = find(ii==1 & ss == N-1);

ic_c = find(ind_c==ic); % finding the position in C which corresponds to intial conditions in S
figure(2)
plot(0:maxhouses,g(ic_c,:))
title('PDF of the number of houses infected from 1 infected individual, $\textit{g(m)}$','Interpreter','latex')
ylabel('Probability','Interpreter','latex')
xlabel('Number of households, $m$','Interpreter','latex')

