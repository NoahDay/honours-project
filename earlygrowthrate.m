% abs.gov.au Aus ave hhsize = 2.6, using parameters based on basic
% parameter Ross 
clear

beta = 6; % internal infection
N = 2; % size of household
houses = 2;
gamma = 1; % recovery
alpha = 1; % external infection
[Q, x] = Qgen(houses, beta, gamma);

ss = x(1,:); % number susceptible for each state
ii = x(2,:); % number infected for each state

n = sum(1:N+1);
ind_c = find(ii>0); % index for C - transient states
Q_c = full(Q(ind_c,ind_c)); % transient states only

f = ii(ind_c)';

% solving e(i) = E[Gamma|X(0) = i] ----------------------------------------
e_temp = -Q_c\f; 
e = zeros(n,1); 
e(ind_c) = e_temp;

ic = find(ii==1 & ss == N-1); % intial conditions
R_star = e(ic);


% finding r ---------------------------------------------------------------
% using proposition 1 we can solve this expectation of a path integral like
% we did for y

ic_c = find(ind_c==ic); 
% finding the position in C which corresponds to intial conditions in S

% only search over possible positive values of early growth
rmin = 0;
rmax = gamma*(R_star - 1);




f2 = -alpha*ii(ind_c)';
z = Q_c\f2;
Rstar = z(ic_c)

rr = fzero(@(r) rfun(r,Q_c,f2,ic_c,ind_c) ,[rmin,rmax])


% with household distribution ---------------------------------------------
% Rstar = sum_k pi_k E(\int \a I(X_k(t))dt

Rstarall = zeros(n,1);
Rstarall(ind_c) = z;

Rstarhh = sum(pi'.*Rstarall)

h = [0.244, 0.334, 0.162, 0.159, 0.067, 0.034];
k = 1:6;
pi = h.*k./sum(k.*h); % household size distribution, probability that a 
% randomly selected person from the population belongs to a household
pi_c = pi(ind_c)/ sum(pi(ind_c)); % renormalise for just the states in c


rrhh = fzero(@(r) rfunhet(r,Q_c,f2,ind_c, pi_c) ,[rmin,rmax])


Td = log(2)/rrhh % doubling time