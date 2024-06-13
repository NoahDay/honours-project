% abs.gov.au Aus ave hhsize = 2.6, using parameters based on basic
% parameter Ross 
clear

beta = 6; % internal infection
N = 1; % size of household
houses = 10;
maxhouses=50;
gamma = 1; % recovery
mu = 0.5; % waning immunity
alpha = 1; % external infection
tmax = 20;
[Q, x] = Qgen(N, beta, gamma);

n = sum(1:N+1);
I_0 = find(x(2,:)==0); % finding the absorbing states

Q_c = Q; % transient states only

for i = N+1:-1:1
    Q_c(I_0(i),:) = [];
    Q_c(:,I_0(i)) = [];
end

f = x(2,:)';
for i = length(I_0):-1:1
    f(I_0(i)) = [];
end

% solving E[Gamma|X(0) = i] -----------------------------------------------

 e_temp = -Q_c\f; 

 e = zeros(n,1); 

 pos = find(x(2,:)~=0);

 for i = 1:length(e_temp)
     e(pos(i)) = e_temp(i);
 end

 %e(8) is e_1
ss = x(1,:); % number susceptible for each state
ii = x(2,:); % number infected for each state
ind_c = find(ii==1 & ss == N-1) % intial conditions
R_star = e(ind_c)

initial_i = find(ii==1 & ss == N-1);
Cind = find(ii>0);
initial_c = find(ind_c==initial_i); % finding the position in C which corresponds to intial conditions in S


% figure(3)
% plot(1:M,R_star)
% title('$R_*$ for a range of household sizes','Interpreter','latex')
% xlabel('Size of household','Interpreter','latex')
% ylabel('$R_*$','Interpreter','latex')

% figure(4)
% plot(1:M,1./R_star)
% title('$\alpha_{critical}$ for a range of household sizes','Interpreter','latex')
% xlabel('Size of household','Interpreter','latex')
% ylabel('$\alpha_{critical}$','Interpreter','latex')



% finding r
% using proposition 1 we can solve this expectation of a path integral like
% we did for y


% Only search over possible positive values of early growth
rmin = 0;
rmax = gamma*(R_star - 1);


f2 = -alpha*ii(Cind);
z = Q_c\f2;
Rstar = z(initial_c);

%rr = fzero( @(r)get_val(r,Qc,f2,Iind,Cind), [rmin, rmax]);

rr = fzero( @(r) rfun(r,Q_c,f2,initial_c,Cind),[rmin,rmax])