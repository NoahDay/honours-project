clear all
% approximating g(m) with single sized households
S_mean = [];
I_mean = [];

total_I = [];

    fprintf('progress')
for i = 1:1000
    
    beta = 0; % internal infection
    N = 3; % size of household
    houses = 10;
    maxhouses=30;
    gamma = 1; % recovery
    mu = 0.5; % waning immunity
    alpha = 1; % external infection
    tmax = 20;


    [total, X_out, t_out, gross,X] = SIRsim(N,beta,gamma,alpha,tmax);
    
    
    total_I = [total_I, gross];
    S_mean(i) = X_out(1,end);
    I_mean(i) = X_out(2,end);
    [m ~] = size(X_out);
    hh_final(i) = m/2;
    hh_infected(i) = max(gross);
    if mod(i,10) == 0
        str = i;
        disp(str);
    end

end

bins = [0:30];
figure(1)
idx = hh_infected ~= 0;
g = hh_infected-1;
[N,edges] = histcounts(g, bins,'Normalization','pdf');
plot(edges(2:end), N);
title('PDF of the number of houses infected from 1 infected individual, simulation','Interpreter','latex')
ylabel('Probability','Interpreter','latex')
xlabel('Number of households, $m$','Interpreter','latex')


% approximating R*, take sum ( g(m)*m)

R_star = sum(N.*edges(1:end-1))-1 % minus one as this is including the initally infected one
