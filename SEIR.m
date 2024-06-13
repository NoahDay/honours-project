function [total, X_out, t_out, gross_I, gross_R] = SEIR(N,beta,gamma,sigma, alpha, max)
%SIRsim Simulates a household SIR model with three events, intenal infection,
%recovery and external infection. To allow for infinite possible households
%once the third event occurs a new household is added to the model (of size 3).
%   TOTAL provides the number of susceptible, infected at any given event
%   time (total(1) = susceptible, total(2) = infected)
%   X_out records the number of susceptible and infected in each household.
%   with each household assigned 2 rows in the matrix
%   t_out is a vector that holds the corresponding event times
%   gross_I is the gross number of people that have been effected at each
%   event time
%   gross_R is the gross number of recovered people at each event time



no_events = 4;

% initial conditions
X = [N-1; 0; 1];  % one infected case, total cases, hh1
t = 0;
gross_E = 0;
gross_I = 1; %gross number of infected
gross_R = 0; % gross number of recovered
hh_dist = N;
    
total = X;
X_out = X;
t_out = 0;

Ind = 1;

% household distribution h(k) is household of size k

h = [0.244, 0.334, 0.162, 0.159, 0.067, 0.034];
k = 1:6;
pi = h.*k./sum(k.*h); % household size distribution, probability that a 
% randomly selected person from the population belongs to a household

while t < 1%total(3,end) > 0 % while there exists some infected individuals
    
    % step 1. Calculate the rates of each event given the current state.
    
    for hh = 1:length(X)/3 % transition rates for each house
        if hh_dist(hh) == 1
            a(1+no_events*(hh-1)) = 0;
        else
            a(1+no_events*(hh-1)) = beta*X(3*hh)*X(3*hh-2)/(hh_dist(hh)-1); % infection- logistic growth beta*I*(N-I)/(N-1)
        end
        a(2+no_events*(hh-1)) = sigma*X(3*hh-1);        % latent progression
        a(3+no_events*(hh-1)) = gamma*X(3*hh);           % recovery
        a(4+no_events*(hh-1)) = alpha*X(3*hh);         % external infection - alphaI or epsilon S??
    end
    
    a0 = sum(a); % rate until the first event occurs - min(A1,A2,..)

    
    % step 2. Calculate the time to the next event.
    
    t = t - log(rand)/a0; 
    
   
    % step 3. Update the state.
    
    r = rand*a0; 
    
    % finding which event occurred

    pos = find(r<cumsum(a));
    e = pos(1);
    
    
    household = 1; % initialise each iteration
    event = e;
    while event > no_events
        event = event - no_events;
        household = household + 1; % which household the event takes place at
    end 
    
    S_row = household*3-2;
    E_row = household*3-1;
    I_row = household*3;
    
    
    if event == 1
        % infection of household
        X(S_row) = X(S_row) - 1;
        X(E_row) = X(E_row) + 1;
        gross_E = [gross_E gross_E(end)+1];
        gross_I = [gross_I gross_I(end)];
        gross_R = [gross_R gross_R(end)];
    elseif event == 2
        % latent progression
        X(E_row) = X(E_row) + 1;
        X(I_row) = X(I_row) + 1;
        gross_E = [gross_E gross_E(end)+1];
        gross_I = [gross_I gross_I(end)];
        gross_R = [gross_R gross_R(end)];
    elseif event == 3
        % recovery of household
        X(I_row) = X(I_row) - 1;
        gross_E = [gross_E gross_E(end)+1];
        gross_I = [gross_I gross_I(end)];
        gross_R = [gross_R gross_R(end)+1];
    elseif event == 4
        % external infection from all households - dont track which
        % individual spread it -ISNT APPLYING TO ALREADY INFECTED HH BUT
        % THIS MIGHT BE NEGLIBLE AS THE % OF THESE IS TINY
        length_X = length(X);
        pi_pos = find(rand<cumsum(pi));
        hh_size = pi_pos(1);
        X(length_X+1) = hh_size-1;
        X(length_X+2) = 1;
        X(length_X+3) = 0;
        hh_dist = [hh_dist, hh_size];
        gross_E = [gross_E, gross_E(end)];
        gross_I = [gross_I, gross_I(end)+1];
        gross_R = [gross_R, gross_R(end)];
    end
    

        
    [length_X_out, ~] = size(X_out);
    if length(X) ~= length_X_out % incase of different dimensions
        X_out(length(X),:) = 0;
    end
        
    % record the time and state after each jump
    X_out = [X_out, X];
    t_out = [t_out, t];
    
    count = [sum(X(1:3:length(X))); sum(X(2:3:length(X)));sum(X(3:3:length(X)))];
    
    total = [total count];
        

    if t > max
        break
    end
end
end
