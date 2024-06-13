function [total, X_out, t_out, gross_I, gross_R, hh_dist] = SIRsimhet(N,beta,gamma,epsilon, max)
%SIRsim Simulates a household SIR model with three events, intenal infection,
%recovery and external infection. To allow for infinite possible households
%once the third event occurs a new household is added to the model with a 
%specified distribution from 2016 https://profile.id.com.au/australia/household-size.
%   TOTAL provides the number of susceptible, infected at any given event
%   time (total(1) = susceptible, total(2) = infected)
%   X_out records the number of susceptible and infected in each household.
%   with each household assigned 2 rows in the matrix
%   t_out is a vector that holds the corresponding event times
%   gross_I is the gross number of people that have been effected at each
%   event time
%   gross_R is the gross number of recovered people at each event time



no_events = 3;

% initial conditions
X = [N-1; 1];  % one infected case, total cases, hh1
t = 0;
gross_I = 1; %gross number of infected
gross_R = 0; % gross number of recovered
hh_dist = [];
    
total = X;
X_out = X;
t_out = 0;

% household distribution h(k) is household of size k

h = [0.244, 0.334, 0.162, 0.159, 0.067, 0.034];
k = 1:6;
pi = h.*k./sum(k.*h); % household size distribution, probability that a 
% randomly selected person from the population belongs to a household

while total(2,end) > 0 % while there exists some infected individuals
    
    household = 1; % initialise each iteration
    
    % step 1. Calculate the rates of each event given the current state.
    
    for hh = 1:length(X)/2 % transition rates for each house
        a(1+3*(hh-1)) = beta*X(2*hh)*X(2*hh-1)/(N-1); % infection- logistic growth beta*I*(N-I)/(N-1)
        a(2+3*(hh-1)) = gamma*X(2*hh);           % recovery
        a(3+3*(hh-1)) = epsilon*X(2*hh);         % external infection - X(2) or X(4)??
    end
    
    a0 = sum(a); % rate until the first event occurs - min(A1,A2,..)

    
    % step 2. Calculate the time to the next event.
    
    t = t - log(rand)/a0; 
    
   
    % step 3. Update the state.
    
    r = rand*a0; 
    
    % defining the intervals of events
    
    for i = 1:length(a) 
        if i == 1
            A(1,1) = 0;
            A(1,2) = a(1);
        else
            A(i,1) = sum(a(1:i-1));
            A(i,2) = sum(a(1:i));
        end
        
    end


    for ii = 1:length(A) % number of intervals
      if r > A(ii,1) && r < A(ii,2) == 1
            e = ii;
      end
    end

    event = e;
    while event > 3
        event = event - 3;
        household = household + 1; % which household the event takes place at
    end 
    
    S_row = household*2-1;
    I_row = household*2;
    
    
    if event == 1
        % infection of household
        X(S_row) = X(S_row) - 1;
        X(I_row) = X(I_row) + 1;
        gross_I = [gross_I gross_I(end)+1];
        gross_R = [gross_R gross_R(end)];
    elseif event == 2
        % recovery of household
        X(I_row) = X(I_row) - 1;
        gross_I = [gross_I gross_I(end)];
        gross_R = [gross_R gross_R(end)+1];
    elseif event == 3
        % external infection from all households - dont track which
        % individual spread it
        length_X = length(X);
        pi_pos = find(rand<cumsum(pi));
        hh_size = pi_pos(1);
        X(length_X+1) = hh_size-1;
        X(length_X+2) = 1;
        hh_dist = [hh_dist hh_size];
        gross_I = [gross_I gross_I(end)+1];
        gross_R = [gross_R gross_R(end)];
    end
    

        
    [length_X_out, ~] = size(X_out);
    if length(X) ~= length_X_out % incase of different dimensions
        X_out(length(X),:) = 0;
    end
        
    % record the time and state after each jump
    X_out = [X_out, X];
    t_out = [t_out, t];
    
    count = [sum(X(1:2:length(X))); sum(X(2:2:length(X)))];
    
    total = [total count];
        

    if t > max
        break
    end
end


end

