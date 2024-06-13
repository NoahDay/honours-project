clear
N = 3;
beta = 6;
mu = 6;
alpha = 1;
gamma = 1;
eps = 0.5;
startI = 0.5;
mmax=10;

 for m=2:mmax
    if mu == Inf % then it is just an sis model as people move from recovered
        % to susceptible instantly
        [Q, x] = QSIS(N, beta, gamma,0);
    else
        [Q, x] = QSIRS(N, beta, gamma, mu, eps);
    end

    ss = x(1,:); % number susceptible for each state
    ii = x(2,:); % number infected for each state

    ic = find(ss == N-1 & ii == 1);
    q = full(Q);

    [V, D] = eig(q','nobalance'); % D has the eigenvalues on the diagonal and V has
    % the corresponding column vectors for each eigen value
    [val,ind] = min(abs(diag(D)));
    V_vec = V(:,ind);
    pi = V(:,ind)./sum(V(:,ind)); % this is the stationary distribution corresponding to an 
    % eigenvalue of 0 which is required Q*pi = 0, and the inital conditions


    mI=sum(V_vec.*ii')/(N*sum(V_vec)) % I(X_eps(inf))/N

    if mI < (eps/alpha)
        eps = eps - 0.5^ mmax;
    elseif mI >= (eps/alpha)
        eps = eps + 0.5^m;
    else
        display('Failure in mean')
    end
 end
eps