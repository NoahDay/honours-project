

beta = 6; % internal infection
N = 2; % size of household
gamma = 1; % recovery
alpha = 1; % external infection
sigma = 20;
% X(1): Susceptible- S
% X(2): Latent - E
% X(3): Infected - I


% parameters
% beta - transmission within household 
% sigma - latent to infectious
% gamma - recovery of each infected individual

% mu = beta*S*I/(N-1); % infection- logistic growth beta*I*(N-I)/(N-1)
% lambda = gamma*I;       % recovery

no_stages = 3;

% Assigning coordinates to states
temp = [];

xxx = 0:N-1;
yyy = 0:N-1;
zzz = 0:N-1;
[XXX, YYY, ZZZ] = meshgrid(xxx,yyy, zzz);
[len,~] = size(XXX);
for i = 1:l^no_stages
    temp(:,i) = [XXX(i);YYY(i); ZZZ(i)];
end

coords = [N*eye(no_stages), temp];
pos = find(sum(coords) <= N);
x = coords(:,pos); % [S, E, I]


ss = x(1,:);
ee = x(2,:);
ii = x(3,:);

n = sum(length(x)); 

q = sparse(n,n);
% Entering recovery and infection rates into Q

% infection
for i=1:n
    if ss(i) > 0 && ii(i) > 0 % checking that there is someone to be infected and someone to spread the disease
        search = (ss == ss(i) - 1 & ee == ee(i) + 1 & ii ==  ii(i));
        b_pos = find(search);
        q(i, b_pos) = beta*ss(i)*ii(i)/(N-1);% + ss(i)*epsilon
    end
end

% latent progression

for i = 1:n
    if ee(i) > 0
        search = (ss == ss(i) & ee == ee(i) - 1 & ii ==  ii(i) + 1);
        sig_pos = find(search);
        q(i,sig_pos) = 2*sigma*ee(i);
    end
end

% recovery
for i = 1:n
    if ii(i) > 0
        search = (ss == ss(i) & ee == ee(i) & ii ==  ii(i) - 1);
        g_pos = find(search);
        q(i,g_pos) = gamma*ii(i);
    end
end


% 
% % external infection 
% for i=1:n
%     if ss(i) > 0 && ii(i) >= 0 % checking that there is someone to be infected and someone to spread the disease
%         search = x == [ss(i)-1; ii(i)+1];
%         pos = sum(search)==2;
%         q(i, find(pos)) = q(i, find(pos)) + ss(i)*epsilon;
%     end
% end


for i = 1:n
    q(i,i) = -sum(q(i,:));
end






