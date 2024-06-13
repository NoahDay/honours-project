clear 

S_mean = [];
I_mean = [];
total_I = [];
total_t = [];
maximum = 5;
delta = linspace(0,maximum,1000);
it = 3;


for i = 1:it
    
    N = 4;

    beta = 3;
    gamma = 0.5;
    alpha = 1;
    
    [total, X_out, t_out, gross_I, gross_R] = SIRsim(N,beta,gamma,alpha,maximum);
    
    for j = 1:length(delta)
        pos = find(delta(j) <= t_out);
        [~, n] = size(pos);
        if n == 0
           break;
        end
        E_X(i,j) = gross_I(pos(1));
        E_R(i,j) = gross_R(pos(1));
    end
    total_t(i) = t_out(end);
    
    if i ==1
        x1 = [gross_I;t_out;X_out];
    elseif i == 2
            x2 = [gross_I;t_out];
        elseif i == 3
            x3 = [gross_I;t_out;X_out];
    end
end

maxx = max(max(E_X));
[x,y]=find(E_X==maxx);


figure(1)
clf
hold on
plot(delta(1:length(E_X)),E_X(1:3,:),'--',delta(1:length(E_X)),mean(E_X),delta(1:length(E_X)),mean(E_R))
%plot(delta(1:length(E_X)),E_X)
ylabel('Number of infected')
xlabel('time')
str = sprintf('Sample paths of the simulation, %d iterations', it);
title(str)
legend('Sample 1','Sample 2','Sample 3','Mean infected', 'Mean recovered')
hold off

figure(2)
plot(delta(1:length(E_X)),E_X(mode(x),:))
str2 = sprintf('Maximum from %d simulations', it);
title(str2)

figure(3)
clf
hold on
plot(x1(2,:),x1(1,:),x2(2,:),x2(1,:),x3(2,:),x3(1,:))
str3 = sprintf('True stochastic paths');
title(str3)
ylabel('Number of infected')
xlabel('time')
legend('Sample 1','Sample 2','Sample 3')
hold off