%OPTION PARAMETERS%
tic
S = 100; % option buy price
K = 100; % option strike price
r = 0.05; % interest rate
T=1; % period to expiry
sigma=0.25; % price volatility
NSteps=24; % number of simulations in [1,T]
%Simulating with Halton random numbers%
N=5000;
G = zeros(1,N);
for i=1:N
    q = qrandstream('halton',NSteps,'Skip',1e3,'Leap',1e2);
    RandMat = qrand(q,i);
    z_RandMat = norminv(RandMat,0,1);
    dt = T/NSteps;
    C = (r-0.5*sigma^2)*dt + sigma*sqrt(dt)*z_RandMat;
    Paths = cumsum([log(S)*ones(i,1),C],2);
    % Underlying asset price path
    SPaths = exp(Paths);
    Payoff = zeros(i,1);
    Payoff = max(0, mean(SPaths(:,2:(NSteps+1)),2)- K);
    G(i) = mean(exp(-r*T) * Payoff);
end 
figure(1)
plot(G,'-b')
grid minor
grid on
title('Convergence Diagram (Asian Option, Halton Random Numbers)')
ylabel('Option Price (Numerical)')
xlabel('No. of Simulations')
lgd=sprintf('$ %0.5f',G(end));
legend({lgd})
saveas(gcf,'Asian_MC_Halton','png')

%Simulating with Sobol random numbers%
N=5000;
G = zeros(1,N);
for i=1:N
    q = qrandstream('sobol',NSteps,'Skip',1e3,'Leap',1e2);
    RandMat = qrand(q,i);
    z_RandMat = norminv(RandMat,0,1);
    dt = T/NSteps;
    C = (r-0.5*sigma^2)*dt + sigma*sqrt(dt)*z_RandMat;
    Paths = cumsum([log(S)*ones(i,1),C],2);
    % Underlying asset price path
    SPaths = exp(Paths);
    Payoff = zeros(i,1);
    Payoff = max(0, mean(SPaths(:,2:(NSteps+1)),2)- K);
    G(i) = mean(exp(-r*T) * Payoff);
end 
figure(2)
plot(G,'-b')
grid minor
grid on
title('Convergence Diagram (Asian Option, Sobol Random Numbers)')
ylabel('Option Price (Numerical)')
xlabel('No. of Simulations')
lgd=sprintf('$ %0.5f',G(end));
legend({lgd})
saveas(gcf,'Asian_MC_Sobol','png')

%Simulating with quasi-random numbers%
N=5000;
G = zeros(1,N);
for i=1:N
    z_RandMat = norminv(rand(i,NSteps),0,1);
    dt = T/NSteps;
    C = (r-0.5*sigma^2)*dt + sigma*sqrt(dt)*z_RandMat;
    Paths = cumsum([log(S)*ones(i,1),C],2);
    % Underlying asset price path
    SPaths = exp(Paths);
    Payoff = zeros(i,1);
    Payoff = max(0, mean(SPaths(:,2:(NSteps+1)),2)- K);
    G(i) = mean(exp(-r*T) * Payoff);
end 
figure(3)
plot(G,'-b')
grid minor
grid on
title('Convergence Diagram (Asian Option, \sim IG(0,1))')
ylabel('Option Price (Numerical)')
xlabel('No. of Simulations')
lgd=sprintf('$ %0.5f',G(end));
legend({lgd})
saveas(gcf,'Asian_MC_Inv_Norm','png')
toc
