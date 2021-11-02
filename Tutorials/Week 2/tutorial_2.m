%|--------------------------|%
%| FMAT3888 Tutorial Week 2 |%
%| Author: Vishaal Lingam   |%
%| Date: 17-08-2021         |%
%|--------------------------|%

% Q1.
% Evaluating the intergral of sin(x) from zero to pi
n=7; % maximum exponenet and number of trials
J = linspace(1,n,n); % vector of exponents
MU = zeros(1,length(J)); % zeros for final evaluation of average

% Looping over 'n' seperate trials
for i=1:length(J)
    N=10^J(i); % number of simulations
    Y=rand(1,N)*pi; % generate N iid copies of U(0,pi)
    mu=0;
    for k=1:N    
        mu=mu+sin(Y(k)); % sum of all random results
    end
    MU(i)=mu/N*pi; % evaluation of average over all simulations
end 

% Creating a plot to visualise the convergece towards the real solution as n>>10
figure(1);
plot(J,MU,'b','LineWidth',2);
hold on;
yline(2,'-r','LineWidth',2);
hold off;
ylabel('Final Result')
title([' Evaluating $ \int_{0}^{\pi} \sin(x)dx $ by Monte-Carlo'],'interpreter','latex')
xlabel('Number of Simulations (10^x)')
legend({'Computed Value','Accepted Value'})
saveas(gcf,'Q1MonteCarlo','png')

% Q2 a)
% Evaluating the Black-Scholes model of an European option 
S0=100; % initial underlying stock price
K=100; % strike price of option
T=1; % time scale
r=0.07; % interest rate
sigma=0.25; % volatility 
N=10^6; % Number of simulations
Y=randn(1,N);
Z=zeros(1,N); % will be used for price for S_T
W=zeros(1,N);

% Looping over given number of simulations
for k=1:N
    Z(k)=S0*exp((r-0.5*sigma^2)*T + (sigma*sqrt(T)*Y(k)));
    W(k)=max(Z(k)-K,0);
end 

price=exp(-r*T)*sum(W)/N;
figure(2);
histogram(W,100)
display(price)
saveas(gcf,'Q2aBlackScholes','png')

% Q2 b)
% Evaluating the Black-Scholes model of European option and plot the convergence diagram for succesive Monte Carlo simulaitons
S0=100; % initial underlying stock price
K=100; % strike price of option
T=1; % time scale
r=0.07; % interest rate
sigma=0.25; % volatility 
N=10^6; % Number of simulations
Y=randn(1,N);
Z=zeros(1,N); % will be used for price for S_T
W=zeros(1,N);

for k=1:N
    Z(k)=S0*exp((r-0.5*sigma^2)*T + (sigma*sqrt(T)*Y(k)));
    W(k)=max(Z(k)-K,0);
end

U=zeros(1,N); % for \bar f 
for k=2:N
    U(k)=U(k-1)+W(k);
end 

for k=1:N
    U(k)=exp(-r*T)*U(k)/k;
end

figure(3);
plot(1:N,U)
xlabel('Number of simulations');
ylabel('Price of call option');
title('Convergence diagram');
saveas(gcf,'Q2bBlackScholes','png')

% Q2 c)
% Redoing parts a) and c) with antithetic variables
S=100;
K=100;
r=0.07;
sig=0.25;
T=1;
N=10^3;
Y=randn(1,N);
Z=zeros(1,N);%used for f(X) roughly speaking
V=zeros(1,N);%used for f(-X)
W=zeros(1,N);%(f(x)+f(-x))/2
for k=1:N    
    Z(k)=S*exp((r-0.5*sig^2)*T+sig*sqrt(T)*Y(k));    
    V(k)=S*exp((r-0.5*sig^2)*T-sig*sqrt(T)*Y(k));    
    W(k)=0.5*(max(Z(k)-K,0)+max(V(k)-K,0));
end
price=exp(-r*T)*sum(W)/N;
figure(4);
histogram(W,100)
saveas(gcf,'Q2cBlackScholes','png')

% Q2 d)
% Redoing part a) and b) with S0=50
% Evaluating the Black-Scholes model of an European option 
S0=50; % initial underlying stock price
K=100; % strike price of option
T=1; % time scale
r=0.07; % interest rate
sigma=0.25; % volatility 
N=10^6; % Number of simulations
Y=randn(1,N);
Z=zeros(1,N); % will be used for price for S_T
W=zeros(1,N);
% Looping over given number of simulations
for k=1:N
    Z(k)=S0*exp((r-0.5*sigma^2)*T + (sigma*sqrt(T)*Y(k)));
    W(k)=max(Z(k)-K,0);
end 
price=exp(-r*T)*sum(W)/N;
figure(5);
histogram(W,100)
display(price)
saveas(gcf,'Q2d1BlackScholes','png')

% Evaluating the Black-Scholes model of European option and plot the convergence diagram for succesive Monte Carlo simulaitons
S0=50; % initial underlying stock price
K=100; % strike price of option
T=1; % time scale
r=0.07; % interest rate
sigma=0.25; % volatility 
N=10^6; % Number of simulations
Y=randn(1,N);
Z=zeros(1,N); % will be used for price for S_T
W=zeros(1,N);
for k=1:N
    Z(k)=S0*exp((r-0.5*sigma^2)*T + (sigma*sqrt(T)*Y(k)));
    W(k)=max(Z(k)-K,0);
end
U=zeros(1,N); % for \bar f 
for k=2:N
    U(k)=U(k-1)+W(k);
end 
for k=1:N
    U(k)=exp(-r*T)*U(k)/k;
end
figure(6);
plot(1:N,U)
xlabel('Number of simulations');
ylabel('Price of call option');
title('Convergence diagram');
saveas(gcf,'Q2d2BlackScholes','png')

%Q2 e)
% Redoing part d) with importance sampling
S=50;
K=100;
r=0.07;
sig=0.25;
T=1;
beta=-(log(K/S)-(r-0.5*sig^2)*T)\sig\sqrt(T);
N=10^4;
Y=randn(1,N);
W=zeros(1,N);%for g_beta(x)
for k=1:N    
    W(k)=max(S*exp((r-0.5*sig^2)*T+sig*sqrt(T)*(Y(k)-beta))-K,0)*exp(beta*Y(k)-0.5*beta^2);
end
price=exp(-r*T)*sum(W)/N;
figure(7);
histogram(W,50)

S=50;
K=100;
r=0.07;
sig=0.25;
T=1;
beta=-(log(K/S)-(r-0.5*sig^2)*T)\sig\sqrt(T);
N=10^4;
Y=randn(1,N);
W=zeros(1,N);%for g_beta(x)
for k=1:N    
    W(k)=max(S*exp((r-0.5*sig^2)*T+sig*sqrt(T)*(Y(k)-beta))-K,0)*exp(beta*Y(k)-0.5*beta^2);
end
U=zeros(1,N);
U(1)=W(1);
for k=2:N    
    U(k)=U(k-1)+W(k);
end
for k=1:N    
    U(k)=exp(-r*T)*U(k)/k;
end
figure(8);
plot(1:N,U)
xlabel('Number of simulations')
ylabel('Price of call option')
title('Convergence diagram with importance sampling')
hold on
