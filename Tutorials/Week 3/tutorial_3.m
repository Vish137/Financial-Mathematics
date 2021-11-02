%|--------------------------|%
%| FMAT3888 Tutorial Week 3 |%
%| Author: Vishaal Lingam   |%
%| Date: 24-08-2021         |%
%|--------------------------|%

%Q1 a).
% Euler discretization: S_{t+dt} = S_{t} + rS_{t}dt + σS_{t}sqrt{dt}B_{t}
% Milstein discretization: S_{t + dt} = S_{t} + rS_{t}dt + σS_{t}sqrt{dt}B_{t} + 0.5(σ^2)*dt*(Z^2 - 1)

%Q1 b).
T = 1;
r = 0.07;
sigma = 0.25;
S0 = 100;
K = [75,85,95];
N = 10000;
M = 1000;
dt = T/M;
mu = 0;

for i=1:N
    S = zeros(1,M+1);
    xi = randn(1,M);
    S(1) = S0;
    % Appending S with discretized Black-Scholes eqn
    for k =1:M
        S(k+1) = S(k) + r*S(k)*dt + sigma*S(k)*sqrt(dt)*xi(k);
    end

    % Applying the condition on expectation value
    if min(S) > K(2)
        mu = mu + 1;
    else
        mu = mu + 0;
    end

end
mu = exp(-1*r*T)*mu/N;
disp(mu)
% Prints ~0.74 for K(1), ~0.50 for K(2), ~0.18 for K(3) as expected


%Q1 c).
% Milstein discretization: S_{t + dt} = S_{t} + rS_{t}dt + σS_{t}sqrt{dt}Z + 0.5(σ^2)*dt*(Z^2 - 1)
T = 1;
r = 0.07;
sigma = 0.25;
S0 = 100;
K = [75,85,95];
N = 10000;
M = 1000;
dt = T/M;
mu = 0;

for i=1:N
    S = zeros(1,M+1);
    xi = randn(1,M);
    S(1) = S0;
    % Appending S with Milstein discretized Black-Scholes eqn
    for k =1:M
        S(k+1) = S(k) + r*S(k)*dt + sigma*S(k)*sqrt(dt)*xi(k) + 0.5*(sigma^2)*dt*(xi(k)^2 - 1);
    end

    %Applying the condiotion on expectation value
    if min(S) > K(2)
        mu = mu + 1;
    else
        mu = mu + 0;
    end

end
mu = exp(-1*r*T)*mu/N;
disp(mu)
%Prints ~0.75 for K(1), ~0.51 for K(2), ~0.19 for K(3)


%Q2 a).
% Euler discretization:
% S_{t + dt} = S_{t} + rS_{t}dt + ρ sqrt{Y_{t}dt}S_{t}B_{t} + sqrt{1 - ρ^2} sqrt{Y_{t}dt}S_{t}Z_{t}
% Y_{t + dt} = Y_{t} + (θ + κY_{t})dt + β sqrt{Y_{t}dt}B_{t}

S0 = 100;
Y0 = 0.2;
theta = 0.3;
kap = 2;
beta = 0.25;
r = 0.05;
rho = 0.5;
T = 1;
M=10000;
dt = T/M;

%Plotting
S = zeros(1,M+1);
Y = zeros(1,M+1);
xi=randn(1,M);
yi=randn(1,M);
S(1)=S0;
Y(1)=Y0;
for k=1:M
    S(k+1) = S(k) + r*S(k)*dt + rho*sqrt(Y(k)*dt)*S(k)*yi(k)+ sqrt(1 - rho^2)*sqrt(Y(k)*dt)*S(k)*xi(k);
    Y(k+1) = Y(k) + (theta + kap*Y(k))*dt + beta*sqrt(Y(k)*dt)*yi(k);
end
figure(1)
plot(0:dt:T,S)
title('S(t)')
xlabel('Time')
ylabel('S(t)')

figure(2)
plot(0:dt:T,Y)
title('Y(t)')
xlabel('Time')
ylabel('Y(t)')

%Q2 b). & c).
S0=100;
Y0=0.2;
theta=0.3;
kappa=2;
beta=0.25;
r=0.05;
rho=0.5;
rrho=sqrt(1-rho^2);
T=1;
N=5000;
M=5000;
dt=T/M;
K=100;
a=0;
b=0;
for n=1:N 
    xi=randn(1,M); 
    eta=randn(1,M); 
    S=zeros(1,M+1); 
    Y=zeros(1,M+1); 
    S(1)=S0; 
    Y(1)=Y0; 
    for m=1:M     
        S(m+1)=S(m)+r*S(m)*dt + sqrt(Y(m))*S(m)*(rho*sqrt(dt)*xi(m)+rrho*sqrt(dt)*eta(m));     
        Y(m+1)=Y(m)+(theta-kappa*Y(m))*dt+beta*sqrt(Y(m))*sqrt(dt)*xi(m); 
    end
    a=a+max(S(M+1)-K,0);%for european call 
    b=b+max(max(S)-K,0);%for lookback
end
call=exp(-r*T)*a/N;
lookback=exp(-r*T)*b/N;




