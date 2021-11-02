S0=100;
r=0.07;
sig=0.25;
T=1;
N=5000;
M=5000;
dt=T/M;
L=85;
mu=0;
W=zeros(1,M+1);
for i=1:N
    S=zeros(1,M+1);
    xi=randn(1,M);
    S(1)=S0;
    for k=1:M
        S(k+1)=S(k)+S(k)*r*dt+S(k)*sig*sqrt(dt)*xi(k);
    end
    if min(S)<=L
        a=0;
    else
        a=1;
    end
    %mu=mu+a;
    W(1)=1;
    W(i+1)=a;
end
sum(W)
U1=zeros(1,N+1);%for \bar f
U1(1)=W(1);
for k=2:N    
    U1(k+1)=U1(k)+W(k+1);
end

for k=1:N
    U1(k+1)=exp(-r*T)*U1(k+1)/(k+1);
end 
price=exp(-r*T)*sum(W)/N
plot(1:N+1,U1)
xlabel('Number of simulations')
ylabel('Price of barrier option')
title('Convergence diagram use Eluer method')
legend('Eluer method')
hold on;

clear all

r=0.07;
sig=0.25;
S0=100;
T=1;L=85;
N=5000;
M=5000;
dt=T/M;
W=zeros(1,M+1);
for i=1:N    
    S=zeros(1,M+1);    
    xi=randn(1,M);    
    S(1)=S0;    
    for k=1:M        
        S(k+1)=S(k)+S(k)*r*dt+S(k)*sig*sqrt(dt)*xi(k)+sig^2*S(k)*dt/2*(xi(k)^2-1);    
    end
    if min(S)<=L        
        a=0;    
    else
        a=1;    
    end
    W(1)=1;
    W(i+1)=a;
end
sum(W)
U=zeros(1,N+1);%for \bar f
U(1)=W(1);
for k=2:N    
    U(k+1)=U(k)+W(k+1);
end

for k=1:N
    U(k+1)=exp(-r*T)*U(k+1)/(k+1);
end 
price=exp(-r*T)*sum(W)/N
plot(1:N+1,U)
xlabel('Number of simulations')
ylabel('Price of barrier option')
title('Convergence diagram use Milstein method')
legend('Milstein method')
hold on;