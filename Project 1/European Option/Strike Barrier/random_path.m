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
U=zeros(1,N+1);%for \bar f
U(1)=W(1);
for k=2:N    
    U(k+1)=U(k)+W(k+1);
end

for k=1:N
    U(k+1)=exp(-r*T)*U(k+1)/(k+1);
end 
price=exp(-r*T)*sum(W)/N
U
sum(W)
plot(1:N+1,U)
xlabel('Number of simulations')
ylabel('Price of barrier option')
title('Convergence diagram use Eluer method')
hold on;
