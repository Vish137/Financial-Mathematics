S=100;
r=0.07;
sig=0.25;
T=1;
N=5000;
M=5000;
L=95;
Y=randn(1,N);
Z=zeros(1,N);%will be used for price for S_T
W=zeros(1,N);

for i=1:M       
    for k=1:N        
        Z(k)=S*exp((r-0.5*sig^2)*T+sig*sqrt(T)*Y(k));    
    end
    
    if min(Z(i))<=L        
        a=0;    
    else
        a=1;    
    end
    
    W(i)=a;    
end
U=zeros(1,N);%for \bar f
U(1)=W(1);
for k=2:N    
    U(k)=U(k-1)+W(k);
end
for k=1:N
    U(k)=exp(-r*T)*U(k)/k;
end 
price=exp(-r*T)*sum(W)/N
plot(1:N,U)
xlabel('Number of simulations')
ylabel('Price of barrier option')
title('Convergence diagram using Motne Carlo')
hold on;
