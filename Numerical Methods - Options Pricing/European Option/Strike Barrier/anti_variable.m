S=100;
K=100;
r=0.07;
sig=0.25;
T=1;
N=5000;
M=5000;
Y=randn(1,N);
Z=zeros(1,N);%used for f(X) roughly speaking
V=zeros(1,N);%used for f(-X)
W=zeros(1,N);%(f(x)+f(-x))/2
M=zeros(1,N);
for k=1:N    
    Z(k)=S*exp((r-0.5*sig^2)*T+sig*sqrt(T)*Y(k));    
    V(k)=S*exp((r-0.5*sig^2)*T-sig*sqrt(T)*Y(k));    
    W(k)=0.5*(max(Z(k)-K,0)+max(V(k)-K,0));
end
for i=1:M       
    for k=1:N        
        Z(k)=S*exp((r-0.5*sig^2)*T+sig*sqrt(T)*Y(k));
        V(k)=S*exp((r-0.5*sig^2)*T-sig*sqrt(T)*Y(k)); 
    end  
    if min(Z(i))<=L        
        a=0;    
    else
        a=1;    
    end  
    if min(V(i))<=L
        b=0;    
    else
        b=1; 
    end
    M(i)=b;
    W(i)=a;    
end
U=zeros(1,N);%for \bar f
U(1)=0.5*(W(1)+M(1));
for k=2:N    
    U(k)=U(k-1)+0.5*(W(k)+M(k));
end
for k=1:N
    U(k)=exp(-r*T)*U(k)/k;
end 
price=exp(-r*T)*sum((W+M)*0.5)/1
plot(1:N,U)
xlabel('Number of simulations')
ylabel('Price of barrier option')
title('Convergence diagram antithetic variable')
hold on;