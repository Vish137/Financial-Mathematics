clear vars
%Option parameters%
sigma=0.2;
r=0.07;
K=100;
T=1;
q=2*r/sigma^2;

%Numerical parameters%
M=1000;
N=200;
dt=sigma^2*T/2/M;
dx=6*sigma*sqrt(T)/N;
lambda=dt/dx^2;

%Creating the x-axis with limits
x=zeros(1,N+1);
for n=1:N+1    
    x(n)=-3*sigma*sqrt(T)+(n-1)*dx;
end

%Solving for u(0,x)
u=zeros(1,N+1);
for n=1:N+1    
    u(n)=max(exp((q+1)/2*x(n))-exp((q-1)/2*x(n)),0);
end

%Solving for u(t + dt)
for m=2:M+1    
    w=zeros(1,N+1); 
    %Boundary conditions
    w(1)=0;    
    w(N+1)=exp((q+1)*x(N+1)/2+(q+1)^2*(m-1)*dt/4)-exp((q-1)*x(N+1)/2+(q-1)^2*(m-1)*dt/4);    
    alpha=(1+2*lambda)*ones(1,N-1);    
    b=zeros(1,N-1);    
    b(1)=lambda*w(1);
    b(N-1)=lambda*w(N+1);    
    b=b+u(2:N); 
    %Create diagonal matrix
    for n=2:N-1        
        alpha(n)=alpha(n)-lambda^2/alpha(n-1);% alpha hat        
        b(n)=b(n)+lambda*b(n-1)/alpha(n-1);% b hat    
    end
    w(N)=b(N-1)/alpha(N-1);    
    for n=N-2:-1:1        
        w(n+1)=(b(n)+lambda*w(n+2))/alpha(n);    
    end
    u=w;
end

S=zeros(1,N+1);
V=zeros(1,N+1);
for n=1:N+1    
    S(n)=K*exp(x(n));    
    V(n)=u(n)*K*exp(-(q-1)/2*x(n)-(q+1)^2/4*sigma^2*T/2) - (S(n) - K*exp(-r*T));
end

figure(1)
plot(S,V,'-r')
grid minor
title('European Black Scholes PDE - Implicit Scheme')
xlabel('Initial Price, S_0')
ylabel('Option Price, V(t,S)')
saveas(gcf,'Eu_Put_Implicit','png')