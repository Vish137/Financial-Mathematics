%|--------------------------|%
%| FMAT3888 Tutorial Week 5 |%
%| Author: Vishaal Lingam   |%
%| Date: 14-09-2021         |%
%|--------------------------|%


%Q1.
%OPTION PARAMETERS
sigma=0.5;
r=0.1;
K=100;
T=1;

%NUMERICAL PARAMETERS
q=2*r/sigma^2;
M=1000;
N=1000;
dt=sigma^2*T/2/M;
dx=6*sigma*sqrt(T)/N;
lambda=dt/dx^2;

%Creating x-axis
x=zeros(1,N+1);
for n=1:N+1    
    x(n)=-3*sigma*sqrt(T)+(n-1)*dx;
end

%Solving u(t + dt)
u=zeros(1,N+1);
for n=1:N+1    
    u(n)=max(exp((q-1)/2*x(n))-exp((q+1)/2*x(n)),0);
end
%Implement Crank-Nicolson 
for m=2:M+1   
    w=zeros(1,N+1); 
    %Establish boundary conditions
    w(1)=exp((q-1)*x(1)/2+(q-1)^2*(m-1)*dt/4)-exp((q+1)*x(1)/2+(q+1)^2*(m-1)*dt/4);   
    w(N+1)=0;   
    %Create tridiagonal matrix
    alpha=(1+lambda)*ones(1,N-1);    
    b=zeros(1,N-1);    
    for n=1:N-1        
        b(n)=lambda/2*u(n)+(1-lambda)*u(n+1)+lambda/2*u(n+2);    
    end
    b(1)=b(1)+lambda/2*w(1);
    b(N-1)=b(N-1)+lambda/2*w(N+1);    
    for n=2:N-1        
        alpha(n)=alpha(n)-lambda^2/4/alpha(n-1);        
        b(n)=b(n)+lambda/2*b(n-1)/alpha(n-1);    
    end
    w(N)=b(N-1)/alpha(N-1);
    for n=N-2:-1:1        
        w(n+1)=(b(n)+lambda/2*w(n+2))/alpha(n);    
    end
    u=w;
end
%Calculate payoff
S=zeros(1,N+1);
V=zeros(1,N+1);
for n=1:N+1    
    S(n)=K*exp(x(n));    
    V(n)=u(n)*K*exp(-(q-1)/2*x(n)-(q+1)^2/4*sigma^2*T/2);
end
figure(1)
plot(S,V,'blue');
title('European Option - Crank-Nicolson Metod')
xlabel('Stock Price, $S$', 'Interpreter','latex')
ylabel('Option Price, $V(0,S)$','Interpreter','latex')
grid minor
saveas(gcf,'Q1a_Eu_Put_CN','png')

clear vars
%Q2.
%OPTION PARAMETERS
sigma=0.5;
r=0.1;
K=100;
T=1;

%NUMERICAL PARAMETERS
q=2*r/sigma^2;
M=1000;
N=1000;
dt=sigma^2*T/2/M;
dx=6*sigma*sqrt(T)/N;
lambda=dt/dx^2;
c=zeros(1,M);

x=zeros(1,N+1);
for n=1:N+1    
    x(n)=-3*sigma*sqrt(T)+(n-1)*dx;
end
u=zeros(1,N+1);
for n=1:N+1    
    u(n)=max(exp((q-1)/2*x(n))-exp((q+1)/2*x(n)),0);
end
for m=2:M+1    
    w=zeros(1,N+1);    
    w(1)=exp((q-1)*x(1)/2+(q-1)^2*(m-1)*dt/4)-exp((q+1)*x(1)/2+(q+1)^2*(m-1)*dt/4);    
    w(N+1)=0;    
    alpha=(1+lambda)*ones(1,N-1);    
    b=zeros(1,N-1);    
    for n=1:N-1        
        b(n)=lambda/2*u(n)+(1-lambda)*u(n+1)+lambda/2*u(n+2);    
    end
    b(1)=b(1)+lambda/2*w(1);
    b(N-1)=b(N-1)+lambda/2*w(N+1);    
    for n=2:N-1        
        alpha(n)=alpha(n)-lambda^2/4/alpha(n-1);        
        b(n)=b(n)+lambda/2*b(n-1)/alpha(n-1);    
    end
    w(N)=b(N-1)/alpha(N-1);    
    for n=N-2:-1:1        
        w(n+1)=(b(n)+lambda/2*w(n+2))/alpha(n);    
    end
    
    %Check if striking out is better
    g=zeros(1,N+1);   
    for n=1:N+1        
        g(n)=exp((q+1)^2/4*(m-1)*dt)*max(exp((q-1)/2*x(n))-exp((q+1)/2*x(n)),0);    
    end
    u=max(w,g); 
    
    %Finding the exercise boundary
    n=1;   
    while u(n)==g(n)        
        n=n+1;    
    end
    c(m-1)=x(n-1);
end
S=zeros(1,N+1);
V=zeros(1,N+1);
for n=1:N+1    
    S(n)=K*exp(x(n));    
    V(n)=u(n)*K*exp(-(q-1)/2*x(n)-(q+1)^2/4*sigma^2*T/2);
end
figure(2)
plot(S,V,'red');
title('American Put - Crank-Nicolson')
xlabel('Stock Price, $S$', 'Interpreter','latex')
ylabel('Option Price, $V(0,S)$','Interpreter','latex')
grid minor
saveas(gcf,'Q1a_Am_Put_CN','png')

figure(3)
dtt=2*dt/sigma^2;
plot(T-dtt:-dtt:0,K*exp(c),'Color',[0 0.5 0]);
title('American Put Exercise Boundary')
xlabel('Contract Period, $T$', 'Interpreter','latex')
ylabel('Stock Price, $S$','Interpreter','latex')
grid minor
saveas(gcf,'Q1a_Am_Put_EB','png')