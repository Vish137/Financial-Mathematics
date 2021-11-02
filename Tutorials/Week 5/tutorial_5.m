%|--------------------------|%
%| FMAT3888 Tutorial Week 5 |%
%| Author: Vishaal Lingam   |%
%| Date: 7-09-2021         |%
%|--------------------------|%

%Q1 a)
%Using the implicit method to solve the European BS PDE
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
    V(n)=u(n)*K*exp(-(q-1)/2*x(n)-(q+1)^2/4*sigma^2*T/2);
end
figure(1)
plot(S,V,'-b')
grid minor
title('European Black Scholes PDE - Implicit Scheme')
xlabel('Initial Price, S_0')
ylabel('Option Price, V(t,S)')
saveas(gcf,'q1a','png')

clear all
%Q1 b)
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
x=zeros(1,N+1);
for n=1:N+1    
    x(n)=-3*sigma*sqrt(T)+(n-1)*dx;
end
u=zeros(1,N+1);
for n=1:N+1    
    u(n)=max(exp((q+1)/2*x(n))-exp((q-1)/2*x(n)),0);
end
for m=2:M+1    
    w=zeros(1,N+1);    
    w(1)=0;    
    w(N+1)=exp((q+1)*x(N+1)/2+(q+1)^2*(m-1)*dt/4)-exp((q-1)*x(N+1)/2+(q-1)^2*(m-1)*dt/4);    
    alpha=(1+lambda)*ones(1,N-1);    
    b=zeros(1,N-1);% for B*u_m+b_{m+1}+b_m    
    for n=1:N-1        
        b(n)=lambda/2*u(n)+(1-lambda)*u(n+1)+lambda/2*u(n+2);    
    end
    b(1)=b(1)+lambda/2*w(1);
    b(N-1)=b(N-1)+lambda/2*w(N+1);    
    for n=2:N-1 %solve Mx=b for tridiagonal matrix M. alpha=1+lambda, beta=gamma=-lambda/2        
        alpha(n)=alpha(n)-lambda^2/4/alpha(n-1);        
        b(n)=b(n)+lambda/2*b(n-1)/alpha(n-1);    
    end
    w(N)=b(N-1)/alpha(N-1);
    for n=N-2:-1:1        
        w(n+1)=(b(n)+lambda/2*w(n+2))/alpha(n);    
    end
    u=w;
end
S=zeros(1,N+1);V=zeros(1,N+1);
for n=1:N+1    
    S(n)=K*exp(x(n));    
    V(n)=u(n)*K*exp(-(q-1)/2*x(n)-(q+1)^2/4*sigma^2*T/2);
end
figure(2)
plot(S,V,'-b')
grid minor
title('European Black Scholes PDE - NC Scheme')
xlabel('Initial Price, S_0')
ylabel('Option Price, V(t,S)')
saveas(gcf,'q1b','png')

clear all
%Q2 a)
dx=2^(-6);
lambda=[2/5,2/6,2/7];
dt=lambda(1)*dx^2;
M=1/dt;N=2/dx;
u=zeros(M+1,N+1);
x=zeros(1,N+1);
for n=1:N+1    
    x(n)=-1+(n-1)*dx;    
    u(1,n)=x(n)^2;
end
u(:,1)=ones(1,M+1);
u(:,N+1)=ones(1,M+1);
for m=1:M    
    p1=dt/dx^2*((m-1)^2*dt^2+0.5)+dt/2/dx;    
    p2=dt/dx^2*((m-1)^2*dt^2+0.5)-dt/2/dx;    
    for n=2:N        
        u(m+1,n)=p1*u(m,n+1)+p2*u(m,n-1)+(1-p1-p2)*u(m,n);    
    end
end
figure(3)
subplot(3,1,1)
mesh(-1:dx:1,0:dt:1,u);
title('\lambda = 2/5')
xlabel('x');
ylabel('t');
dt=lambda(2)*dx^2;
M=1/dt;N=2/dx;
u=zeros(M+1,N+1);
x=zeros(1,N+1);
for n=1:N+1    
    x(n)=-1+(n-1)*dx;    
    u(1,n)=x(n)^2;
end
u(:,1)=ones(1,M+1);
u(:,N+1)=ones(1,M+1);
for m=1:M    
    p1=dt/dx^2*((m-1)^2*dt^2+0.5)+dt/2/dx;    
    p2=dt/dx^2*((m-1)^2*dt^2+0.5)-dt/2/dx;    
    for n=2:N        
        u(m+1,n)=p1*u(m,n+1)+p2*u(m,n-1)+(1-p1-p2)*u(m,n);    
    end
end
figure(3)
subplot(3,1,2)
mesh(-1:dx:1,0:dt:1,u);
title('\lambda = 2/6')
xlabel('x');
ylabel('t');
dt=lambda(3)*dx^2;
M=1/dt;N=2/dx;
u=zeros(M+1,N+1);
x=zeros(1,N+1);
for n=1:N+1    
    x(n)=-1+(n-1)*dx;    
    u(1,n)=x(n)^2;
end
u(:,1)=ones(1,M+1);
u(:,N+1)=ones(1,M+1);
for m=1:M    
    p1=dt/dx^2*((m-1)^2*dt^2+0.5)+dt/2/dx;    
    p2=dt/dx^2*((m-1)^2*dt^2+0.5)-dt/2/dx;    
    for n=2:N        
        u(m+1,n)=p1*u(m,n+1)+p2*u(m,n-1)+(1-p1-p2)*u(m,n);    
    end
end
figure(3)
subplot(3,1,3)
mesh(-1:dx:1,0:dt:1,u);
title('\lambda = 2/7')
xlabel('x');
ylabel('t');
saveas(gcf,'q2bc','png')

