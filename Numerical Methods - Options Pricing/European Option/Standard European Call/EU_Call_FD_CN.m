%Finite Difference - European Call - Crank-Nicolson Method
sigma=0.2;r=0.07;K=100;T=1;q=2*r/sigma^2;M=1000;N=200;
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
    b(1)=b(1)+lambda/2*w(1);b(N-1)=b(N-1)+lambda/2*w(N+1);    
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
S=zeros(1,N+1);
V=zeros(1,N+1);
for n=1:N+1    
    S(n)=K*exp(x(n));    
    V(n)=u(n)*K*exp(-(q-1)/2*x(n)-(q+1)^2/4*sigma^2*T/2);
end
plot(S,V);