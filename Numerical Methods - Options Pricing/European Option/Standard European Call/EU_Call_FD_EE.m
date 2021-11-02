%EU Call - Finite Difference - Explicit Euler Scheme - Black Scholes PDE
sigma=0.2;r=0.07;K=100;T=1;q=2*r/sigma^2;M=1000;N=100;
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
    w=zeros(1,N+1);%for u(t+dt)    
    w(1)=0;    
    w(N+1)=exp((q+1)*x(N+1)/2+(q+1)^2*(m-1)*dt/4)-exp((q-1)*x(N+1)/2+(q-1)^2*(m-1)*dt/4);    
    for n=2:N        
        w(n)=lambda*u(n-1)+(1-2*lambda)*u(n)+lambda*u(n+1);    
    end%get u(t+dt)    
    u=w;
end
S=zeros(1,N+1);
V=zeros(1,N+1);
for n=1:N+1    
    S(n)=K*exp(x(n));    
    V(n)=u(n)*K*exp(-(q-1)/2*x(n)-(q+1)^2/4*sigma^2*T/2);
end
plot(S,V);
