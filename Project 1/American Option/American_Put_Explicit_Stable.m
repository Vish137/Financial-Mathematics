%American Call Option - Explicit Finite Difference Scheme

sigma=0.2;r=0.07;K=100;T=1;q=2*r/sigma^2;M=30000;N=1000;
dt=sigma^2*T/2/M;
dx=6*sigma*sqrt(T)/N;
lambda=dt/dx^2;
x=zeros(1,N+1);
lambda

for n=1:N+1    
    x(n)=-3*sigma*sqrt(T)+(n-1)*dx;
end
u=zeros(1,N+1);
for n=1:N+1    
    u(n)=max(exp((q-1)/2*x(n))-exp((q+1)/2*x(n)),0);
end

%t-loop%
for m = 2:M+1
    w = zeros(1,N+1);
    %Compute Boundary Values%
    w(1) = exp((q-1)*x(1)/2+(q-1)^2*(m-1)*dt/4)-exp((q+1)*x(1)/2+(q+1)^2*(m-1)*dt/4);
    w(N+1) = 0;
    for n = 2:N
        w(n)=lambda*u(n-1)+(1-2*lambda)*u(n)+lambda*u(n+1);
    end
    g = zeros(1,N+1);
    for n = 1:N+1 %Step 2, different from European option
        g(n)=exp((q+1)^2/4*(m-1)*dt)*max(exp((q-1)/2*x(n))-exp((q+1)/2*x(n)),0);
    end
    u = max(w,g);
    n = 1;
    while u(n) == g(n)
        n = n+1;
    end
end
S=zeros(1,N+1);
V=zeros(1,N+1);
for n = 1:N+1
    S(n) = K*exp(x(n));
    V(n) = u(n)*K*exp(-(q-1)/2*x(n)-(q+1)^2/4*sigma^2*T/2);
end
plot(S,V,'red');
title ('American Put Option - Explicit Method')
xlabel ('Initial Price, S_0')
ylabel ('Option Price, $(V_0, S_0)$','Interpreter','latex')
grid minor
legend({'\lambda = 0.4630','Interpreter','latex'})