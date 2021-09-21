S=100;K=100;r=0.1;sig=0.3;T=1;
N=10000;
Y=randn(1,N);
Z=zeros(1,N);%for S_T
W=zeros(1,N);%for (S_T-K)^+
for k=1:N
    Z(k)=S*exp((r-0.5*sig^2)*T+sig*sqrt(T)*Y(k));
    W(k)=max(Z(k)-K,0);
end
U=zeros(1,N);%for \bar f
U(1)=W(1);
for k=2:N
    U(k)=U(k-1)+W(k);
end
for k=1:N
    U(k)=exp(-r*T)*U(k)/k;
end
plot(1:N,U)

hold on;


Y=randn(1,N);
Z=zeros(1,N);
V=zeros(1,N);
W=zeros(1,N);
for k=1:N
    Z(k)=S*exp((r-0.5*sig^2)*T+sig*sqrt(T)*Y(k));
    V(k)=S*exp((r-0.5*sig^2)*T-sig*sqrt(T)*Y(k));
    W(k)=0.5*(max(Z(k)-K,0)+max(V(k)-K,0));
end
U=zeros(1,N);
U(1)=W(1);
for k=2:N
    U(k)=U(k-1)+W(k);
end
for k=1:N
    U(k)=exp(-r*T)*U(k)/k;
end
plot(1:N,U)
hold on;

beta=-(log(K/S)-(r-0.5*sig^2)*T)\sig\sqrt(T);
N=10^4;
Y=randn(1,N);
W=zeros(1,N);%for g_beta(x)
for k=1:N
    W(k)=max(S*exp((r-0.5*sig^2)*T+sig*sqrt(T)*(Y(k)-beta))-K,0)*exp(beta*Y(k)-0.5*beta^2);
end
U=zeros(1,N);
U(1)=W(1);
for k=2:N
    U(k)=U(k-1)+W(k);
end
for k=1:N
    U(k)=exp(-r*T)*U(k)/k;
end
plot(1:N,U)
xlabel('Number of simulations')
ylabel('Price of call option')
title('Convergence diagram with MC and  ANTI and importance sampling')
hold on
