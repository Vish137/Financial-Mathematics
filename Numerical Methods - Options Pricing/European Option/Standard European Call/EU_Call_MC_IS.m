%European Call Option - Importance Sampling
S=50;K=100;r=0.07;sig=0.25;T=1;
beta=-(log(K/S)-(r-0.5*sig^2)*T)\sig\sqrt(T);
N=10^4;
Y=randn(1,N);
W=zeros(1,N);%for g_beta(x)
for k=1:N    
    W(k)=max(S*exp((r-0.5*sig^2)*T+sig*sqrt(T)*(Y(k)-beta))-K,0)*exp(beta*Y(k)-0.5*beta^2);
end
price=exp(-r*T)*sum(W)/N
histogram(W,50)