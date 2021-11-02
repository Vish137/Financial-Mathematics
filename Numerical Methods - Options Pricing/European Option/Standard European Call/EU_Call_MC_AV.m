%European Call Option - Monte Carlo - Antithetic Variable Method
S=100;K=100;r=0.07;sig=0.25;T=1;N=10^6;
Y=randn(1,N);
Z=zeros(1,N);%used for f(X) roughly speaking
V=zeros(1,N);%used for f(-X)
W=zeros(1,N);%(f(x)+f(-x))/2
for k=1:N    
    Z(k)=S*exp((r-0.5*sig^2)*T+sig*sqrt(T)*Y(k));    
    V(k)=S*exp((r-0.5*sig^2)*T-sig*sqrt(T)*Y(k));    
    W(k)=0.5*(max(Z(k)-K,0)+max(V(k)-K,0));
end
price=exp(-r*T)*sum(W)/N
histogram(W,100)