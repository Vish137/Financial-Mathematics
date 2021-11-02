%European Call Option - Monte Carlo - Algorithm
S=100;K=100;r=0.07;sig=0.25;T=1;N=10^6;
Y=randn(1,N);
Z=zeros(1,N);%will be used for price for S_T
W=zeros(1,N);
for k=1:N    
    Z(k)=S*exp((r-0.5*sig^2)*T+sig*sqrt(T)*Y(k));%S_T for each simulation    
    W(k)=max(Z(k)-K,0);%(S_T-K)^+
end
price=exp(-r*T)*sum(W)/N
histogram(W,100)