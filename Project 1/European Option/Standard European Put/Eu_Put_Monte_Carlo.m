%OPTION PARAMETERS%
S=100;K=100;r=0.07;sig=0.25;T=1;N=10^4;

%STANDARD MONTE CARLO APPROACH
Y=randn(1,N);
Z=zeros(1,N);%for S_T
W=zeros(1,N);%for (S_T-K)^+
for k=1:N    
    Z(k)=S*exp((r-0.5*sig^2)*T+sig*sqrt(T)*Y(k));    
    W(k)=max(-Z(k)+K,0);
end
U=zeros(1,N);
U(1)=W(1);
for k=2:N    
    U(k)=U(k-1)+W(k);
end
stderror=zeros(1,N);
for k=1:N    
    U(k)=exp(-r*T)*U(k)/k;
    stderror(k)=std(W)/sqrt(k); %Check the standard error of the payoff at each timestep
end
figure(1)
plot(1:N,U,'-b','DisplayName','Std. Monte Carlo')
xlabel('Number of simulations')
ylabel('Price of European Put Option')
title('Monte Carlo Convergence Diagram')
%grid minor
hold on

figure(2)
plot(1:N,stderror,'-r','LineWidth',1.5,'DisplayName','Std. Monte Carlo')
xlabel('Number of simulations')
ylabel('Standard Error, $\frac{S}{\sqrt{n}}$','Interpreter','latex','FontSize',14)
title('Standard Error')
%grid minor
xlim([-100 inf])
hold on

clear all
%OPTION PARAMETERS%
S=100;K=100;r=0.07;sigma=0.25;T=1;N=10^4;
%IMPORTANCE SAMPLING MONTE CARLO APPROACH
beta=-log(K/S)-(r-0.5*sigma^2)*T\sigma\sqrt(T);
Y=randn(1,N);
W=zeros(1,N);
for k=1:N
    W(k)=max(S*exp((r-0.5*sigma^2)*T+sigma*sqrt(T)*(Y(k)-beta))-K,0)*exp(beta*Y(k)-0.5*beta^2);
end 
U=zeros(1,N);
U(1)=W(1);
for k=2:N
    U(k)=U(k-1)+W(k);
end
stderror=zeros(1,N);
for k=1:N
    U(k)=exp(-r*T)*U(k)/k;
    stderror(k)=std(W)/sqrt(k);
end
figure(1)
plot(1:N,U,'DisplayName','Imp. Sam. Monte Carlo')
hold on

figure(2)
plot(1:N,stderror,'DisplayName','Imp. Sam. Monte Carlo')
hold on

clear all
%OPTION PARAMETERS%
S=100;K=100;r=0.07;sigma=0.25;T=1;N=10^4;
%VARIANCE REDUCTION MONTE CARLO APPROACH
N=10^4;
Y=randn(1,N); %brownian motion random term
Z=zeros(1,N); %to compute Y2
V=zeros(1,N); %to compute Y1
W=zeros(1,N); %to compute the payoff as an unbiased estimator

for k=1:N
    Z(k)=S*exp((r-0.5*sigma^2)*T+sigma*sqrt(T)*Y(k));
    V(k)=S*exp((r-0.5*sigma^2)*T+sigma*sqrt(T)*Y(k));
    W(k)=0.5*(max(K-Z(k),0)+max(K-V(k),0));
end 
U=zeros(1,N);
U(1)=W(1);
for k=2:N
    U(k)=U(k-1)+W(k);
end 
stderror=zeros(1,N);
for k=1:N
    U(k)=exp(-r*T)*U(k)/k;
    stderror(k)=std(W)/sqrt(k);
end 

figure(1)
plot(1:N,U,'DisplayName','Antithetic var')
hold on
figure(1)
yline(BSExact(S,K,r,sigma,T),'-.r','','LineWidth',0.5,'DisplayName','Black-Scholes Exact')
legend
hold off
saveas(gcf,'EU_Put_MC_Convergence','png')

figure(2)
plot(1:N,stderror,'--b','DisplayName','Antithetic var')
legend
hold off
saveas(gcf,'EU_Put_MC_stderr','png')