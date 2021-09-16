S=100;K=100;r=0.07;sig=0.25;T=1;N=10^4;
Y=randn(1,N);
Z=zeros(1,N);%for S_T
W=zeros(1,N);%for (S_T-K)^+
for k=1:N    
    Z(k)=S*exp((r-0.5*sig^2)*T+sig*sqrt(T)*Y(k));    
    W(k)=max(-Z(k)+K,0);
end
U=zeros(1,N);%for \bar f
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
plot(1:N,U,'-b')
xlabel('Number of simulations')
ylabel('Price of European Put Option')
title('Monte Carlo Convergence Diagram')
grid minor
xlim([-100 inf])
saveas(gcf,'EU_Put_MC_Convergence','png')

figure(2)
plot(1:N,stderror,'-r','LineWidth',1.5)
xlabel('Number of simulations')
ylabel('Standard Error, $\frac{S}{\sqrt{n}}$','Interpreter','latex','FontSize',14)
title('Standard Error')
grid minor
xlim([-100 inf])
saveas(gcf,'EU_Put_MC_stderr','png')