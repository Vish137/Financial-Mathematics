%TEST 1 - Function of Strike%
%Option Parameters%
N=100;
K=100;
r=0.07;
sigma=1.25;

%Grid%
T=1;
m=zeros(1,N);
n=zeros(1,N);
for i=1:N
    [m(i),n(i)] = BSExact(S,i,r,sigma,T);
end 
%Plotting both option prices
figure(1);
plot(m,'-r','LineWidth',1.5);
hold on;
plot(n,'-b','LineWidth',1.5);
grid minor
title('Black Scholes Pricing (r = 0.07, \sigma = 1.25)')
xlabel('Strike Price, K')
ylabel('Option Price, C')
legend({'put','call'})
saveas(gcf,'BS_option_test1','png')
hold off;

%Plotting call price 
figure(2);
plot(n,'-b','LineWidth',1.5);
grid minor
title('Black Scholes - Call Option (r = 0.07, \sigma = 1.25)')
xlabel('Strike Price, K')
ylabel('Option Price, C')
saveas(gcf,'BS_call_option_test1','png')

%Plotting put price 
figure(3);
plot(m,'-r','LineWidth',1.5);
grid minor
title('Black Scholes - Put Option (r = 0.07, \sigma = 1.25)')
xlabel('Strike Price, K')
ylabel('Option Price, C')
saveas(gcf,'BS_put_option_test1','png')
