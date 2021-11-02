%OPTION PARAMETERS
T=1;
sigma =0.2;
r =0.07;
K =100;
q=2*r/ sigma ^2;

%NUMERICAL PARAMETERS
M =1000;
N =100;
dt= sigma ^2* T /(2* M);
dx =6* sigma * sqrt (T)/N;
lambda =dt/dx^2;
x= zeros (1,N+1);
for n=1:N+1
 x(n)= -3* sigma*sqrt (T)+(n-1)*dx;
end

u= zeros (1,N+1);
for n=1:N+1
    u(n)=max(exp ((q+1) /2*x(n))-exp ((q-1) /2*x(n)) ,0);
end

for m=2:M+1
   v= zeros (1,N+1);
   v(1) =0;
   v(N+1)=exp ((q+1)*x(N+1) /2+(q+1) ^2*(m-1)*dt /4) -exp ((q-1)*x(N+1) /2+(q-1) ^2*(m-1) *dt /4) ;
   for n=2:N
       v(n)= lambda *u(n-1) +(1-2* lambda )*u(n)+ lambda *u(n+1);
   end
    u=v;
end
S= zeros (1,N+1);
V= zeros (1,N+1);
for n=1:N+1
    S(n)=K*exp(x(n));
    V(n)=u(n)*K*exp (-(q -1) /2*x(n) -(q+1) ^2/4* sigma ^2*T/2) - (S(n)-K*exp(-r*T));
end

%PLOTTING
figure(1)
plot (S,V,'LineWidth',1,'Color','blue');
title ('European Put Option - Explicit Method')
xlabel ('Initial Price, S_0')
ylabel ('Option Price, $(V_0, S_0)$','Interpreter','latex')
grid minor
legend({'\lambda = 0.1389','Interpreter','latex'})
saveas(gcf,'Eu_Put_Explicit_Stable','png')

clear vars
%OPTION PARAMETERS
T=1;
sigma =0.2;
r =0.07;
K =100;
q=2*r/ sigma^2;

%NUMERICAL PARAMETERS
M = 1000;
N = 200;
dt= sigma^2* T /(2* M);
dx =6* sigma * sqrt (T)/N;
lambda =dt/dx^2;

%Create x-axis
x= zeros (1,N+1);
for n=1:N+1
 x(n)= -3* sigma*sqrt (T)+(n-1)*dx;
end

%Solve for payoff
u= zeros (1,N+1);
for n=1:N+1
    u(n)=max(exp ((q+1) /2*x(n))-exp ((q-1) /2*x(n)) ,0);
end

%Implement the explicit method 
for m=2:M+1
   v= zeros (1,N+1);
   v(1) =0;
   v(N+1)=exp ((q+1)*x(N+1) /2+(q+1) ^2*(m-1)*dt /4) -exp ((q-1)*x(N+1) /2+(q-1) ^2*(m-1) *dt /4) ;
   for n=2:N
       v(n)= lambda *u(n-1) +(1-2* lambda )*u(n)+ lambda *u(n+1);
   end
    u=v;
end

%Calculate option price
S= zeros (1,N+1);
V= zeros (1,N+1);
for n=1:N+1
    S(n)=K*exp(x(n));
    V(n)=u(n)*K*exp (-(q -1) /2*x(n) -(q+1) ^2/4* sigma ^2*T/2) - (S(n)-K*exp(-r*T));
end

%PLOTTING
figure(2)
plot (S,V,'LineWidth',1,'Color','blue');
title ('European Put Option - Explicit Method')
xlabel ('Initial Price, S_0')
ylabel ('Option Price, $(V_0, S_0)$','Interpreter','latex')
grid minor
legend({'\lambda = 0.5556','Interpreter','latex'})
saveas(gcf,'Eu_Put_Explicit_Unstable','png')

