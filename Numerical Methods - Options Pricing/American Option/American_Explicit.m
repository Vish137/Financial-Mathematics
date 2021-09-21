%OPTION PARAMETERS%
r=0.05;
sigma=0.25;
Nt=2000;
Ns=200;
Smax=200;
Smin=0;
T=1;
K=100;

dt = T/Nt; %time steps
ds=(Smax-Smin)/Ns; %price steps
V(1:Ns+1,1+Ns+1)=0;
S=Smin+(0:Ns)*ds;
t=(0:Nt)*dt;
V(1:Ns+1,1)=max(S-K,0);

%Boundary conditions
V(1,1:Nt+1)=0;
V(Ns+1,1:Nt+1)=Smax*exp(-t)-K*exp(-r*t);

%Explicit algorithm
for j=1:Nt
    for n=2:Ns
        V(n,j+1)=0.5*dt*((sigma^2)*(n^2)-r*n)*V(n-1,j)+(1-dt*((sigma^2)*(n^2)+r))*V(n,j)+0.5*dt*((sigma^2)*(n^2)+(r)*n)*V(n+1,j);
    end
    V(1:Ns+1,j+1)=max(V(1:Ns+1,j+1),V(1:Ns+1,1)); %check if striking out is better
end
plot(S,V(:,1),'-k','LineWidth',1.5)
hold on
plot(S,V(:,Nt/2),'-.r','LineWidth',1)
hold on
plot(S,V(:,Nt+1),'-b','LineWidth',1.5)
hold off
title('American Call Option (Explicit Method)')
xlabel('S')
ylabel('V(S,t)')
legend({'t=1','t=N/2','t=N+1'},'Location','northwest')
saveas(gcf,'american_explicit','png')
