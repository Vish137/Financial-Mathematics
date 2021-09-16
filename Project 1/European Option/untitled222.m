%explicit VS implicit (how dt affects them)

T=1;N=10;dt=0.001;M=T/dt;u0=1;u0p=2;
u=zeros(1,M+1);
u(1)=u0;
u(2)=u0+u0p*dt;
for i=3:M+1
    u(i)=(2+(1+(i-2)^2*dt^2)*dt-dt^2)/(1+(1+(i-2)^2*dt^2)*dt)*u(i-1)-1/(1+(1+(i-2)^2*dt^2)*dt)*u(i-2);
end
u(M+1)
plot(0:dt:1,u);
hold on;
T=1;N=10;dt=2^(-N);M=T/dt;u0=1;u0p=2;
w=zeros(1,M+1);
w(1)=u0;
w(2)=u0+u0p*dt;
for i=3:M+1
    w(i)=(2-(1+(i-2)^2*dt^2)*dt-dt^2)*w(i-1)+((1+(i-2)^2*dt^2)*dt-1)*w(i-2);
end
w(M+1)
plot(0:dt:1,w);
hold on;