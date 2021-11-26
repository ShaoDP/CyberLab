close all;
 
L1=-5;
L2=5;
L=L2-L1;

T=L*1/1000;

xf=L1:T:L2;
figure(1);
for i=1:1:9
    gs=-0.5*(xf+2-(i-1)*0.5).^2;
    uf=exp(gs);
	hold on;
	plot(xf,uf);
end
xlabel('x');ylabel('Membership function');

figure(2);
subplot(311);
plot(t,x(:,1),'r',t,x(:,2),'b');
ylabel('tracking of x');
subplot(312);
plot(t,y(:,1),'r',t,y(:,2),'b');
ylabel('tracking of y');
subplot(313);
plot(t,psi(:,1),'r',t,psi(:,2),'b');
xlabel('time/s');ylabel('tracking of \psi');

figure(3);
subplot(311);
plot(t,tol(:,1),'k');hold on;
ylabel('\tau_x/N');
subplot(312);
plot(t,tol(:,2),'k');hold on;
ylabel('\tau_y/N');
subplot(313);
plot(t,tol(:,3),'k');hold on;
xlabel('time/s');ylabel('\tau_\psi/Nm');

figure(4);
subplot(311);
plot(t,u(:,1),'r',t,u(:,2),'b');
ylabel('tracking of x');
subplot(312);
plot(t,v(:,1),'r',t,v(:,2),'b');
ylabel('tracking of y');
subplot(313);
plot(t,r(:,1),'r',t,r(:,2),'b');
xlabel('time/s');ylabel('tracking of \psi');

figure(5);
xe = x(:,1)-x(:,2);
ye = x(:,1)-x(:,2);
psie = x(:,1)-x(:,2);

subplot(311);
plot(t,xe,'r');hold on;
ylabel('x_e/m');
subplot(312);
plot(t,ye,'r');hold on;
ylabel('x_e/m');
subplot(313);
plot(t,psie,'r');hold on;
xlabel('time/s');ylabel('\psi_e/deg');
