clc;clear;close all;

t_final = 600;           % final simulation time (sec)
t_rudderexecute = 10;    % time rudder is executed (sec)
h = 0.1;                 % sampling time (sec)

N = round(t_final/h);               % number of samples
xout = zeros(N+1,17);                % memory allocation

x = zeros(12,1);
x(9) = 100;

disp('Simulating...')

for i=1:N+1
    time = (i-1)*h;

    %  ui = [ delta_r delta_s n ]'  where
    %
    %    delta_r:   rudder angle (rad)
    %    delta_s:   aft stern plane (rad)
    %    n:         propeller revolution (rpm)
    delta_r =  30*pi/180;
    delta_s =  50*pi/180;
    n = 1500;

    ui = [delta_r delta_s n]';

% state vector: x = [ u v w p q r x y z phi theta psi ]' and speed U in m/s
    [k11,U] = remus100(x,           ui);       % ship model
    [k12,U] = remus100(x+0.5*h*k11, ui);       % ship model
    [k13,U] = remus100(x+0.5*h*k12, ui);       % ship model
    [k14,U] = remus100(x+    h*k13, ui);       % ship model

    x = x + h * (k11 + 2*k12 + 2*k13 + k14) / 6;
    x(12) = pipi(x(12));

    xout(i,:) = [time,x(1:12)',U,ui']; 
end

t     = xout(:,1);
u     = xout(:,2); 
v     = xout(:,3);         
w     = xout(:,4); 
p     = xout(:,5);
q     = xout(:,6);
r     = xout(:,7)*180/pi;
x     = xout(:,8);
y     = xout(:,9);
z     = xout(:,10);
phi   = xout(:,11)*180/pi;
theta = xout(:,12)*180/pi;
psi   = xout(:,13)*180/pi;
U     = xout(:,14);

% plots
figure(1)
plot3(x,y,z,'linewidth',2),grid,axis('equal'),xlabel('x-position'),ylabel('y-position'),ylabel('z-position')
title('position')
figure(2)
subplot(311),plot(t,u,'linewidth',2),grid,ylabel('speed u (m/s)')
subplot(312),plot(t,v,'linewidth',2),grid,ylabel('speed v (m/s)')
subplot(313),plot(t,w,'linewidth',2),grid,ylabel('speed w (m/s)')
figure(3)
subplot(311),plot(t,p,'linewidth',2),grid,ylabel('speed p (deg/s)')
subplot(312),plot(t,q,'linewidth',2),grid,ylabel('speed q (deg/s)')
subplot(313),plot(t,r,'linewidth',2),grid,ylabel('speed r (deg/s)')
figure(4)
subplot(311),plot(t,x,'linewidth',2),grid,ylabel('pos x (m)')
subplot(312),plot(t,y,'linewidth',2),grid,ylabel('pos y (m)')
subplot(313),plot(t,z,'linewidth',2),grid,ylabel('pos z (m)')
figure(5)
subplot(311),plot(t,phi,  'linewidth',2),grid,ylabel('angle \phi (deg)')
subplot(312),plot(t,theta,'linewidth',2),grid,ylabel('angle \theta (deg)')
subplot(313),plot(t,psi,  'linewidth',2),grid,ylabel('angle \psi (deg)')