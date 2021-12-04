function [t,u,v,r,x,y,psi,U] = pull_out(x,ui,U0,t_final,t_rudderexecute,h)
% ZIGZAG      [t,u,v,r,x,y,psi,U] = zigzag(ship,x,ui,t_final,t_rudderexecute,h,maneuver)
%             performs the zig-zag maneuver, see ExZigZag.m
%
% Inputs :
% 'ship'          = ship model. Compatible with the models under .../gnc/VesselModels/
% x               = initial state vector for ship model
% ui              = [delta,:] where delta=0 and the other values are non-zero if any
% t_final         = final simulation time
% t_rudderexecute = time control input is activated
% h               = sampling time
% maneuver        = [rudder angle, heading angle]. Default 20-20 deg that is: maneuver = [20, 20] 
%                    rudder is changed to maneuver(1) when heading angle is larger than maneuver(2)
%
% Outputs :
% t               = time vector
% u,v,r,x,y,psi,U = time series
%

if nargin~=6, error('number of inputs must be 6'); end
if t_final<t_rudderexecute, error('t_final must be larger than t_rudderexecute'); end
if nargin==5, h = 0.1; end

N = round(t_final/h);               % number of samples
xout = zeros(N+1,9);                % memory allocation
U = U0;

disp('Simulating...')

u_ship=ui;

for i=1:N+1
    time = (i-1)*h;
    psi = x(6)*180/pi;
    r   = x(3);

    if round(time) < t_rudderexecute
        u_ship(1) = 20*pi/180; 
    else
        u_ship(1) = 0;
    end
 
    [k11,U] = mariner(x,           u_ship);       % ship model
    [k12,U] = mariner(x+0.5*h*k11, u_ship);       % ship model
    [k13,U] = mariner(x+0.5*h*k12, u_ship);       % ship model
    [k14,U] = mariner(x+    h*k13, u_ship);       % ship model

    x = x + h * (k11 + 2*k12 + 2*k13 + k14) / 6;
    
    xout(i,:) = [time,x(1:6)',U,u_ship(1)]; 
end

for i=1:N+1
    time = (i-1)*h;
    psi = x(6)*180/pi;
    r   = x(3);
    
    if round(time) < t_rudderexecute
        u_ship(1) = -20*pi/180; 
    else
        u_ship(1) = 0;
    end
 
    [k11,U] = mariner(x,           u_ship);       % ship model
    [k12,U] = mariner(x+0.5*h*k11, u_ship);       % ship model
    [k13,U] = mariner(x+0.5*h*k12, u_ship);       % ship model
    [k14,U] = mariner(x+    h*k13, u_ship);       % ship model

    x = x + h * (k11 + 2*k12 + 2*k13 + k14) / 6;
    
    xout1(i,:) = [time,x(1:6)',U,u_ship(1)]; 
end


% time-series
t     = xout(:,1);
u     = xout(:,2); 
v     = xout(:,3);         
r1    = xout(:,4)*180/pi; 
r2    = xout1(:,4)*180/pi; 
x     = xout(:,5);
y     = xout(:,6);
psi   = xout(:,7)*180/pi;
U     = xout(:,8);
delta_c = xout(:,9)*180/pi;

% plots
figure(1)
subplot(111),
plot(t,r1,'b','linewidth',2)
hold on
plot(t,r2,'b','linewidth',2)
plot(t,0*zeros(length(t),1),'r','linewidth',2)
hold off
xlabel('time (s)'),title('yaw rate r (deg/s)'),grid


