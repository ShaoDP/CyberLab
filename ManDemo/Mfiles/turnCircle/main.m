clc;clear;close all;

t_final = 800;           % final simulation time (sec)
t_rudderexecute = 100;    % time rudder is executed (sec)
h = 0.1;                 % sampling time (sec)

disp('Turning Circle test for the Mariner class vessel')

% Mariner class cargo ship, cruise speed U0 = 7.7 m/s (see mariner.m)
x  = zeros(7,1);   % x  = [ u v r x y psi delta ]' (initial values)
ui = -35*pi/180;   % delta_c = -delta_R at time t = t_rudderexecute
U0 = 7.7175;

[t,u,v,r,x,y,psi,U] = turnCircle(x,ui,U0,t_final,t_rudderexecute,h,[20,20]);
