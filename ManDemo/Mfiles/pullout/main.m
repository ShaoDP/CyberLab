clc;clear;close all;

t_final = 1200;           % final simulation time (sec)
t_rudderexecute = 600;    % time rudder is executed (sec)
h = 0.1;                 % sampling time (sec)

disp('Pullout maneuver for the Mariner class vessel')

% 20-20 zigzag maneuver for the Mariner class cargo ship
% cruise speed U0 = 7.7 m/s (see mariner.m)
x  = zeros(7,1);   % x  = [ u v r x y psi delta ]' (initial values)
ui = 0;            % delta_c = 0  for time t < t_rudderexecute
U0 = 7.7175;

[t,u,v,r,x,y,psi,U] = pull_out(x,ui,U0,t_final,t_rudderexecute,h);
