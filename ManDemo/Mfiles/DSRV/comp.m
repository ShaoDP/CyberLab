clear; close all; clc;

%% sample data
h = 0.1;
tSpan = 600;
Num = tSpan / h;
to_deg = 180.0/pi;

%% save data
t_data = zeros(Num, 1);
x1_data = zeros(Num, 6);
x2_data = x1_data;

%% vessel 1
% x = [ w q x z theta ]' for a deep submergence rescue vehicle (DSRV) 
x1 = zeros(5, 1);

%% vessel2
x2 = zeros(5, 1);

%% input
delta = 10.0 / to_deg;

%% RK
for i = 1:Num
    % model 1
%     [xdot,U] = DSRV(x,u)
    [k11, U] = DSRV(x1, delta);
    [k12, U] = DSRV(x1+0.5*h*k11, delta);
    [k13, U] = DSRV(x1+0.5*h*k12, delta);
    [k14, U] = DSRV(x1+h*k13, delta);
    
    x1 = x1 + h * (k11 + 2*k12 + 2*k13 + k14) / 6;
    
    % model 2
%     [xdot,U] = DSRV(x,u)
    [k11, U] = DSRV1(x2, delta);
    [k12, U] = DSRV1(x2+0.5*h*k11, delta);
    [k13, U] = DSRV1(x2+0.5*h*k12, delta);
    [k14, U] = DSRV1(x2+h*k13, delta);
    
    x2 = x2 + h * (k11 + 2*k12 + 2*k13 + k14) / 6;
    % save
    t_data(i) = i * h;
    x1_data(i, :) = [x1' U];
    x2_data(i, :) = [x2' U];
end

%% plot curves
figure
plot(t_data, x1_data(:, 1), t_data, x2_data(:, 1));
title("speed w");xlabel("t/s");ylabel("w/[m/s]")
legend("model 1", "model 2");

figure
plot(t_data, x1_data(:, 2)*to_deg, t_data, x2_data(:, 2)*to_deg);
title("speed q");xlabel("t/s");ylabel("q/[deg/s]")
legend("model 1", "model 2");

figure
plot(t_data, x1_data(:, 3), t_data, x2_data(:, 3));
title("pos x");xlabel("t/s");ylabel("x/m")
legend("model 1", "model 2");

figure
plot(t_data, x1_data(:, 4), t_data, x2_data(:, 4));
title("pos z");xlabel("t/s");ylabel("z/m")
legend("model 1", "model 2");

figure
plot(t_data, x1_data(:, 5)*to_deg, t_data, x2_data(:, 5)*to_deg);
title("angle \theta");xlabel("t/s");ylabel("\theta/[deg]")
legend("model 1", "model 2");

figure
plot(t_data, x1_data(:, 6), t_data, x2_data(:, 6));
title("speed U");xlabel("t/s");ylabel("U/[m/s]")
legend("model 1", "model 2");

figure
plot(x1_data(:, 3), x1_data(:, 4), x2_data(:, 3), x2_data(:, 4)); 
title("position");xlabel("x/m");ylabel("z/m")
legend("model 1", "model 2");
