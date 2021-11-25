function x_dot = shipPlant2(x, tau)
%% vessel 2 model
% mass
m11  = 1.44846E+7;
m22  = 2.54376E+7;
m33  = 5.0972430E+10;

% drag
Xuu = 3.6000E+04;
Yvv = 4.3700E+05;
Nrr = 4.1623750E+010;
M = [m11  0.0  0.0;
     0.0  m22  0.0;
     0.0  0.0  m33]  * 10;
D =  [Xuu  0.0  0.0;
      0.0  Yvv  0.0;
      0.0  0.0  Nrr] * 20;
psi = x(3);
J_eta = zeros(3,3);
J_eta(1,1) = cos(psi);
J_eta(1,2) = -sin(psi);
J_eta(2,1) = sin(psi);
J_eta(2,2) = cos(psi);
J_eta(3,3) = 1;

nu=[x(4);x(5);x(6)];
eta_dot = J_eta*nu;
nu_dot  = M \ (tau - D*nu.*abs(nu));
x_dot = [
    eta_dot;
    nu_dot
    ];
end

