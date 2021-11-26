%double links adptive fuzzy control by backstepping
function [sys,x0,str,ts]= Shipplant1 (t,x,u,flag)
switch flag,
  case 0,
    [sys,x0,str,ts]=mdlInitializeSizes;
  case 1,
    sys=mdlDerivatives(t,x,u);
  case 3,
    sys=mdlOutputs(t,x,u);
  case {2, 4, 9 }
    sys = [];
  otherwise
    error(['Unhandled flag = ',num2str(flag)]);
end
function [sys,x0,str,ts]=mdlInitializeSizes
sizes = simsizes;
sizes.NumContStates  = 6;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 6;
sizes.NumInputs      = 3;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 0;
sys=simsizes(sizes);
x0=[0.0 0.0 0.0 0 0 0];
str=[];
ts=[];
function sys=mdlDerivatives(t,x,u)
%   x=[x,y,psi,u,v,r]'
%%              ship data            %%%%%%%%
m = 4.3e5*4.5833^3;         %4.3e5kg
x_g = -0.0137*4.5833;       %-0.0137m
Iz = 5.4e6*4.5833^5;        %kg.s^2.m^2
rho = 1025;                 %kg.m^-3
L = 48*4.5833;              %48m

rho_a = 1.2;    %空气密度
rho_w = rho;%水密度
rho_s = 7.85e3; %钢密度
g = 9.81;       %重力加速度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  zero dimensional hydrodynamic coefficients  %%%
X_du_neg = -9.469e-4;
Y_dv_neg = -8.316e-3;
N_dr_neg = -3.647e-4;
Y_dr_neg = -8.078e-4;
X_auu_neg = -1.2e-3;
Y_avv_neg = -3.264e-3;
Y_arv_neg = -2.251e-3;
Y_avr_neg = -2.251e-3;
Y_arr_neg = -4.63e-3;
N_avv_neg = 5.3e-4;
N_arv_neg = -1.427e-4;
N_avr_neg = -1.427e-4;
N_arr_neg = -4.63e-3;

%% dimensional hydrodynamic coefficients
X_du = X_du_neg*0.5*rho*L^3;
Y_dv = Y_dv_neg*0.5*rho*L^3;
N_dr = N_dr_neg*0.5*rho*L^5;
Y_dr = Y_dr_neg*0.5*rho*L^4;
X_auu = X_auu_neg*0.5*rho*L^2;
Y_avv = Y_avv_neg*0.5*rho*L^2;
Y_arv = Y_arv_neg*0.5*rho*L^3;
Y_avr = Y_avr_neg*0.5*rho*L^3;
Y_arr = Y_arr_neg*0.5*rho*L^4;
N_avv = N_avv_neg*0.5*rho*L^3;
N_arv = N_arv_neg*0.5*rho*L^4;
N_avr = N_avr_neg*0.5*rho*L^4;
N_arr = N_arr_neg*0.5*rho*L^5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% tranformation matrix
psi = x(3); 

J_eta = zeros(3,3);
J_eta(1,1) = cos(psi);
J_eta(1,2) = -sin(psi);
J_eta(2,1) = sin(psi);
J_eta(2,2) = cos(psi);
J_eta(3,3) = 1;
nu=[x(4);x(5);x(6)];%inv(J_eta)*
u_s = nu(1); v = nu(2); r = nu(3);
S_rou=zeros(3,3);
S_rou(1,2)=-r;S_rou(2,1)=r;
d_J=J_eta*S_rou;
 
%% system inertia matrix
M = zeros(3,3);
M(1,1) = m-X_du;%m-X_du
M(2,2) = m-Y_dv;%m-Y_dv
M(2,3) = m*x_g-Y_dr;
M(3,2) = M(2,3);
M(3,3) = Iz-N_dr;

% change M
M = 10 * M;

M_eta=inv(M)'*M*inv(M);
%% Coriolis-centripetal matrix
C_nu = zeros(3,3);
C_nu(1,3) = -m*(x_g*r+v)+Y_dv*v+Y_dr*r;
C_nu(2,3) = m*u_s-X_du*u_s;
C_nu(3,1) = -C_nu(1,3);
C_nu(3,2) = -C_nu(2,3);
C_eta = inv(M)'*(C_nu-M*inv(J_eta)*d_J)*inv(M);
%% damping matrix
D_nu = zeros(3,3);
D_nu(1,1) = -X_auu*abs(u_s);
D_nu(2,2) = -Y_avv*abs(v)-Y_arv*abs(r);
D_nu(2,3) = -Y_avr*abs(v)-Y_arr*abs(r);
D_nu(3,2) = -N_avv*abs(v)-N_arv*abs(r);
D_nu(3,3) = -N_avr*abs(v)-N_arr*abs(r);
D =1*3.0e+05*[0.50242 0 0; 
    0 2.7229   -43.933;
    0 -43.933 4189.4];
D_nu=(D_nu+D);
D_eta=inv(M)'*D_nu*inv(M);

d=0*[0.25e5*sin(t) 0.25e5*sin(t) 0.25e6*sin(t)]';

tol=u;
d_position=J_eta*nu;
sys(1)=d_position(1);
sys(2)=d_position(2);
sys(3)=d_position(3);
% S=-inv(M_eta)*(C_eta+D_eta)*J_eta*nu+inv(M_eta)*(u-d);
S=-inv(M)*(C_nu+D_nu)*nu+inv(M)*(tol-d);
sys(4)=S(1);
sys(5)=S(2);
sys(6)=S(3);
function sys=mdlOutputs(t,x,u)
sys(1)=x(1);
sys(2)=x(2);
sys(3)=x(3);
sys(4)=x(4);
sys(5)=x(5);
sys(6)=x(6);