function [sys,x0,str,ts]=ssss(t,x,u,flag)
switch flag,
  case 0,
    [sys,x0,str,ts]=mdlInitializeSizes;
  case 1,
    sys=mdlDerivatives(t,x,u);
  case 3,
    sys=mdlOutputs(t,x,u);
  case {2, 4, 9}
    sys = [];
  otherwise
    error(['Unhandled flag = ',num2str(flag)]);
end

function [sys,x0,str,ts]=mdlInitializeSizes
global lamda1 lamda2 k

sizes = simsizes;
sizes.NumContStates  = 3+6;%zeta1 eta2c
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 6;
sizes.NumInputs      = 12;%etad etad_dot eta nu
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 0;
sys=simsizes(sizes);
x0=0.0*ones(9,1);
str=[];
ts=[];

lamda1=0.02*diag([30;30;30]);lamda2=0.2*diag([7e+1;7e1;7e1]);k=10.5;
function sys=mdlDerivatives(t,x,u)
global lamda1 k
x_p=u(7);y_p=u(8);psi=u(9);
dx=u(10);dy=u(11);dr=u(12);
J_eta = zeros(3,3);
J_eta(1,1) = cos(psi);
J_eta(1,2) = -sin(psi);
J_eta(2,1) = sin(psi);
J_eta(2,2) = cos(psi);
J_eta(3,3) = 1;
x1=[x_p;y_p;psi];
x2=J_eta*[dx;dy;dr];

r1=u(1);
dr1=u(4);
% ddr1=u(7);
r2=u(2);
dr2=u(5);
% ddr2=u(8);
r3=u(3);
dr3=u(6);
% ddr3=u(9);
etad=[r1;r2;r3];
detad=[dr1;dr2;dr3];
% ddyd=[ddr1;ddr2;ddr3];

%% modified
for i=1:1:3
    Fai1(i)=x(i);
end
for i=1:1:3
    Fai2(i)=x(i+3);
end
for i=1:1:3
    Fai3(i)=x(i+6);
end

u1(1)=exp(-1/2*((x(1)+1.25)/0.6)^2);
u1(2)=exp(-1/2*(x(1)/0.6)^2);
u1(3)=exp(-1/2*((x(1)-1.25)/0.6)^2);

u2(1)=exp(-1/2*((x(2)+1.25)/0.6)^2);
u2(2)=exp(-1/2*(x(2)/0.6)^2);
u2(3)=exp(-1/2*((x(2)-1.25)/0.6)^2);

u3(1)=exp(-1/2*((x(3)+1.25)/0.6)^2);
u3(2)=exp(-1/2*(x(3)/0.6)^2);
u3(3)=exp(-1/2*((x(3)-1.25)/0.6)^2);

u4(1)=exp(-1/2*((x(4)+1.25)/0.6)^2);
u4(2)=exp(-1/2*(x(4)/0.6)^2);
u4(3)=exp(-1/2*((x(4)-1.25)/0.6)^2);

sum1=0;
for i=1:1:3
        fs1=u1(i)*u2(i)*u3(i);
        sum1=sum1+fs1;
        P1(i)=fs1/(sum1+0.01);
end
P2=P1;P3=P2;


z1=x1-etad;
eta2c0=-lamda1*z1+detad;
% alfa1=-lamda1*z1+dyd;
% z2=x2-alfa1;
zeta1=[x(1);x(2);x(3)];
% eta2c=[x(4);x(5);x(6);x(7);x(8);x(9)];
eta2c1=[x(4);x(5);x(6)];
eta2c2=[x(7);x(8);x(9)];
zeta1_dot=-lamda1*zeta1+(eta2c1-eta2c0);

omega_n=diag([20 20 20]);zeta=diag([5 5 5]);
eta2c_dot=[eta2c2;
            -omega_n*omega_n*eta2c1-2*omega_n*zeta*eta2c2+omega_n*omega_n*eta2c0];

for i=1:1:3
    sys(i)=zeta1_dot(i)*1;
end
for i=1:1:6
    sys(3+i)=eta2c_dot(i)*1;
end

function sys=mdlOutputs(t,x,u)
global lamda1 lamda2
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
psi=u(9);
% nu=[u(7);u(8);u(9)];

J_eta = zeros(3,3);
J_eta(1,1) = cos(psi);
J_eta(1,2) = -sin(psi);
J_eta(2,1) = sin(psi);
J_eta(2,2) = cos(psi);
J_eta(3,3) = 1;
nu=[u(10);u(11);u(12)];%inv(J_eta)*
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
M_eta=inv(J_eta)'*M*inv(J_eta);
%% Coriolis-centripetal matrix
C_nu = zeros(3,3);
C_nu(1,3) = -m*(x_g*r+v)+Y_dv*v+Y_dr*r;
C_nu(2,3) = m*u_s-X_du*u_s;
C_nu(3,1) = -C_nu(1,3);
C_nu(3,2) = -C_nu(2,3);
C_eta = inv(J_eta)'*(C_nu-M*inv(J_eta)*d_J)*inv(J_eta);
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
D_eta=inv(J_eta)'*D_nu*inv(J_eta);

% d=0*[0.25*sin(t) 0.25*sin(t) 0.25*sin(t)]';

x_p=u(7);y_p=u(8);psi=u(9);
dx=u(10);dy=u(11);dr=u(12);
J_eta = zeros(3,3);
J_eta(1,1) = cos(psi);
J_eta(1,2) = -sin(psi);
J_eta(2,1) = sin(psi);
J_eta(2,2) = cos(psi);
J_eta(3,3) = 1;
x1=[x_p;y_p;psi];
x2=J_eta*[dx;dy;dr];
eta=x1;nu=[dx;dy;dr];

r1=u(1);
dr1=u(4);
% ddr1=u(7);
r2=u(2);
dr2=u(5);
% ddr2=u(8);
r3=u(3);
dr3=u(6);
% ddr3=u(9);
etad=[r1;r2;r3];
detad=[dr1;dr2;dr3];
% ddyd=[ddr1;ddr2;ddr3];

zeta1=[x(1);x(2);x(3)];
eta2c=[x(4);x(5);x(6);x(7);x(8);x(9)];
eta2c1=[x(4);x(5);x(6)];
eta2c2=[x(7);x(8);x(9)];

z1=x1-etad;
eta2c0=-lamda1*z1+detad;
z2=x2-eta2c1;
v1=z1-zeta1;
f=inv(M_eta)*C_eta*x2+inv(M_eta)*D_eta*x2;
tauc0=M_eta*(-lamda2*z2+eta2c2+f-v1);

tol=tauc0;
tol=J_eta*tol;

sys(1)=tol(1);
sys(2)=tol(2);
sys(3)=tol(3);
sys(4)=eta2c2(1);
sys(5)=eta2c2(2);
sys(6)=eta2c2(3);