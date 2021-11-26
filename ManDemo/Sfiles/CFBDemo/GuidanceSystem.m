function [sys,x0,str,ts] = spacemodel(t,x,u,flag)

switch flag,
case 0,
    [sys,x0,str,ts]=mdlInitializeSizes;
case 1,
    sys=mdlDerivatives(t,x,u);
case 3,
    sys=mdlOutputs(t,x,u);
case {2,4,9}
    sys=[];
otherwise
    error(['Unhandled flag = ',num2str(flag)]);
end

function [sys,x0,str,ts]=mdlInitializeSizes
global Acoeff Bcoeff Ccoeff eta_d
sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 6;
sizes.NumInputs      = 0;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;
sys = simsizes(sizes);
x0  = [];
str = [];
ts  = [0 0];

A_D=pi/180;
eta_d(1)=4;
eta_d(2)=5;
eta_d(3)=10*A_D;
[Acoeff, Bcoeff, Ccoeff] = RefMod(0,200,[0;0;0],eta_d);
function sys=mdlOutputs(t,x,u)
global Acoeff Bcoeff Ccoeff eta_d
if(t<=200)
        eta_dr = [polyval(Acoeff,t);polyval(Bcoeff,t);polyval(Ccoeff,t)];
        d_eta_dr = [polyval([0 5*Acoeff(1) 4*Acoeff(2) 3*Acoeff(3) 2*Acoeff(4) Acoeff(5)],t);...
            polyval([0 5*Bcoeff(1) 4*Bcoeff(2) 3*Bcoeff(3) 2*Bcoeff(4) Bcoeff(5)],t);...
            polyval([0 5*Ccoeff(1) 4*Ccoeff(2) 3*Ccoeff(3) 2*Ccoeff(4) Ccoeff(5)],t)];
        d_d_eta_dr = [polyval([0 0 20*Acoeff(1) 12*Acoeff(2) 6*Acoeff(3) 2*Acoeff(4)],t);...
            polyval([0 0 20*Bcoeff(1) 12*Bcoeff(2) 6*Bcoeff(3) 2*Bcoeff(4)],t);...
            polyval([0 0 20*Ccoeff(1) 12*Ccoeff(2) 6*Ccoeff(3) 2*Ccoeff(4)],t)];
    else
        eta_dr = eta_d;
        d_eta_dr = [0;0;0];
        d_d_eta_dr = [0;0;0];
end
    

% qd=2*sin(0.2*pi*t);
% dqd=2*0.2*pi*cos(0.2*pi*t);
% ddqd=-2*0.2*0.2*pi*pi*sin(0.2*pi*t);

sys(1)=eta_dr(1);
sys(2)=eta_dr(2);
sys(3)=eta_dr(3);
sys(4)=d_eta_dr(1);
sys(5)=d_eta_dr(2);
sys(6)=d_eta_dr(3);
% sys(7)=d_d_eta_dr(1);
% sys(8)=d_d_eta_dr(2);
% sys(9)=d_d_eta_dr(3);
% sys(1)=2;
% sys(2)=1;
% sys(3)=5*pi/180;
% sys(4)=0;
% sys(5)=0;
% sys(6)=0;