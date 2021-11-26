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
sizes = simsizes;
sizes.NumContStates  = 6;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 3;
sizes.NumInputs      = 9;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;
sys = simsizes(sizes);
x0=[0 0 0 0 0 0];
str=[];
ts=[0 0];
function sys=mdlDerivatives(t,x,u)
Gama1=12;Gama2=12;Gama3=12;
k=10.20;
w =[u(1) u(2) u(3)];
r =[u(4) u(5) u(6)];
ut=[u(7) u(8) u(9)];

if w(1)<0
   Xp1=0;Xn1=1;
else
   Xp1=1;Xn1=0;
end
if w(2)<0
   Xp2=0;Xn2=1;
else
   Xp2=1;Xn2=0;
end
if w(3)<0
   Xp3=0;Xn3=1;
else
   Xp3=1;Xn3=0;
end
X1=[Xp1;-Xn1];
X2=[Xp2;-Xn2];
X3=[Xp3;-Xn3];

d1p=x(1);d1n=x(2);
d2p=x(3);d2n=x(4);
d3p=x(5);d3n=x(6);

D1=[d1p d1n]';
D2=[d2p d2n]';
D3=[d3p d3n]';

dD1=Gama1*X1*r(1)-k*Gama1*D1*abs(r(1));
dD1p=dD1(1);
dD1n=dD1(2);

dD2=Gama2*X2*r(2)-k*Gama2*D2*abs(r(2));
dD2p=dD2(1);
dD2n=dD2(2);

dD3=Gama3*X3*r(3)-k*Gama3*D3*abs(r(3));
dD3p=dD3(1);
dD3n=dD3(2);
sys(1)=dD1p;
sys(2)=dD1n;
sys(3)=dD2p;
sys(4)=dD2n;
sys(5)=dD3p;
sys(6)=dD3n;
function sys=mdlOutputs(t,x,u)
w =[u(1) u(2) u(3)];
r =[u(4) u(5) u(6)];
ut=[u(7) u(8) u(9)];

if w(1)<0
   Xp1=0;Xn1=1;
else
   Xp1=1;Xn1=0;
end
if w(2)<0
   Xp2=0;Xn2=1;
else
   Xp2=1;Xn2=0;
end
if w(3)<0
   Xp3=0;Xn3=1;
else
   Xp3=1;Xn3=0;
end
X1=[Xp1;-Xn1];
X2=[Xp2;-Xn2];
X3=[Xp3;-Xn3];

d1p=x(1);d1n=x(2);
d2p=x(3);d2n=x(4);
d3p=x(5);d3n=x(6);

D1=[d1p d1n]';
D2=[d2p d2n]';
D3=[d3p d3n]';
D1=[6e4 5e4]';
D2=[6e4 5e4]';
D3=[6e4 5e4]';

sys(1)=D1'*X1;
sys(2)=D2'*X2;
sys(3)=D3'*X3;