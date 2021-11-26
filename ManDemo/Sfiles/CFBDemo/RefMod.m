function [Acoef, Bcoef, Ccoef] = RefMod(t0,t1,eta0,etad)

x0 = eta0(1);
y0 = eta0(2);
psi0 = eta0(3);
x1 = etad(1);
y1 = etad(2);
psi1 = etad(3);

x0d=0;x1d=0;y0d=0;y1d=0;psi0d=0;psi1d=0;
x0dd=0.0;x1dd=0;y0dd=0.0;y1dd=0;psi0dd=0;psi1dd=0;

AC = [t0^5 t0^4 t0^3 t0^2 t0 1;t1^5 t1^4 t1^3 t1^2 t1 1];
BC = [5*t0^4 4*t0^3 3*t0^2 2*t0 1 0;5*t1^4 4*t1^3 3*t1^2 2*t1 1 0];
CC = [20*t0^3 12*t0^2 6*t0 2 0 0;20*t1^3 12*t1^2 6*t1 2 0 0];
ZERO = zeros(2,6);

y = [x0 x1 y0 y1 psi0 psi1 x0d x1d y0d y1d psi0d psi1d x0dd x1dd y0dd y1dd psi0dd psi1dd]';

A = [AC ZERO ZERO;
        ZERO AC ZERO;
        ZERO ZERO AC;
        BC ZERO ZERO;
        ZERO BC ZERO;
        ZERO ZERO BC;
        CC ZERO ZERO;
        ZERO CC ZERO;
        ZERO ZERO CC;];
    
x=inv(A)*y;

Acoef = x(1:6);
Bcoef = x(7:12);
Ccoef = x(13:18);