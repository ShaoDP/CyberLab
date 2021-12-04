function [M,C,D] = clarke83(U,L,B,T,Cb,R66,xg,T_surge)
% [M,N] = clarke83(U,L,B,T,Cb,R66,xg,T_surge) computes the system matrices 
% of a linear maneuvering model based on Clarke et al. (1983). The  
% hydrodynamic derivatives are based on multiple  linear regression from two 
% sets of model tests. The first data set (Yv, Yr, Nv, Nr) is obtained from 
% rotating arm model experiments, while the second data set 
% (Yvdot, Yrdot, Nvdot, Nrdot, Yv, Yr, Nv, Nr) was obtained from a PMM  model.
% The surge model is approximated by Xudot = -0.1 * m. The time constant 
% in surge is optionally with default value T_surge = L such that 
% Xu = -(m-Xudot) / T_surge.
%
% Outputs: 3x3 model matrices M and N in surge, sway and yaw
%      .
%    M nu + N(U) nu = tau,     where N(U) = C(U) + D
% 
% corresponding to the linear maneuvering model
% 
%  (m - Xudot) udot - Xu u                            = (1-t) T
%  (m - Yvdot) vdot + (m - Yrdot)  rdot - Yv v - Yr r = Yd delta
%  (m - Yvdot) vdot + (Iz - Nrdot) rdot - Nv v - Nr r = Nd delta
%
% Note that the coefficients Yv, Yr, Nv and Nr in the N(U) matrix includes 
% linear damping D and the linearized Coriolis and centripetal matrix C(U).
%
% Inputs:
%
%  U:   speed (m/s)
%  L:   length (m)
%  B:   beam (m)
%  T:   draft (m)
%  Cb:  block coefficient (-), Cb = V / (L*B*T) where V is the displaced volume
%  R66: radius of gyration in yaw (smaller vessels R66 ≈ 0.25L, tankers R66 ≈ 0.27L)
%  xg:  x-coordinate of the CG
%  T_surge: (optionally) time constant in surge (default: T_surge = L)
%
% Reference:  CLARKE, D., GEDLING, P. and HINE. G. (1983). The application of 
% manoeuvring criteria in hull design using linear thory. Trans.  R. lnsm nav. 
% Archit.  125, 45-68. 
% 
% Author:    Thor I. Fossen
% Date:      22 Oct 2020
% Revisions: 14 Jun 2021 - Removed the C matrix and introduced N(U)

% Rigid body parameters
rho = 1025;                     % density of water
V = Cb * L * B * T;             % volume displacment
m = rho * V;                    % mass
Iz = m * R66^2 + m * xg^2;      % moment of inerta about the CO

MRB = [ m   0       0           % rigid-body inertia matrix
        0   m       m*xg
        0   m*xg    Iz      ];

% Nondimenisonal hydrodynamic derivatives in surge
Xudot = -0.1 * m;
if (nargin == 7)
    T_surge = L; 
end
UU = sqrt(U(1)^2 + U(2)^2) + 0.001; % avoid singularity for U = 0;
Xu = -((m-Xudot)/T_surge) / (0.5 * rho * L^2 * UU);
Xudot = Xudot / (0.5 * rho * L^3);

% Nondimenisonal hydrodynamic derivatives in sway and yaw
% from Clarke et al. (1983)
S = pi * (T/L)^2;                 % scale factor

Yvdot = -S * ( 1 + 0.16 * Cb * B/T - 5.1 * (B/L)^2 );
Yrdot = -S * ( 0.67 * B/L - 0.0033 * (B/T)^2 );
Nvdot = -S * ( 1.1 * B/L - 0.041 * (B/T) );
Nrdot = -S * ( 1/12 + 0.017 * Cb * (B/T) - 0.33 * (B/L) );
Yv = -S * ( 1 + 0.4 * Cb * (B/T) );
Yr = -S * ( -1/2 + 2.2 * (B/L) - 0.08 * (B/T) );
Nv = -S * ( 1/2 + 2.4 * (T/L) );
Nr = -S * ( 1/4 + 0.039 * (B/T) - 0.56 * (B/L) );

% Nondimenisonal hydrodynamic matrices 
MA_prime = [ -Xudot   0        0
              0      -Yvdot   -Yrdot
              0      -Nvdot   -Nrdot ];

N_prime = [ -Xu   0  0
             0  -Yv -Yr
             0  -Nv -Nr ];
 
% Dimensional model (Fossen 2021, Appendix D)   
T    = diag([1 1 1/L]);
Tinv = diag([1 1 L]);

MA = (0.5 * rho * L^3) * Tinv^2 * (T * MA_prime * Tinv);
D =  (0.5 * rho * L^2 * UU) * Tinv^2 * (T * N_prime * Tinv);
C = zeros(3, 3);
M = MRB + MA;       % system inertia matrix

 
 
