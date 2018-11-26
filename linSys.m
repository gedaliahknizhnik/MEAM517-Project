function [A,B] = linSys(w_star, u_star, params)
% LINSYS(w_star, u_star, params) returns the matrices A and B of the
%   linearized system dw_e = A*w_e + B*u_e, linearized about the desired
%   trajectory w_star and u_star. The matrix forms are generated
%   symbolically by SYMSOLVEEQSMOTION and the copied here for numerical
%   evaluation.
% ----------------------------------------------------------------------
% INPUTS:
%   w_star - the current desired state
%   u_star - the current desired control input
%   params - a structure of problem parameters (masses, inertias, lengths,
%       etc.)
% OUTPUTS:
%   A - the matrix dfdw evaluated at (w_star,u_star)
%   B - the matrix dfdu evaluated at (w_star,u_star)

 x = w_star(1);  y = w_star(2);  th = w_star(3);  ph = w_star(4);
dx = w_star(5); dy = w_star(6); dth = w_star(7); dph = w_star(8);

% Parameters
m = params.m;
B = params.B;
C = params.C;
a = params.a;

% Additional force variables (currently all zero except Tapp=u)
[Fx,Fy,Fxw,Fyw,Tresc,Tresb] = forces(w_star,u_star,params);
k1    = params.k1;

%dfdw
A = [ 0, 0, 0, 0, 1, 0, 0, 0;
      0, 0, 0, 0, 0, 1, 0, 0;
      0, 0, 0, 0, 0, 0, 1, 0;
      0, 0, 0, 0, 0, 0, 0, 1;
      0, 0,  (a*(u_star*cos(th) - Tresc*cos(th) + Fyw*a*cos(2*th) - Fxw*a*sin(2*th)))/C - (Fx*sin(2*th) - Fyw*cos(2*th) - Fy*cos(2*th) + Fxw*sin(2*th) + dth*dx*m*cos(2*th) + dth*dy*m*sin(2*th))/m, 0, -(dth*sin(2*th))/2,    -dth*sin(th)^2, - (dx*sin(2*th))/2 - dy*sin(th)^2, 0;
      0, 0,  (Fx*cos(2*th) + Fxw*cos(2*th) + Fy*sin(2*th) + Fyw*sin(2*th) + dth*dy*m*cos(2*th) - dth*dx*m*sin(2*th))/m + (a*(u_star*sin(th) - Tresc*sin(th) + Fxw*a*cos(2*th) + Fyw*a*sin(2*th)))/C, 0,      dth*cos(th)^2, (dth*sin(2*th))/2,   (dy*sin(2*th))/2 + dx*cos(th)^2, 0;
      0, 0,  (a*(Fxw*cos(th) + Fyw*sin(th)))/C, 0, 0, 0, 0, 0;
      0, 0, -(a*(Fxw*cos(th) + Fyw*sin(th)))/C, 0, 0, 0, 0, 0];
 
 
%dfdu 
B = [ 0;
      0;
      0;
      0;
      (a*sin(th))/C;
     -(a*cos(th))/C;
     -1/C;
      1/B + 1/C];
  
end