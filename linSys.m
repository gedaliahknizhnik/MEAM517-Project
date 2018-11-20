function [A,B] = linSys(xd, ud, params)

 x = xd(1);  y = xd(2);  th = xd(3);  ph = xd(4);
dx = xd(5); dy = xd(6); dth = xd(7); dph = xd(8);

% Parameters
m = params.m;
B = params.B;
C = params.C;
a = params.a;

% Additional force variables (currently all zero except Tapp=u)
Tresb = params.Tresb;
Tresc = params.Tresc;
Fx    = params.Fx; 
Fy    = params.Fy; 
Fxw   = params.Fxw; 
Fyw   = params.Fyw;
k1    = params.k1;

 F = -0.2*sqrt(dx^2+dy^2);
 Fx = F*cos(th);
 Fy = F*sin(th);
 Tresb = -0.5*dth;
 Tresc = -0.5*dth;

%dfdx
A = [ 0, 0, 0, 0, 1, 0, 0, 0;
      0, 0, 0, 0, 0, 1, 0, 0;
      0, 0, 0, 0, 0, 0, 1, 0;
      0, 0, 0, 0, 0, 0, 0, 1;
      0, 0,  (a*(ud*cos(th) - Tresc*cos(th) + Fyw*a*cos(2*th) - Fxw*a*sin(2*th)))/C - (Fx*sin(2*th) - Fyw*cos(2*th) - Fy*cos(2*th) + Fxw*sin(2*th) + dth*dx*m*cos(2*th) + dth*dy*m*sin(2*th))/m, 0, -(dth*sin(2*th))/2,    -dth*sin(th)^2, - (dx*sin(2*th))/2 - dy*sin(th)^2, 0;
      0, 0,  (Fx*cos(2*th) + Fxw*cos(2*th) + Fy*sin(2*th) + Fyw*sin(2*th) + dth*dy*m*cos(2*th) - dth*dx*m*sin(2*th))/m + (a*(ud*sin(th) - Tresc*sin(th) + Fxw*a*cos(2*th) + Fyw*a*sin(2*th)))/C, 0,      dth*cos(th)^2, (dth*sin(2*th))/2,   (dy*sin(2*th))/2 + dx*cos(th)^2, 0;
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