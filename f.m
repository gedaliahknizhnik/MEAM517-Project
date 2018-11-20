function f = f(w,u,params)

% Extract state variables
x = w(1);   y = w(2);  th = w(3);  ph = w(4);
dx = w(5); dy = w(6); dth = w(7); dph = w(8);

% Parameters
m  = params.m;
B  = params.B;
C  = params.C;
a  = params.a;
k1 = params.k1;

% Additional force variables (currently all zero except Tapp=u)
[Fx,Fy,Fxw,Fyw,Tresc,Tresb] = forces(w,u);

% Regular
% ddx = (Fx + Fxw + (Fy*sin(2*th))/2 + (Fyw*sin(2*th))/2 - Fx*sin(th)^2 - Fxw*sin(th)^2 - (dth*dx*m*sin(2*th))/2 - dth*dy*m*sin(th)^2)/m + ((Fyw*a^2*m*sin(2*th))/2 - Fxw*a^2*m*sin(th)^2 + u*a*m*sin(th) - Tresc*a*m*sin(th))/(C*m);
% ddy = (Fy + Fyw - Fy*cos(th)^2 - Fyw*cos(th)^2 + (Fx*sin(2*th))/2 + (Fxw*sin(2*th))/2 + dth*dx*m*cos(th)^2 + (dth*dy*m*sin(2*th))/2)/m - (Fyw*a^2*m*cos(th)^2 - (Fxw*a^2*m*sin(2*th))/2 + u*a*m*cos(th) - Tresc*a*m*cos(th))/(C*m);

% Adjusting for integration error
ddx = (dy*k1*sin(2*th))/2 - (Fxw*(sin(th)^2 - 1))/m - (dth*dx*sin(2*th))/2 - dth*dy*sin(th)^2 - (Fx*(sin(th)^2 - 1))/m - dx*k1*sin(th)^2 + (Fy*sin(2*th))/(2*m) + (Fyw*sin(2*th))/(2*m) + (Fyw*a^2*sin(2*th))/(2*C) - (Fxw*a^2*sin(th)^2)/C - a*dth*k1*sin(th) + (u*a*sin(th))/C - (Tresc*a*sin(th))/C;
ddy = (Fy + Fyw - Fy*cos(th)^2 - Fyw*cos(th)^2 + (Fx*sin(2*th))/2 + (Fxw*sin(2*th))/2 + dth*dx*m*cos(th)^2 + (dth*dy*m*sin(2*th))/2 - dy*k1*m*cos(th)^2  + (dx*k1*m*sin(2*th))/2 + a*dth*k1*m*cos(th))/m - (Fyw*a^2*m*cos(th)^2 - (Fxw*a^2*m*sin(2*th))/2 + u*a*m*cos(th) - Tresc*a*m*cos(th))/(C*m);

% Regardless
ddth = -(u - Tresc + Fyw*a*cos(th) - Fxw*a*sin(th))/C;
ddph = (u - Tresc + Fyw*a*cos(th) - Fxw*a*sin(th))/C + (C*u + C*Tresb)/(B*C);

f = [dx;dy;dth;dph;ddx;ddy;ddth;ddph];

%[w,u*ones(size(w)),f]
%pause

end