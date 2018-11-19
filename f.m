function f = f(w,u,params)

% Extract state variables
x = w(1);   y = w(2);  th = w(3);  ph = w(4);
dx = w(5); dy = w(6); dth = w(7); dph = w(8);

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

F = -3*sqrt(dx^2+dy^2);
Fx = F*cos(th);
Fy = F*sin(th);

ddx = (Fx + Fxw + (Fy*sin(2*th))/2 + (Fyw*sin(2*th))/2 - Fx*sin(th)^2 - Fxw*sin(th)^2 - (dth*dx*m*sin(2*th))/2 - dth*dy*m*sin(th)^2)/m + ((Fyw*a^2*m*sin(2*th))/2 - Fxw*a^2*m*sin(th)^2 + u*a*m*sin(th) - Tresc*a*m*sin(th))/(C*m);
ddy = (Fy + Fyw - Fy*cos(th)^2 - Fyw*cos(th)^2 + (Fx*sin(2*th))/2 + (Fxw*sin(2*th))/2 + dth*dx*m*cos(th)^2 + (dth*dy*m*sin(2*th))/2)/m - (Fyw*a^2*m*cos(th)^2 - (Fxw*a^2*m*sin(2*th))/2 + u*a*m*cos(th) - Tresc*a*m*cos(th))/(C*m);
ddth = -(u - Tresc + Fyw*a*cos(th) - Fxw*a*sin(th))/C;
ddph = (u - Tresc + Fyw*a*cos(th) - Fxw*a*sin(th))/C + (C*u + C*Tresb)/(B*C);

f = [dx;dy;dth;dph;ddx;ddy;ddth;ddph];

%[w,u*ones(size(w)),f]
%pause

end