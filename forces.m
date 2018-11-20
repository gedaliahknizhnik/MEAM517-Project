function [Fx,Fy,Fxw,Fyw,Tresc,Tresb] = forces(w,u)

    % Extract state variables
     x = w(1);  y = w(2);  th = w(3);  ph = w(4);
    dx = w(5); dy = w(6); dth = w(7); dph = w(8);

    Fxw   = 0; 
    Fyw   = 0;

    F = -0.2*sqrt(dx^2+dy^2);
    Fx = F*cos(th);
    Fy = F*sin(th);

    Tresb = -0.5*dth;
    Tresc = -0.5*dth;

end