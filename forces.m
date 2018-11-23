function [Fx,Fy,Fxw,Fyw,Tresc,Tresb] = forces(w,u,params)

    % Extract state variables
     x = w(1);  y = w(2);  th = w(3);  ph = w(4);
    dx = w(5); dy = w(6); dth = w(7); dph = w(8);

    Fxw   = 0; 
    Fyw   = 0;

    F = 0;%-0.1*sqrt(dx^2+dy^2);
    Fx = F*cos(th);
    Fy = F*sin(th);
    
    % Cross-wind
%     Fy = -0.1;
%     Fx = 0.1;

    Tresb = 0;%-0.5*dph;   
    Tresc = 0;%-0.5*dth;

end