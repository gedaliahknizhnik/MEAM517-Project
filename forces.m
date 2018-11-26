function [Fx,Fy,Fxw,Fyw,Tresc,Tresb] = forces(w,u,params)
% FORCES(w,u,params) returns the external forces acting on the Chaplygin
%   Beanie, including dissipative/frictional forces.
% ----------------------------------------------------------------------
% INPUTS:
%   w      - the current state
%   u      - the current control input
%   params - a structure of problem parameters (masses, inertias, lengths,
%       etc.)
% OUTPUTS:
%   Fx    - the x-axis force acting on the COM
%   Fy    - the y-axis force acting on the COM
%   Fxw   - the x-axis force acting on the wheel
%   Fyw   - the y-axis force acting on the wheel
%   Tresc - the resistive torque acting on C only.
%   Tresb - the resistive torque acting on B only.

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