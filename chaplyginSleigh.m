function [dw] =  chaplyginSleigh(t,w,params)
% CHAPLYGINSLEIGH(t,w,params) is an ODE45 wrapper function for the dynamics
%   f, which takes in the time and state and outputs the derivatives.
%   Modified from the HW 5 code creat by Mathew Halm.
% ------------------------------------------------------------------------
% INPUTS:
%   t  - the current time
%   w  - the state vector (x,y,th,ph,dx,dy,dth,dph)
%   params - a structure of problem parameters (masses, inertias, lengths,
%   etc.)
%       Tapp = [t, T] - the torques returned by SNOPT
% OUTPUTS
%   dw - the derivative of the state vector (dx, dy, dth, dph, ddx, ddy,
%       ddth, ddph)

%% Select control input
TappVals = params.Tapp;

% Linearly interpolate between them
for ii = 1:size(TappVals,1)
    if (TappVals(ii,1) > t)
        break;
    end
end

dt = TappVals(ii,1)-TappVals(ii-1,1);
dTapp = TappVals(ii,2)-TappVals(ii-1,2);
Tapp = TappVals(ii-1,2) + (dTapp)*(t-TappVals(ii-1,1))/dt;

%% Calculate dynamics
dw = f(w,Tapp,params);

waitbar(t/params.maxTime,params.waitBar,'Solving ODE')

end