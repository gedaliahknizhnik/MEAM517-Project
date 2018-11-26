function [dLdt] = Ldot(t, L, t_stars, w_stars, u_stars, polys, params)
% LDOT(t, L, Q, R, constants) calculates the time derivative of L(t). This
%   is used to solve the time-varying Riccati equation and derive the LQR
%   controller gains. This function is designed to be called by an ODE
%   solver and so treats L as a vector.
% -------------------------------------------------------------------------
% INPUTS:
%   t       - current time
%   L       - current value of L(t) passed as a column vector (for ODE45)
%   t_stars - a vector of times corresponding to the SNOPT outputs w_stars
%       and u_stars
%   w_stars - a matrix of SNOPT state outputs
%   u_stars - a vector of SNOPT control outputs
%   polys   - the polynomial trajectory output by SNOPT (currently unused)
%   params  - a structure of problem parameters (masses, inertias, lengths,
%       etc.) as well as R and Q
% OUTPUTS:
%   dLdt - the derivative part of the Riccati equation, reshaped into a
%       column vector.
%
% MODIFIED from MEAM517 HW 5 solution code. Credit to Mat Halm

Q = params.Q;
R = params.R;

L = reshape(L, params.nw, params.nw);

[xdNow,udNow] = getDes(t,t_stars,w_stars,u_stars,polys);

[A,B] = linSys(xdNow,udNow,params);

dLdt = -1/2*Q*(L^(-1))' - A'*L + 1/2*L*L'*B*R^(-1)*B'*L;

dLdt = reshape(dLdt,params.nw^2,1);

%frac = (params.maxTime - t)/params.maxTime;
%waitbar(frac,params.waitBar,'Solving Riccati equation ODE')

end