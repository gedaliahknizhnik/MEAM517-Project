function [dLdt] = Ldot(t, L, Q, R, td, xd, ud, polys, params)
% LDOT(t, L, Q, R, constants) calculates the time derivative of L(t).
%    Input parameters:
%    t - current time
%    L - current value of L(t)
%    Q - state cost matrix
%    R - input cost matrix
%    constants - struct containing inertial and geometric constants for the
%        quadrotor system, contained in constants.g, contstants.I, etc.

% modify this line to calculate \dot L(t).
%dLdt = zeros(6,6);

L = reshape(L,params.nx,params.nx);


% xd = x_d(t);
% ud = u_d(t);

[xdNow,udNow] = getDes(t,td,xd,ud,polys);

[A,B] = linSys(xdNow,udNow,params);

dLdt = -1/2*Q*(L^(-1))' - A'*L + 1/2*L*L'*B*R^(-1)*B'*L;

dLdt = reshape(dLdt,params.nx^2,1);

end