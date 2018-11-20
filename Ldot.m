function [dLdt] = Ldot(t, L, Q, R, params)
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

xd = x_d(t);
ud = u_d(t);

A1 = [0,0, -(1/constants.m)*cos(xd(3))*(ud(1)+ud(2));
    0, 0, -(1/constants.m)*sin(xd(3))*(ud(1)+ud(2));
    0, 0, 0];

A = [zeros(3,3), eye(3);
     A1, zeros(3,3)];
 
B1 = [-(1/constants.m)*sin(xd(3)), -(1/constants.m)*sin(xd(3));
       (1/constants.m)*cos(xd(3)),  (1/constants.m)*cos(xd(3));
       constants.a/constants.I, -constants.a/constants.I];
   
B = [zeros(3,2);B1];

dLdt = -1/2*Q*(L^(-1))' - A'*L + 1/2*L*L'*B*R^(-1)*B'*L;  

end