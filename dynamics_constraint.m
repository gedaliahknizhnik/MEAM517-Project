function h_i = dynamics_constraint(w_i, u_i, w_ip1, u_ip1, params)
% DYNAMICS_CONSTRAINT(w_i, u_i, w_{i+1}, u_{i+1}, params) computes the
%   vector collocation constraint on the midpoint of a cubic spline:
%       h_i = \dot{s}_i(\Delta t/2) - f(s_i(\Delta t/2), 0.5(u_{i+1}+u{i})) = 0
% ------------------------------------------------------------------------
% INPUTS:
%   w_i    - w_i, the state at the start of the interval.
%   u_i    - u_i, the input at the start of the interval.
%   w_ip1  - w_{i+1}, the state at the end of the interval.
%   u_ip1  - u_{i+1}, the input at the end of the interval.
%   params - a structure of problem parameters (masses, inertias, lengths,
%       etc.)
% OUTPUTS
%   h_i - quantity from above expresion that should be 0
%
% MODIFIED from MEAM517 HW 5 solution code. Credit to Mat Halm


f_i = f(w_i, u_i, params);
f_ip1 = f(w_ip1, u_ip1, params);

s_i_mid = 0.5 * (w_i + w_ip1) - (params.dt / 8) * (f_ip1 - f_i);

s_i_dot_mid = (1.5 / params.dt) * (w_ip1 - w_i) - 0.25 * (f_i + f_ip1);

u_mid = 0.5*(u_i + u_ip1);

h_i = s_i_dot_mid - f(s_i_mid, u_mid, params);

end