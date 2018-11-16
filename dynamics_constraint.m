function h_i = dynamics_constraint(x_i, u_i, x_ip1, u_ip1, dt, params)
%DYNAMICS_CONSTRAINT(xm, um, xp, up, dt) computes the vector contstraint
%   h_i = \dot{s}_i(\Delta t/2) - f(s_i(\Delta t/2), 0.5(u_{i+1}+u{i})) = 0
%
%   @param x_i: x_i, the state at the start of the interval.
%   @param u_i: u_i, the input at the start of the interval.
%   @param x_ip1: x_{i+1}, the state at the end of the interval.
%   @param u_ip1: u_{i+1}, the input at the end of the interval.
%   @param dt: \Delta t, the duration of the interval
%
%   @outpu h_i: quantity from above expresion that should be 0

% MODIFIED - Credit to Mat Halm


f_i = f(x_i, u_i, params);
f_ip1 = f(x_ip1, u_ip1, params);

s_i_mid = 0.5 * (x_i + x_ip1) - (dt / 8) * (f_ip1 - f_i);

s_i_dot_mid = (1.5 / dt) * (x_ip1 - x_i) - 0.25 * (f_i + f_ip1);

%size(s_i_mid)
%size(s_i_dot_mid)


u_mid = 0.5*(u_i + u_ip1);

h_i = s_i_dot_mid - f(s_i_mid, u_mid, params);

%[x_i, u_i*ones(size(x_i)),x_ip1,u_ip1*ones(size(x_i)), h_i]

end