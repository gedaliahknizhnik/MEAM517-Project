function [h_i,dH_i] = dynamics_constraint_with_derivative(x_i, u_i, x_ip1, u_ip1, dt, params)
% DYNAMICS_CONSTRAINT_WITH_DERIVATIVE(x_i, u_i, x_ip1, u_ip1, dt) returns
%   and computes the gradient of the vector constraint asssociated with
%   dynamics_constraint(x_i, u_i, x_ip1, u_ip1, dt).
%
%   @param x_i: x_i, the state at the start of the interval.
%   @param u_i: u_i, the input at the start of the interval.
%   @param x_ip1: x_{i+1}, the state at the end of the interval.
%   @param u_ip1: u_{i+1}, the input at the end of the interval.
%   @param dt: \Delta t, the duration of the interval
%
%   @output h_i: constraint value from dynamics_constraint
%   @output dH_i: jacobian of h_i evaluated at(x_i, u_i, x_ip1, u_ip1).

  h_i = dynamics_constraint(x_i, u_i, x_ip1, u_ip1, dt, params);

  % use numerical derivatives to compute dH
  % dH = [dh/dx0 dh/du0 dh/dx1 dh/du1]
  % where the partial derivatives are written (dh/dx0)_ij = dh_i/dx0_j
  delta = 1e-8;
  dH_i = zeros(numel(x_i), 2*(numel(x_i)+numel(u_i)));
  for j=1:numel(x_i)
      dx = zeros(size(x_i));
      %dx = zeros(numel(x_i),1);
      dx(j) = delta;
      dHx_i_j = dynamics_constraint(x_i + dx, u_i, x_ip1, u_ip1, dt, params) - h_i;
      dHx_ip1_j = dynamics_constraint(x_i, u_i, x_ip1 + dx, u_ip1, dt, params) - h_i;
      dH_i(:,j) = dHx_i_j/delta;
      dH_i(:,j + numel(x_i) + numel(u_i)) = dHx_ip1_j/delta;
  end
  
  for j=1:numel(u_i)
      du = zeros(numel(u_i),1);
      du(j) = delta;
      dHu_i_j = dynamics_constraint(x_i, u_i + du, x_ip1, u_ip1, dt, params) - h_i;
      dHu_ip1_j = dynamics_constraint(x_i, u_i, x_ip1, u_ip1 + du, dt, params) - h_i;
      dH_i(:,j + numel(x_i)) = dHu_i_j/delta;
      dH_i(:,j + numel(x_i) + numel(u_i) + numel(x_ip1)) = dHu_ip1_j/delta;
  end

end