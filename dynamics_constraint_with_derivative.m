function [h_i,dH_i] = dynamics_constraint_with_derivative(w_i, u_i, w_ip1, u_ip1, params)
% DYNAMICS_CONSTRAINT_WITH_DERIVATIVE(w_i, u_i, w_ip1, u_ip1, params) returns
%   and computes the gradient of the vector constraint asssociated with
%   dynamics_constraint(w_i, u_i, w_ip1, u_ip1, params).
% ------------------------------------------------------------------------
% INPUTS:
%   w_i    - w_i, the state at the start of the interval.
%   u_i    - u_i, the input at the start of the interval.
%   w_ip1  - w_{i+1}, the state at the end of the interval.
%   u_ip1  - u_{i+1}, the input at the end of the interval.
%   params - a structure of problem parameters (masses, inertias, lengths,
%     etc.)
% OUTPUTS:
%   h_i: constraint value from dynamics_constraint
%   dH_i: jacobian of h_i evaluated at(w_i, u_i, w_ip1, u_ip1).
%
% MODIFIED from MEAM517 HW 5 solution code. Credit to Mat Halm

  h_i = dynamics_constraint(w_i, u_i, w_ip1, u_ip1, params);

  % use numerical derivatives to compute dH
  % dH = [dh/dx0 dh/du0 dh/dx1 dh/du1]
  % where the partial derivatives are written (dh/dx0)_ij = dh_i/dx0_j
  delta = 1e-8;
  dH_i = zeros(numel(w_i), 2*(numel(w_i)+numel(u_i)));
  for j=1:numel(w_i)
      dx = zeros(size(w_i));
      %dx = zeros(numel(x_i),1);
      dx(j) = delta;
      dHx_i_j = dynamics_constraint(w_i + dx, u_i, w_ip1, u_ip1, params) - h_i;
      dHx_ip1_j = dynamics_constraint(w_i, u_i, w_ip1 + dx, u_ip1, params) - h_i;
      dH_i(:,j) = dHx_i_j/delta;
      dH_i(:,j + numel(w_i) + numel(u_i)) = dHx_ip1_j/delta;
  end
  
  for j=1:numel(u_i)
      du = zeros(numel(u_i),1);
      du(j) = delta;
      dHu_i_j = dynamics_constraint(w_i, u_i + du, w_ip1, u_ip1, params) - h_i;
      dHu_ip1_j = dynamics_constraint(w_i, u_i, w_ip1, u_ip1 + du, params) - h_i;
      dH_i(:,j + numel(w_i)) = dHu_i_j/delta;
      dH_i(:,j + numel(w_i) + numel(u_i) + numel(w_ip1)) = dHu_ip1_j/delta;
  end

end