function [g,dG] = trajectory_cost(z, params)
% TRAJECTORY_COST(z, params) computes the cost and cost jacobian. The cost
%   is assumed to be the integral of the control effort over the time
%   interval, approximated by a trapezoidal sum.
% -------------------------------------------------------------------
% INPUTS:
%   z      - vector of decision variables containing the w_i and u_i.
%   params - a structure of problem parameters (masses, inertias, lengths,
%       etc.)
% OUTPUTS:
%   g    - total accrued cost.
%   dG_i - jacobian of total  accrued cost.
%
% MODIFIED from MEAM517 HW 5 solution code. Credit to Mat Halm

nw = params.nw; nu = params.nu; N = params.N; dt = params.dt;

g = 0;
dG = zeros(N*(nw + nu),1);

% Calculate the cost of the trajectory
for i=1:(N-1)
   u_i_inds = (1:nu) + nw * i + nu * (i - 1);
   u_ip1_inds = (1:nu) + nw * (i+1) + nu * i;
   g = g + dt*0.5*norm(z(u_i_inds))^2 + dt*0.5*norm(z(u_ip1_inds))^2;
   dG(u_i_inds) = dG(u_i_inds) + dt*z(u_i_inds);
   dG(u_ip1_inds) = dG(u_ip1_inds) + dt*z(u_ip1_inds);
    
end


end