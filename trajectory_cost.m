function [g,dG] = trajectory_cost(z, N, nx, nu, dt)
%TRAJECTORY_COST(z) computes the cost and cost jacobian.
%   @param z: vector of decision variables containing the x_i and u_i.
%   @param N: number of sample points
%   @param nx: dimension of state vector, x
%   @param nu: dimension of input vector, u
%   @param dt: \Delta t, the inter-sample interval duration

%   @output g: total accrued cost.
%   @output dG_i: jacobian of total  accrued cost.
%
% MODIFIED - Credit to Mat Halm

g = 0;
dG = zeros(N*(nx + nu),1);

for i=1:(N-1)
   u_i_inds = (1:nu) + nx * i + nu * (i - 1);
   u_ip1_inds = (1:nu) + nx * (i+1) + nu * i;
   g = g + dt*0.5*norm(z(u_i_inds))^2 + dt*0.5*norm(z(u_ip1_inds))^2;
   dG(u_i_inds) = dG(u_i_inds) + dt*z(u_i_inds);
   dG(u_ip1_inds) = dG(u_ip1_inds) + dt*z(u_ip1_inds);
    
end


end