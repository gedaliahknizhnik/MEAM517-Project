function [ceq, dCeq] = all_constraints(z, params)
% ALL_CONSTRAINTS(z, N, nx, nu, dt, params) returns all the vector
%   contstraints on the system and their gradients via finite differencing.
%   Modified from the HW 5 code creat by Mathew Halm.
% ------------------------------------------------------------------------
% INPUTS:
%   z  - the final state vector (including all xs and us)    
%   nx - the size of the state vector x
%   nu - the size of the control vector u
%   N  - the number of segments
%   dt - the time between segments (constant)
%   params - a structure of problem parameters (masses, inertias, lengths,
%       etc.)
% OUTPUTS
%   ceq  - a vector containing all the collocation constraint values (should
%           equal 0)
%   dCeq - a sparse matrix containing the gradient of the constraint vector
%           ceq with respect to the state vector z. Returned as a sparse
%           matrix since each constraint depends only on x_i, x_(i+1), u_i,
%           and u_(i+1).
%
% MODIFIED from MEAM517 HW 5 solution code. Credit to Mat Halm

[ceq, dCeq] = dynamics_constraints(z, params);

dCeq = sparse(dCeq);

end