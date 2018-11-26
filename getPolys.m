function [polys] = getPolys(w_stars,u_stars,params)
% GETPOLYS(w_stars,u_stars,params) computes the polynomial trajectories
%   that span the SNOPT outputs w_stars. This corresponds to the
%   collocation splines used in dynamics_contraint().
% ---------------------------------------------------------------------
% INPUTS: 
%   w_stars - a matrix of SNOPT output states
%   u_stars - a vector of SNOPT output controls
%   params  - a structure of problem parameters (masses, inertias, lengths,
%       etc.) 
% OUTPUTS:
%   polys - a matrix of polynomial coefficients (in descending order of
%       power)

    nw = params.nw; dt = params.dt;
    polys = zeros(size(w_stars,1)-1,4,nw);

    for ii=1:(size(w_stars,1)-1)     
        % Extract points
        w_i   = w_stars(ii,:);
        w_ip1 = w_stars(ii+1,:);
        f_i   = f(w_stars(ii,:),u_stars(ii,:),params)';
        f_ip1 = f(w_stars(ii+1,:),u_stars(ii+1,:),params)';
        
        % Create the polynomials (in descending order of power)
        polys(ii,4,:) = w_i;
        polys(ii,3,:) = f_i;
        polys(ii,2,:) = 3/dt^2*(w_ip1-w_i) - 1/dt*(f_ip1 + 2*f_i);
        polys(ii,1,:) = -2/dt^3*(w_ip1-w_i) + 1/dt^2*(f_i + f_ip1);
    end
    
    
end