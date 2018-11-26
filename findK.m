function [K,u] = findK(t, w, ts, Ls, w_stars, u_stars, polys, params)
% FINDK(t,w, ts, Ls, w_stars, u_stars, polys, params) finds the controller
%   gain matrix K and the optimal control u = u_star + K(w_star - w) for
%   the TVLQR controller given by the ODE solution structure Ls.
% ------------------------------------------------------------------------
% INPUTS:
%   t       - the current time in the integration.
%   w       - the current state.
%   ts      - a vector of times corresponding to the SNOPT solutions w_stars
%       and u_stars.
%   w_stars - a matrix of the state solutions to the optimal trajectory
%       returned by SNOPT.
%   u_stars - a vector of the control solutions to the optimal trajectory
%       returned by SNOPT.
%   polys   - the polynomials representing the trajectory returned by SNOPT
%       (currently unused)
%   params  - a structure of problem parameters (masses, inertias, lengths,
%       etc.)
% OUTPUTS:
%   K - the matrix of controller gains
%   u - the optimal control input u = u_star + K(w_star - w)

    % Evaluate current optimal trajectory point.
    [w_star,u_star] = getDes(t,ts,w_stars,u_stars,polys);

    % Linearize system about desired trajectory
    [~,B] = linSys(w_star,u_star,params);
 
    % Evaluate L at the current time and then reshape it into a matrix
    % (passed as a vector used by the ODE structure).
    L = reshape(deval(Ls,t),params.nw,params.nw);
    % Calculate the cost-to-go.
    S = L*L';
    % Calculate the optimal controller gain
    K = params.R^(-1)*B'*S;
    % Calculate the optimal controller input.
    u = u_star + K*(w_star-w);
end